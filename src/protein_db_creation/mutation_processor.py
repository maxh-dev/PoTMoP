import json
from tqdm import tqdm
from pyensembl import EnsemblRelease

from src.logger import logger
from src.protein_db_creation.mutation_processing.maf_parser import MAFParser
from src.protein_db_creation.mutation_processing.vcf_parser import VCFParser
from src.protein_db_creation.mutation_processing.mutation_protein_converter import MutationToProteinConverter
from src.util.fasta_handler import FastaReader
from src.util.general import load_or_create_json, remove_duplicate_sequences

from src.util.protein_nameing import create_fasta_header

class MutationProcessor:
    def __init__(self, release_number: str, mutation_info_json: str, mutation_directory: str, fasta_file: str, protein_reference_json: str):
        """
        Initialize the MutationProcessor with required parameters for release number, mutation directory,
        fasta file, and protein reference JSON.
        
        :param release_number: Release version for Ensembl data (e.g., "104").
        :param mutation_info_json: 
        :param mutation_directory: Directory path containing mutation files (VCF, MAF).
        :param fasta_file: Path to the FASTA file for storing mutation sequences.
        :param protein_reference_json: Path to the JSON file containing protein reference data.
        """
        self.release_number = release_number

        self.mutation_info_json = mutation_info_json
        with open(self.mutation_info_json, 'r') as f:
            self.mutation_info_dict = json.load(f)

        self.mutation_directory = mutation_directory
        self.fasta_file = fasta_file

        self.protein_reference_json = protein_reference_json
        self.protein_reference_dict = load_or_create_json(protein_reference_json)

        self.ensembl_data = EnsemblRelease(self.release_number)

        self.mutations_dict = {}
        self.mutations_list = []

        maf_parser = MAFParser(self.mutation_directory)
        self.mutations_dict.update(maf_parser.get_variants())
        
        #TODO: finish when needed
        # vcf_parser = VCFParser(self.mutation_directory)
        # self.mutations_dict.update(vcf_parser.get_variants())

    def append_entry(self, header: str, sequence: str) -> None:
        """
        Appends a new entry to the FASTA file with the given header and sequence.
        
        :param header: Header for the FASTA entry.
        :param sequence: Sequence data to be written under the header.
        """
        entry_string = ""
        entry_string += f"\n{header}"
        entry_string += f"{sequence[:-1]}"

        with open(self.fasta_file, 'a') as f:
            f.write(entry_string)

    def append_mutations_list_to_fasta(self) -> None:
        """
        Appends all entries from mutations_list to the FASTA file.
        Each entry contains a header and a protein sequence.
        """
        fr = FastaReader(self.fasta_file)
        fr.read_fasta()
        sequence_to_header_dict = fr.get_sequence_to_header_dict()

        for elem in self.mutations_list:
            transcript_id, mutation_str, mutated_protein = elem
            mutation_id = str(self.mutation_info_dict["mutation_type_to_id"][mutation_str])
            mutation_dict = self.mutation_info_dict["mutation_information"][mutation_id]
            if not mutation_id in mutation_dict["mutation_db_id"]:
                continue
            if mutation_dict["sequence"] in sequence_to_header_dict:
                logger.info(f"Skipping {mutation_str} already present.")
                continue
            else:
                sequence_to_header_dict[mutation_dict["sequence"]] = None
            gene_name = self.ensembl_data.transcript_by_id(mutation_dict["transcript_id"]).gene_name
            self.append_entry(create_fasta_header(mutation_dict["mutation_db_id"], gene_name), mutation_dict["sequence"])

            self.protein_reference_dict[mutation_dict["mutation_db_id"]] = {
                "uniprot_id": None,
                "mutation": mutation_str
            }

        logger.info(f"Saving protein_reference to {self.protein_reference_json}")
        with open(self.protein_reference_json , "w") as file:
            json.dump(self.protein_reference_dict, file, indent=4)

    def insert_new_mutation(self, transcript_id, mutation_str,mutated_protein):
        last_mutation_id = self.mutation_info_dict["mutation_count"]
        new_mutation_id = last_mutation_id + 1
        self.mutation_info_dict["mutation_count"] = new_mutation_id
        self.mutation_info_dict["sequences_to_id"][mutated_protein] = [new_mutation_id]
        self.mutation_info_dict["mutation_type_to_id"][mutation_str] = new_mutation_id
        self.mutation_info_dict["mutation_information"][str(new_mutation_id)] = {
            "transcript_id": transcript_id,
            "mutation_db_id": f"{transcript_id}_MUT{new_mutation_id}",
            "mutation_str": mutation_str,
            "sequence": mutated_protein
        }
        
    def insert_existing_mutation(self, transcript_id, mutation_str, mutated_protein, existing_mutation_db_id):
        last_mutation_id = self.mutation_info_dict["mutation_count"]
        new_mutation_id = last_mutation_id + 1
        self.mutation_info_dict["mutation_count"] = new_mutation_id

        self.mutation_info_dict["sequences_to_id"][mutated_protein].append(new_mutation_id)
        self.mutation_info_dict["mutation_type_to_id"][mutation_str] = new_mutation_id
        self.mutation_info_dict["mutation_information"][str(new_mutation_id)] = {
            "transcript_id": transcript_id,
            "mutation_db_id": existing_mutation_db_id,
            "mutation_str": mutation_str,
            "sequence": mutated_protein
        }
    

    def insert_new_mutations(self):
        for elem in self.mutations_list:
            transcript_id, mutation_str, mutated_protein = elem

            if mutated_protein in self.mutation_info_dict["sequences_to_id"]:
                mutation_str_present = False 
                for found_mutation_id in self.mutation_info_dict["sequences_to_id"][mutated_protein]:
                    found_mutation_id = str(found_mutation_id)
                    found_mutation_entry = self.mutation_info_dict["mutation_information"][found_mutation_id]
                    mutation_db_id = self.mutation_info_dict["mutation_information"][found_mutation_id]["mutation_db_id"]
                    if mutation_str == found_mutation_entry["mutation_str"]:
                        mutation_str_present = True
                        continue
                if mutation_str_present:
                    continue
                #if we arrive at this point, the mutation is not present yet 
                self.insert_existing_mutation(transcript_id, mutation_str, mutated_protein, mutation_db_id)
            else: 
                self.insert_new_mutation(transcript_id, mutation_str, mutated_protein)

        with open(self.mutation_info_json , "w") as file:
            json.dump(self.mutation_info_dict, file, indent=4)


    def translate_mutation(self, chromosome: str, mutation_position_genomic: str, ref: str, alt: str):
        """
        Process a single mutation by creating an instance of MutationToProteinConverter, 
        computing mutated sequences, and storing the results.
        
        :param chromosome: Chromosome name or number (e.g., "chr1").
        :param mutation_position_genomic: Genomic position of the mutation (e.g., "123456").
        :param ref: Reference allele for the mutation.
        :param alt: Alternate allele for the mutation.
        """
        mutation_instance = MutationToProteinConverter(
            ensembl_data = self.ensembl_data,
            chromosome = chromosome,
            mutation_position = int(mutation_position_genomic),
            ref_allele = ref,
            alt_allele = alt
        )
        
        mutation_instance.compute_mutated_sequences()
        result_list = mutation_instance.get_mutated_proteins()
        self.mutations_list.extend(result_list)

    def translate_mutations(self) -> None:
        """
        Process all mutations in the mutations_dict by creating instances of MutationToProteinConverter,
        calling compute_mutated_sequences(), and converting the mutations to proteins.
        
        For each mutation, the REF and ALT alleles are checked for SNP (single nucleotide polymorphism).
        Only valid SNPs are processed.
        """
        i = 0
        for mutation_key, mutation_data in tqdm(self.mutations_dict.items(), desc="Processing Mutations"):
            if mutation_data["REF"] != mutation_data["ALT_1"]:
                self.translate_mutation(mutation_data["Chromosome"], mutation_data["PositionGenomic"], mutation_data["REF"], mutation_data["ALT_1"])
            if mutation_data["REF"] != mutation_data["ALT_2"]:
                self.translate_mutation(mutation_data["Chromosome"], mutation_data["PositionGenomic"], mutation_data["REF"], mutation_data["ALT_2"])
        