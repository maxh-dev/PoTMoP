from pyensembl import EnsemblRelease
from Bio.Seq import Seq

from src.logger import logger


class MutationToProteinConverter:
    def __init__(self, ensembl_data: EnsemblRelease, chromosome: str, mutation_position: int, ref_allele: str, alt_allele: str) -> None:
        """
        Initialize MutationToProteinConverter with the required data for mutation processing.
        
        :param ensembl_data: EnsemblRelease object containing gene data.
        :param chromosome: Chromosome name or number where mutation occurs (e.g., "chr1").
        :param mutation_position: Genomic position of the mutation (e.g., 123456).
        :param ref_allele: Reference allele for the mutation.
        :param alt_allele: Alternate allele for the mutation.
        """
        self.ensembl_data = ensembl_data
        self.chromosome = chromosome
        self.mutation_position_genomic = mutation_position
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele


        self.matches = []
        self.mutated_sequences = []

    def find_matching_transcripts(self) -> None:
        """
        Find all matching transcripts at the specified mutation position and chromosome.
        Matches are stored in self.matches
        """
        genes = self.ensembl_data.genes_at_locus(contig=self.chromosome, position=self.mutation_position_genomic)

        if not genes:
            return

        for gene in genes:
            transcripts = gene.transcripts
            if not transcripts:
                continue

            for transcript in transcripts:
                if transcript.start <= self.mutation_position_genomic <= transcript.end:
                    try:
                        # Convert from an absolute chromosomal position to the offset into this transcripts spliced mRNA.
                        mutation_position_transcript = transcript.spliced_offset(self.mutation_position_genomic)

                        if not transcript.start_codon_complete:
                            logger.debug(f"For {transcript.transcript_id} no start codon present.")
                            continue

                        if not transcript.start_codon_spliced_offsets[2] < mutation_position_transcript < transcript.stop_codon_spliced_offsets[0]:
                            logger.debug(f"For {transcript.transcript_id} mutation_position not inside coding region")
                            continue

                        transcript_sequence = transcript.coding_sequence
                        if transcript_sequence[mutation_position_transcript - 1] != self.ref_allele:
                            logger.debug(f"Reference allele does not match for transcript {transcript.transcript_id}. Skipping.")
                            continue

                        self.matches.append((transcript, mutation_position_transcript))
                    except:
                        logger.debug("absolute chromosomal position is outside any exons. Skipping.")
                        continue

    def get_coding_sequences(self) -> None:
        """
        Generate mutated sequences for all matching transcripts and store them in self.mutated_sequences.
        """
        for transcript, mutation_position_transcript in self.matches:
            transcript_sequence = transcript.coding_sequence
            mutated_sequence = (
                    transcript_sequence[:mutation_position_transcript - 1] +
                    self.alt_allele +
                    transcript_sequence[mutation_position_transcript:]
            )
            self.mutated_sequences.append((transcript, mutated_sequence, mutation_position_transcript))

    def compute_mutated_sequences(self) -> None:
        """
        Find matching transcripts and compute the mutated sequences for the mutation.
        """
        self.find_matching_transcripts()
        if not self.matches:
            logger.debug("No matches found.")
            return
        self.get_coding_sequences()
        logger.debug(f"For {self.chromosome}:{self.mutation_position_genomic} {self.ref_allele}/{self.alt_allele}: found {len(self.mutated_sequences)} mutated sequences.")

    def get_mutated_proteins(self) -> tuple[list, int]:
        """
        Generate the mutated proteins from the mutated sequences and return them along with the updated mutation counter.
        
        :return: A list of mutated protein sequences and the updated mutations counter.
        """
        result_list = []
        if not self.mutated_sequences:
            logger.debug("No proteins to create.")
            return result_list
        for transcript, mutated_sequence, mutation_position_transcript in self.mutated_sequences:
            mutated_protein = str(Seq(mutated_sequence).translate())
            # if mutation indroduces stop codon, resize the protein
            mutated_protein = f"{mutated_protein.split('*')[0]}*"

            org_protein = str(Seq(transcript.coding_sequence).translate())

            if mutated_protein == org_protein:
                continue 

            protein_mutation_position = None 
            for i in range(len(mutated_protein)):
                if mutated_protein[i] != org_protein[i]:
                    protein_mutation_position = i 
                    protein_aa_ref = org_protein[i]
                    protein_aa_alt = mutated_protein[i]
                    break
            if not protein_mutation_position:
                # Some transcript coding sequences have a prepature stop codon, this can cause this condition to be true
                logger.info(f"Mutation: {transcript.transcript_id}_{mutation_position_transcript}_{self.ref_allele}/{self.alt_allele} could not be processed.")
                continue
            if protein_aa_ref == protein_aa_alt:
                logger.info(f"Mutation: {transcript.transcript_id}_{mutation_position_transcript}_{self.ref_allele}/{self.alt_allele} did not affect the protein.")
                continue
            mutation = f"{self.chromosome}_{self.mutation_position_genomic}|{transcript.transcript_id}_{mutation_position_transcript}_{self.ref_allele}/{self.alt_allele}|{protein_mutation_position}_{protein_aa_ref}/{protein_aa_alt}"
            """header = create_fasta_header(
                gene_id = transcript.gene_id,
                gene_name = transcript.gene_name,
                protein_id = transcript.protein_id,
                transcript_id = id_of_mutation,
                transcript_name = transcript.transcript_name, 
                mutation = mutation
            )"""
            result_list.append((transcript.transcript_id, mutation, mutated_protein))
        return result_list