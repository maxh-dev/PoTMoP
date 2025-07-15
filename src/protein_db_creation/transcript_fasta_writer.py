from tqdm.auto import tqdm

from pyensembl import EnsemblRelease
from Bio.Seq import Seq
import re 

from src.logger import logger
from src.util.protein_nameing import create_fasta_header
from src.util.general import load_or_create_json
import json

class TranscriptFastaWriter:
    def __init__(self, release_number: int, sequence_to_uniprot: str, output_file: str, protein_reference_json: str):
        """
        Initialize the TranscriptFastaWriter class.

        Parameters:
        - release_number (int): Ensembl release number.
        - sequence_to_uniprot (str): Path to JSON mapping sequences to Uniprot.
        - output_file (str): Path to the output FASTA file.
        - protein_reference_json (str): Path to the protein reference JSON file.
        """
        self.release_number = release_number

        with open(sequence_to_uniprot) as f:
            self.sequence_to_uniprot = json.load(f)
        
        self.output_file = output_file
        self.protein_reference_json = protein_reference_json
        self.protein_reference_dict = load_or_create_json(protein_reference_json)
        self.ensembl_data = EnsemblRelease(release_number)

    def is_valid_aa_sequence(self, sequence: str) -> bool:
        """
        Check if the given amino acid sequence is valid.

        A valid sequence:
        - Ends with '*'
        - Contains only standard amino acid characters (ACDEFGHIKLMNPQRSTVWY, case insensitive, excluding '*')

        Parameters:
        - sequence (str): Amino acid sequence to validate.

        Returns:
        - bool: True if valid, False otherwise.
        """    
        if not sequence.endswith("*"):
            return False
        sequence = sequence[:-1]
        return bool(re.fullmatch(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence, re.IGNORECASE))

    def write_transcripts_to_fasta(self) -> None:
        """
        Iterates over transcripts and writes valid sequences to a FASTA file.
        Keeps track of duplicates, invalid sequences, and Uniprot mappings.
        """
        logger.info(f"Writing transcripts to {self.output_file}...")
        output_count = 0
        seen_sequences = {}
        duplicate_count = 0
        invalid_count = 0
        uni_count = 0
        transcript_ids = {}
        with open(self.output_file, 'w') as fasta_file:
            for transcript in tqdm(
                self.ensembl_data.transcripts(), desc="Parsing transcripts to memory..."
            ):
                transcript_id = transcript.transcript_id
                transcript_name = transcript.transcript_name
                gene_id = transcript.gene_id
                gene_name = transcript.gene_name
                protein_id = transcript.protein_id

                if not transcript.is_protein_coding:
                    continue

                if transcript.biotype != "protein_coding":
                    continue

                if not transcript.complete:
                    continue
   
                bps_sequence = transcript.coding_sequence
                        
                if not bps_sequence:
                    continue

                aa_sequence = str(Seq(bps_sequence).translate())

                if not self.is_valid_aa_sequence(aa_sequence):
                    #logger.info(f"Transcript {transcript_id} has invalid AA sequence.")
                    invalid_count += 1
                    continue

                aa_sequence = aa_sequence[:-1]

                if not gene_name:
                    gene_name = f"noname{gene_id}"

                # Skip duplicate sequences
                if aa_sequence in seen_sequences:
                    duplicate_count += 1
                    continue
                seen_sequences[aa_sequence] = True

                output_count += 1

                uniprot_id = None
                uniprot_name = None
                if aa_sequence in self.sequence_to_uniprot:
                    uniprot_id = self.sequence_to_uniprot[aa_sequence]["id"]
                    uniprot_name = self.sequence_to_uniprot[aa_sequence]["name"]
                    uni_count += 1

                # Write in FASTA format
                fasta_file.write(create_fasta_header(transcript_id, gene_name))
                fasta_file.write(f"{aa_sequence}\n")

                self.protein_reference_dict[transcript_id] = {
                    "uniprot_id": uniprot_id,
                    "mutation": None
                }

        logger.info(f"{output_count} transcripts written to {self.output_file}. {duplicate_count} were duplicates. {invalid_count} were invalid. {uni_count} could be mapped to Uniprot.")
        self._save_protein_reference()

    def _save_protein_reference(self) -> None:
        """
        Saves the protein reference dictionary to the JSON file.
        """
        logger.info(f"Saving protein reference to {self.protein_reference_json}")
        with open(self.protein_reference_json, "w") as file:
            json.dump(self.protein_reference_dict, file, indent=4)