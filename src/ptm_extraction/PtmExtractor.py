import pandas as pd
import re
from tqdm import tqdm
import json

from src.logger import logger

from src.util.fasta_handler import FastaReader
from src.util.general import get_transcript_id_from_header

ms_fragger_mods = {"42.0106": "Acetyl",
                   "79.9663": "Phospho",
                   "114.0429": "GG",
                   "14.0157": "Methyl",
                   "0.9840": "Citrullination",
                   "57.0215": "Carbamidomethylation"}

modification_sites = {'Phospho': ['S', 'T', 'Y'],
                          'Acetyl': ['K'],
                          'Methyl': ['K', 'R'],
                          'GG': ['K'],
                          'Citrullination': ['R'],
                          'Deamidated': ['N', 'Q', 'R'],}

class PtmExtractor:
    def __init__(self, experiment_file: str, selected_ptms: list, experiment_groups_file: str, protein_fasta: str, output_path: str):
        """
        Initializes the PtmExtractor class.
        
        :param experiment_file: Path to the file containing experimental data ( only work with combined_modified_peptide.tsv).
        :param selected_ptms: List of selected post-translational modifications that shpuld be extracted.
        :param experiment_groups_file: Path to the file containing the matching between sample name and group. 
        :param protein_fasta: Path to the protein FASTA file.
        :param output_path: Path where output JSON will be saved.
        """
        self.experiment_file = experiment_file
        self.selected_ptms = selected_ptms
        self.experiment_groups_file = experiment_groups_file
        self.groups_df = pd.read_csv(self.experiment_groups_file)

        self.protein_fasta = protein_fasta
        self.output_path = output_path

        self.transcript_mods = {}

    def obtain_proteins(self) -> None:
        """
        Reads protein data from the FASTA file and processes it to get transcript_mods and transcript_sequence.
        """
        fr = FastaReader(self.protein_fasta)
        fr.read_fasta()
        # maps transcript_id to a dict containing all samples which map to empty lists
        self.transcript_mods = {
            get_transcript_id_from_header(elem[0]): {sample: [] for sample in self.groups_df['Sample']} 
            for elem in fr.get_data()
        }
        # maps transcript_id to its sequence
        self.transcript_sequence = {get_transcript_id_from_header(elem[0]):elem[1] for elem in fr.get_data()}


    def process_header(self, line: str) -> tuple[int, int, int, int, int, int, list[int], list[str]]:
        """
        Processes the header line of the experimental data to identify column indices.
        
        :param line: Header line of the experimental data.
        :return: Indices for various fields such as peptide sequence, protein ID, etc.
        """
        exp_idx = []
        exp_names = []
        row = line.split('\t')

        for i, field in enumerate(row):
            if field == "Protein ID":
                prot_accession_idx = i
            elif field == "Peptide Sequence":
                pep_seq_idx = i
            elif field == "Modified Sequence":
                pep_mod_seq_idx = i
            elif field == "Assigned Modifications":
                assigned_mod_idx = i
            elif field == "Protein":
                protein_idx = i
            elif field == "Mapped Proteins":
                mapped_proteins_idx = i
            elif "Intensity" in field:
                if not "MaxLFQ Intensity" in field:
                    exp_idx.append(i)
                    exp_names.append(field.replace(" Intensity", "").strip())
        return prot_accession_idx, pep_seq_idx, pep_mod_seq_idx, assigned_mod_idx, protein_idx, mapped_proteins_idx, exp_idx, exp_names


    def get_exact_indexes(self, mod_sequence: str) -> list[int]:
        """
        Gets the exact indexes of the modifications in the peptide sequence.
        
        :param mod_sequence: The modified peptide sequence.
        :return: List of indices for modifications in the sequence.
        """
        indexes = []
        current_index = 1
        inside_brackets = False
        for i, char in enumerate(mod_sequence):
            if i == 0 and mod_sequence[i] == '[':
                indexes.append(0)
            if char == '[':
                inside_brackets = True
            elif char == ']':
                inside_brackets = False
                continue
            elif not inside_brackets and char.isalpha():
                if i + 1 < len(mod_sequence) and mod_sequence[i + 1] == '[':
                    indexes.append(current_index)
            if not inside_brackets:
                current_index += 1

        return indexes
    
    def get_offset(self, peptide: str, sequence: str) -> int:
        """
        Finds the offset of a peptide in the given sequence.
        
        :param peptide: The peptide sequence to search for.
        :param sequence: The full protein sequence.
        :return: The offset of the peptide in the sequence, or -1 if not found.
        """
        offset = 0
        if peptide in sequence:
            offset = sequence.index(peptide)
            return offset
        else: 
            return -1

    def process_modifications(self, mod_sequence: str, peptide_offset: int, sequence: str) -> list[str]:
        """
        Processes modifications in the peptide sequence and returns the modification strings.
        
        :param mod_sequence: The modified peptide sequence.
        :param peptide_offset: The offset of the peptide in the sequence.
        :param sequence: The full protein sequence.
        :return: A list of modifications in string format.
        """
        all_mods = []
        matches = re.findall(r'\[(-?\d+\.\d+?)\]', mod_sequence)
        peptide = re.sub(r'\[(-?\d+\.\d+?)\]', '', mod_sequence)
        peptide = peptide.replace("n", "")
        mod_sequence = mod_sequence.replace("n", "")
        aa_offsets = self.get_exact_indexes(mod_sequence)
        for i, match in enumerate(matches):
            if match in ms_fragger_mods \
                and ms_fragger_mods[match] in self.selected_ptms \
                and ms_fragger_mods[match] in modification_sites:
                # TODO: check with Chris on how to assign mods here 
                if aa_offsets[i] != 0:
                    modified_aa = peptide[aa_offsets[i]-1]
                else:
                    # skip n-term
                    continue
                    
                mod_string = f"{ms_fragger_mods[match]}({modified_aa})@{peptide_offset + aa_offsets[i]}"

                if modified_aa not in modification_sites[ms_fragger_mods[match]]:
                    #raise ValueError(f"Modification {ms_fragger_mods[match]} at AA {modified_aa} not allowed.")
                    logger.info(f"Modification {ms_fragger_mods[match]} at AA {modified_aa} not allowed.")
                    continue
                all_mods.append(mod_string)
        return all_mods
        
    def process_row(self, line: str, pep_seq_idx: int, pep_mod_seq_idx: int, assigned_mod_idx: int, protein_idx: int,
                    mapped_proteins_idx: int, exp_idx: list[int], exp_names: list[str]) -> None:
        """
        Processes a row of experimental data to identify and store modifications.
        
        :param line: A line from the experimental data.
        :param pep_seq_idx: Index of peptide sequence column.
        :param pep_mod_seq_idx: Index of modified peptide sequence column.
        :param assigned_mod_idx: Index of assigned modifications column.
        :param protein_idx: Index of protein column.
        :param mapped_proteins_idx: Index of mapped proteins column.
        :param exp_idx: List of indices for experimental data columns.
        :param exp_names: List of experiment names corresponding to the intensity columns.
        """
        row = line.split('\t')
        row = [field.strip() for field in row]

        if row[assigned_mod_idx] != '': 
            # peptide is modified 

            # get all associated proteins 
            mapped_proteins = []
            protein = row[protein_idx]
            mapped_proteins.append(get_transcript_id_from_header(protein))
            if row[mapped_proteins_idx] != '':
                further_proteins = row[mapped_proteins_idx].split(", ")
                for protein in further_proteins:
                    mapped_proteins.append(get_transcript_id_from_header(protein))

            # only keep significant ones
            filtered_proteins = [transcript_id for transcript_id in mapped_proteins if transcript_id in self.transcript_sequence]

            if not filtered_proteins:
                return

            # assign mod to each protein
            for transcript_id in filtered_proteins:
                offset = self.get_offset(row[pep_seq_idx], self.transcript_sequence[transcript_id])
                if offset == -1: 
                    self.peptide_not_in_sequence_count += 1
                    continue
                all_mods = self.process_modifications(row[pep_mod_seq_idx], offset, self.transcript_sequence[transcript_id])

                if not all_mods:
                    continue

                for i, idx in enumerate(exp_idx):
                    if row[idx] != "0.0":
                        if exp_names[i] not in self.transcript_mods[transcript_id]:
                            raise KeyError(f"Experiment {exp_names[i]} not found in groups.csv")
                        self.transcript_mods[transcript_id][exp_names[i]].extend(all_mods)
        
    def process_proteins(self, save_to_json: bool = False) -> dict:
        """
        Processes all the proteins in the experimental data and extracts modifications.
        
        :param save_to_json: Flag indicating whether to save the results to a JSON file.
        :return: Dictionary of transcript modifications.
        """
        self.peptide_not_in_sequence_count = 0

        with open(self.experiment_file, 'r', encoding="utf-8") as f:
            for line in tqdm(f, desc="Reading file", unit=" lines"):
                if line.startswith('Peptide Sequence'):
                    prot_accession_idx, pep_seq_idx, pep_mod_seq_idx, assigned_mod_idx, protein_idx, mapped_proteins_idx, exp_idx, exp_names = self.process_header(line)
                else:
                    self.process_row(line, pep_seq_idx, pep_mod_seq_idx, assigned_mod_idx, protein_idx, mapped_proteins_idx, exp_idx, exp_names)

        logger.info(f"{self.peptide_not_in_sequence_count} peptides where not in their mapped protein sequence.")

        if save_to_json:
            with open(self.output_path, 'w') as f:
                json.dump(self.transcript_mods, f, indent = 4)

        return self.transcript_mods