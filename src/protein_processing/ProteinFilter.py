import pandas as pd
import json
from src.util.fasta_handler import FastaReader, FastaWriter
from src.util.general import load_or_create_json

from src.logger import logger

class ProteinFilter():

    def __init__(self, fasta_file: str, filtered_fasta: str, protein_reference_json: str, 
                 log2_fold_change: str, corrected_p_values: str, take_all_proteins: bool, p_thresh: float = None):
        """
        Initializes the SignificantProteinFilter class.
        
        Parameters:
        fasta_file (str): Path to the input FASTA file.
        filtered_fasta (str): Path to the output filtered FASTA file.
        protein_reference_json (str): Path to the protein reference JSON file.
        log2_fold_change (str): Path to the CSV file containing log2 fold change data.
        corrected_p_values (str): Path to the CSV file containing corrected p-values.
        p_thresh (float): Threshold for significance based on corrected p-value.
        """
        self.fasta_file = fasta_file
        self.filtered_fasta = filtered_fasta
        self.protein_reference_json = protein_reference_json
        self.protein_reference_dict = load_or_create_json(protein_reference_json)
        self.log2_fold_change = pd.read_csv(log2_fold_change)
        self.corrected_p_values = pd.read_csv(corrected_p_values)
        self.take_all_proteins = take_all_proteins
        self.p_thresh = p_thresh
        if take_all_proteins:
            logger.info("All proteins will be transferred.")
        else: 
            if not self.p_thresh:
                raise ValueError("p_thresh can't be None, if take_all_proteins is false.") 
            else:
                logger.info(f"Filtering for significant proteins with p_thresh = {self.p_thresh}")

        self.significant_protein_dict = {}

        self.obtain_proteins()

    def obtain_proteins(self) -> None:
        """
        Reads the input FASTA file and stores the data in a list of tuples (header, sequence).
        """
        fr = FastaReader(self.fasta_file)
        fr.read_fasta()
        self.fasta_data = fr.get_data()

    def filter_significant_proteins(self) -> dict[str, dict]:
        """
        Filters proteins based on significance criteria and stores relevant data.
        
        Returns:
        Dict[str, Dict]: Dictionary containing significant protein data.
        """
        merged_data = pd.merge(self.log2_fold_change, self.corrected_p_values, on="Protein ID", how="inner")
        protein_dict = merged_data.set_index("Protein ID")[["corrected_p_value", "log2_fold_change"]].to_dict(orient="index")

        sig_count = 0
        structure_present_count = 0

        filtered_proteins = []
        for header, sequence in self.fasta_data:            
            header_elems = header.split("|")
            transcript_id = header_elems[1]
            gene_id = header_elems[2]

            # check if transcript_id was found in probe
            if not transcript_id in protein_dict:
                continue
            
            if "P" in transcript_id:
                logger.info(transcript_id)
            if len(transcript_id) < 15: 
                logger.info(f"transcript_id {transcript_id} excluded, as its name does not fit the format.")
                continue

            p_val = protein_dict[transcript_id]["corrected_p_value"]
            log2_fold_change = protein_dict[transcript_id]["log2_fold_change"]

            if ((p_val < self.p_thresh) and not self.take_all_proteins) or self.take_all_proteins == True: 
                sig_count += 1
                # save all significant proteins in significant_protein_dict (with uniprot id if present)
                if self.protein_reference_dict[transcript_id]["uniprot_id"]:
                    structure_present_count += 1 
                    self.significant_protein_dict[transcript_id] = {
                        "structure_present": True,
                        "pdb_id": self.protein_reference_dict[transcript_id]["uniprot_id"],
                        "sequence": sequence,
                    }
                else: 
                    self.significant_protein_dict[transcript_id] = {
                        "structure_present": False,
                        "pdb_id": None,
                        "sequence": sequence,
                    }

                # update the protein_reference_dict
                self.protein_reference_dict[transcript_id]["significant"] = {
                        "p_val": p_val,
                        "log2_fold_change": log2_fold_change
                    }
                # append to filtered_proteins fasta     
                filtered_proteins.append((header, sequence))

        logger.info(f"Found {sig_count} significant proteins of which {structure_present_count} structures are present. {sig_count-structure_present_count} structure predictions nessesary.")

        fw = FastaWriter(self.filtered_fasta, filtered_proteins)
        fw.write_fasta()
        logger.info(f"Filterd FASTA was written to {self.filtered_fasta}")

        with open(self.protein_reference_json , "w") as file:
            json.dump(self.protein_reference_dict, file, indent=4)  
        return self.significant_protein_dict
