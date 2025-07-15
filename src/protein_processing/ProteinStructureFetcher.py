import os
import shutil

from structuremap.processing import download_alphafold_cif, download_alphafold_pae
import json 
from src.util.fasta_handler import FastaWriter
from src.util.general import copy_file_to_dir

from src.logger import logger

class ProteinStructureFetcher:
    def __init__(self, protein_dict: dict[str, dict], structure_prediction_length_limit: int, 
                 structure_output_dir: str, structure_custom_output_dir: str, protein_output_path: str, af_fasta_path: str):
        """
        Initializes the ProteinStructureFetcher class.
        
        Parameters:
        protein_dict (Dict[str, Dict]): Dictionary containing protein data of the significant proteins 
            found in SignificantProteinFilter. 
        structure_prediction_length_limit (int): Maximum length for structure prediction.
        structure_output_dir (str): Directory for PDB structure output files.
        structure_custom_output_dir (str): Directory for custom structure output files.
        protein_output_path (str): Path to save the protein output JSON.
        af_fasta_path (str): Path for storing AlphaFold input FASTA files.
        """
        self.protein_dict = protein_dict
        self.structure_prediction_length_limit = structure_prediction_length_limit
        self.structure_output_dir = structure_output_dir
        self.structure_custom_output_dir = structure_custom_output_dir
        self.protein_output_path = protein_output_path
        self.af_fasta_path = af_fasta_path

        self.pdb_present_protein_ids = []

        self.protein_output_dict = {}
        self.protein_output_dict["present"] = {}
        self.protein_output_dict["rerun"] = {}

        self.rerun_alphafold_list = []
         
        self.pdb_cif_dir = os.path.join(self.structure_output_dir, 'cif')
        self.pdb_pae_dir = os.path.join(self.structure_output_dir, 'pae')
        os.makedirs(self.pdb_cif_dir, exist_ok=True)
        os.makedirs(self.pdb_pae_dir, exist_ok=True)

        self.custom_cif_dir = os.path.join(self.structure_custom_output_dir, 'cif')
        self.custom_pae_dir = os.path.join(self.structure_custom_output_dir, 'pae')
        os.makedirs(self.custom_cif_dir, exist_ok=True)
        os.makedirs(self.custom_pae_dir, exist_ok=True)


    def transcript_structure_available(self, transcript_id):
        cif_file_path = os.path.join(self.custom_cif_dir, f"{transcript_id}.cif")
        return os.path.exists(cif_file_path)


    def process_protein_dict(self) -> None:
        """
        Processes the protein dictionary to categorize proteins based on structure availability.
        """
        for transcript_id in self.protein_dict: 
            elem = self.protein_dict[transcript_id]
            
            if elem["structure_present"]:
                # structure_present is true, if the sequence was matched to a uniprot sequence.
                # In that case the Alphafold Protein DB will be checked for the structure files. 
                self.pdb_present_protein_ids.append(elem["pdb_id"])
            elif self.transcript_structure_available(transcript_id): # This checks if the structure of the transcript was predicted in a previous run with Boltz-1. 
                cif_data_path = os.path.join(self.custom_cif_dir, f"{transcript_id}.cif")  
                pae_data_path = os.path.join(self.custom_pae_dir, f"pae_{transcript_id}.hdf")  
                self.protein_output_dict["present"][transcript_id] = {
                    "sequence": elem["sequence"],
                    "pdb_id": None, 
                    "path": (cif_data_path, pae_data_path)
                } 
            else: 
                # Other wise the structure will need to be predicted.  
                self.protein_output_dict["rerun"][transcript_id] = {
                    "sequence": elem["sequence"],
                    "pdb_id": None, 
                    "path": None
                }
                self.rerun_alphafold_list.append((transcript_id, elem["sequence"]))

    def add_downloads_to_output_dict(self, valid_existing_proteins_cif: list[str], invalid_proteins_cif: list[str]) -> None:
        """
        Updates the protein output dictionary with valid and invalid structure downloads.
        
        Parameters:
        valid_existing_proteins_cif (List[str]): List of valid downloaded or existing CIF structures.
        invalid_proteins_cif (List[str]): List of invalid CIF structures.
        """
        #TODO might inlcude pae info as well 

        for key in self.protein_dict: 
            elem = self.protein_dict[key]

            if elem["structure_present"]:
                if elem["pdb_id"] in valid_existing_proteins_cif:
                    cif_data_path = os.path.join(self.pdb_cif_dir, f"{elem['pdb_id']}.cif")  
                    pae_data_path = os.path.join(self.pdb_pae_dir, f"pae_{elem['pdb_id']}.hdf")  
                    self.protein_output_dict["present"][key] = {
                        "sequence": elem["sequence"],
                        "pdb_id": elem["pdb_id"], 
                        "path": (cif_data_path, pae_data_path)
                    }
                elif elem["pdb_id"] in invalid_proteins_cif:
                    self.protein_output_dict["rerun"][key] = {
                        "sequence": elem["sequence"],
                        "pdb_id": elem["pdb_id"], 
                        "path": None
                    }
                else:
                    logger.error(f"The elem {elem['pdb_id']} could not be found in either list.")
            else:
                continue 

    
    def fetch_structures(self) -> dict[str, dict]:
        """
        Fetches protein structure data using AlphaFold and updates the output dictionary.
        
        Returns:
        Dict[str, Dict]: Updated dictionary containing structure information.
        """
        logger.info("Fetching Proteins cif and pae now.")
        valid_proteins_cif, invalid_proteins_cif, existing_proteins_cif = download_alphafold_cif(
            proteins=self.pdb_present_protein_ids, out_folder=self.pdb_cif_dir)
        
        valid_proteins_pae, invalid_proteins_pae, existing_proteins_pae = download_alphafold_pae(
            proteins=self.pdb_present_protein_ids, out_folder=self.pdb_pae_dir)
        
        valid_existing_proteins_cif = valid_proteins_cif + existing_proteins_cif

        self.add_downloads_to_output_dict(valid_existing_proteins_cif, invalid_proteins_cif)

        return self.protein_output_dict

    
    def save_to_json(self) -> None:
        """
        Saves the protein structure reference dictionary to a JSON file.
        """
        logger.info(f"Saving protein structural reference to {self.protein_output_path} now.")
        with open(self.protein_output_path , "w") as json_file:
            json.dump(self.protein_output_dict, json_file, indent=4)

    def create_fasta_for_alphafold_execution(self) -> None:
        """
        Creates a FASTA file for rerun sequences for AlphaFold execution.
        """
        rerun_list = []
        for key in self.protein_output_dict["rerun"]:
            sequence = self.protein_output_dict["rerun"][key]["sequence"]
            rerun_list.append((key, sequence))
        FastaWriter(f"{self.af_fasta_path}/alphafold.fasta", rerun_list).write_fasta()

    def create_fasta_for_boltz_execution(self) -> None:
        """
        Creates individual FASTA files for rerun sequences for Boltz-1 execution.
        """  
        for transcript_id in self.protein_output_dict["rerun"]:
            header = "A|protein|"
            sequence = self.protein_output_dict["rerun"][transcript_id]["sequence"]
            if len(sequence) <= self.structure_prediction_length_limit:
                FastaWriter(f"{self.af_fasta_path}/{transcript_id}.fasta", [(header, sequence)]).write_fasta()
            else: 
                logger.info(f"Protein {transcript_id} skipped, exceeds the structure_prediction_length_limit.")
 