import json 
import os
import tqdm
import numpy as np
import h5py
import json
import matplotlib.pyplot as plt
import h5py
import shutil
from Bio.PDB import PDBParser, MMCIFIO
import glob
import shutil

from src.logger import logger

class StructureIntegrator:
    def __init__(self, af_data_dir, predictions_dir, structure_reference_json_location):
        self.af_data_dir = af_data_dir
        self.predictions_dir = predictions_dir
        self.structure_reference_json_location = structure_reference_json_location

        with open(structure_reference_json_location, 'r') as f:
            self.structure_reference = json.load(f)

    def create_folders(self):
        cif_dir = os.path.join(self.af_data_dir, "cif")
        os.makedirs(cif_dir, exist_ok=True)
        pae_dir = os.path.join(self.af_data_dir, "pae")
        os.makedirs(pae_dir, exist_ok=True)

        return cif_dir, pae_dir
    
    def convert_pdb_to_cif(self, pdb_file, cif_file):
        # Parse the PDB file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)

        # Write the structure to mmCIF format
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(cif_file)

        logger.info(f"PDB file converted to CIF format and saved as {cif_file}")

    
    def af2_pae_json_to_hdf(self, input_path, output_path):
        with open(input_path) as f:
            data = json.load(f)[0]
        pae_matrix = np.array(data['predicted_aligned_error'])

        matrix_flat = np.array(pae_matrix).flatten()

        with h5py.File(output_path, 'w') as f:
            f.create_dataset('dist', data=matrix_flat)


    def load_af2_data(self):
        cif_dir, pae_dir = self.create_folders()

        for af2_dir_name in tqdm.tqdm(sorted(os.listdir(self.predictions_dir))):
            transcript_id = af2_dir_name
            alphafold_model_path = os.path.join(self.predictions_dir, af2_dir_name)
            model_pdb_path = os.path.join(alphafold_model_path, "selected_prediction.pdb")
            pae_json_path = os.path.join(alphafold_model_path, "predicted_aligned_error.json")

            cif_path = os.path.join(cif_dir, f"{transcript_id}.cif")
            self.convert_pdb_to_cif(model_pdb_path, cif_path)

            af2_hdf = os.path.join(pae_dir, f"pae_{transcript_id}.hdf")
            self.af2_pae_json_to_hdf(pae_json_path, af2_hdf)

            self.structure_reference["present"][transcript_id] = self.structure_reference["rerun"][transcript_id].copy()
            self.structure_reference["present"][transcript_id]["path"] = (cif_path, af2_hdf)
            del self.structure_reference["rerun"][transcript_id]

        with open(self.structure_reference_json_location , "w") as file:
            json.dump(self.structure_reference, file, indent=4) 