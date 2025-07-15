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

class BoltzStructureIntegrator:
    def __init__(self, af_data_dir: str, predictions_dir: str, structure_reference_json_location: str):
        """
        Initializes the BoltzStructureIntegrator class.
        
        :param af_data_dir: Directory where AlphaFold data is stored.
        :param predictions_dir: Directory containing structure predictions.
        :param structure_reference_json_location: Path to the structure reference JSON file.
        """
        self.af_data_dir = af_data_dir
        self.predictions_dir = predictions_dir
        self.structure_reference_json_location = structure_reference_json_location

        with open(structure_reference_json_location, 'r') as f:
            self.structure_reference = json.load(f)

    def create_folders(self) -> tuple[str, str]:
        """
        Creates required directories for CIF and PAE files.
        
        :return: Tuple containing paths to CIF and PAE directories.
        """
        cif_dir = os.path.join(self.af_data_dir, "cif")
        os.makedirs(cif_dir, exist_ok=True)
        pae_dir = os.path.join(self.af_data_dir, "pae")
        os.makedirs(pae_dir, exist_ok=True)

        return cif_dir, pae_dir
    

    def boltz_pae_npz_to_hdf(self, input_path: str, output_path: str) -> None:
        """
        Converts a Boltz-1 PAE NPZ file to HDF format.
        
        :param input_path: Path to the NPZ file.
        :param output_path: Path to the output HDF file.
        """
        npz_file = np.load(input_path)
        with h5py.File(output_path, "w") as hdf_file:
            for key in npz_file:
                if key == "pae":
                    hdf_file.create_dataset('dist', data=np.array(npz_file[key]).flatten())


    def load_boltz_data(self) -> None:
        """
        Loads structure data from Boltzmann predictions and updates structure reference JSON.
        """
        cif_dir, pae_dir = self.create_folders()

        for boltz_dir_name in tqdm.tqdm(sorted(os.listdir(self.predictions_dir))):
            transcript_id = boltz_dir_name
            boltz_model_path = os.path.join(self.predictions_dir, boltz_dir_name)
            cifs = glob.glob(os.path.join(boltz_model_path, "*.cif"))
            npzs = glob.glob(os.path.join(boltz_model_path, "*.npz"))
            if len(cifs) != 1 or len(npzs) != 1: 
                logger.info(f"Skipping {transcript_id}, files not complete.")
                continue 
            model_cif_path = cifs[0]
            pae_npz_path = npzs[0]

            cif_path = os.path.join(cif_dir, f"{transcript_id}.cif")
            shutil.copy(model_cif_path, cif_path)

            af2_hdf = os.path.join(pae_dir, f"pae_{transcript_id}.hdf")
            self.boltz_pae_npz_to_hdf(pae_npz_path, af2_hdf)

            self.structure_reference["present"][transcript_id] = self.structure_reference["rerun"][transcript_id].copy()
            self.structure_reference["present"][transcript_id]["path"] = (cif_path, af2_hdf)
            del self.structure_reference["rerun"][transcript_id]

        with open(self.structure_reference_json_location , "w") as file:
            json.dump(self.structure_reference, file, indent=4) 

