import os
import pandas as pd
import re
import json
from imblearn.under_sampling import RandomUnderSampler

from src.logger import logger

class PtmFormatter:
    def __init__(self, ptms_for_proteins_path: str, experiment_groups_file: str, output_dir: str):
        """
        Initializes the PtmFormatter with PTM data, experimental groups, and output directory.
        
        :param ptms_for_proteins_path: Path to JSON file containing PTM data for proteins.
        :param experiment_groups_file: Path to CSV file containing experiment group mappings.
        :param output_dir: Directory where output files will be saved.
        """
        with open(ptms_for_proteins_path) as f:
            self.ptms_for_proteins = json.load(f)
        df = pd.read_csv(experiment_groups_file)
        df = self.perform_sampling(df, "RandomUnderSampler")
        self.all_groups = df['groups'].unique()
        self.file_to_group = dict(zip(df['Sample'], df['groups']))

        self.output_dir = output_dir

    def perform_sampling(self, df, sampling_type):
        logger.info(f"Sample counts before sampling: {df['groups'].value_counts()}")

        if sampling_type == "RandomUnderSampler":
            undersampler = RandomUnderSampler(sampling_strategy='auto', random_state=42)
        else:
            raise ValueError("Only RandomUnderSampler supported.")

        X_resampled, y_resampled = undersampler.fit_resample(df[['Sample']], df['groups'])

        df_resampled = pd.DataFrame({'Sample': X_resampled['Sample'], 'groups': y_resampled})

        logger.info(f"Sample counts after sampling: {df_resampled['groups'].value_counts()}")

        return df_resampled
    
    def extract_ptm_info(self, input_str: str) -> dict[str, str]:
        """
        Extracts PTM type, amino acid, and position from a given input string.
        
        :param input_str: String containing PTM information.
        :return: Dictionary containing 'type', 'aa', and 'pos' keys with corresponding values.
        """
        patterns = {
            "type": r"^([A-Za-z]+)\(",
            "aa": r"\((\w)\)",
            "pos": r"@(\d+)"
        }
        return {key: re.search(pattern, input_str).group(1) for key, pattern in patterns.items()}
    

    def compute_columns(self, supported_ptms: list[str]) -> tuple[list[str], list[str]]:
        """
        Computes column names for the output CSV file based on supported PTMs and experiment groups.
        
        :param supported_ptms: List of supported PTM types.
        :return: Tuple containing a list of column names and a list of PTM group columns.
        """
        columns = ['protein_id', 'AA', 'position']  
        ptm_group_list = []
        for group in self.all_groups:
            for ptm in supported_ptms:
                ptm_group_list.append(f"{group}_{ptm}")
                ptm_group_list.append(f"{group}_{ptm}_distinct")
        columns.extend(ptm_group_list)
        return columns, ptm_group_list 
        

    
    def build_ptm_csv(self, result_dict: dict[str, dict[str, dict[str, int]]]) -> pd.DataFrame:
        """
        Builds a DataFrame from the processed PTM dictionary.
        
        :param result_dict: Processed PTM data.
        :return: DataFrame containing PTM information formatted for CSV output.
        """
        supported_ptms = ['ac', 'p']
        ptms_name_matching = {"Phospho":"p", "Acetyl": "ac"}
        columns, ptm_group_list = self.compute_columns(supported_ptms)  
        data_rows = []

        for protein_id in result_dict:
            for ptm in result_dict[protein_id]:
                ptm_info = self.extract_ptm_info(ptm)
                data = {'protein_id': protein_id, 'AA': ptm_info['aa'], 'position': ptm_info['pos']}
                for ptm_group in ptm_group_list:
                    data[ptm_group] = 0
                
                if ptm_info["type"] in ptms_name_matching:
                    ptm_type = ptms_name_matching[ptm_info["type"]]
                else:
                    logger.error(f"PTM Type {ptm_info['type']} not found. Skipping")
                    continue
                for group in result_dict[protein_id][ptm]["all"]:
                    ptm_id = f"{group}_{ptm_type}"
                    data[ptm_id] = result_dict[protein_id][ptm]["all"][group]
                for group in result_dict[protein_id][ptm]["distinct"]:
                    ptm_id = f"{group}_{ptm_type}_distinct"
                    data[ptm_id] = result_dict[protein_id][ptm]["distinct"][group]
                data_rows.append(data)

        df = pd.DataFrame(data_rows, columns=columns)
        result_df = df.groupby(['protein_id', 'AA', 'position'], as_index=False)[ptm_group_list].sum()
        return result_df
    
    def process_ptm_dictionary(self) -> None:
        """
        Processes the PTM dictionary, converts it into a DataFrame, and saves it as a CSV file.
        """
        result_dict = {}
        for protein_id in self.ptms_for_proteins:
            temp_dict = {}
            for subject in self.ptms_for_proteins[protein_id]:
                if subject not in self.file_to_group:
                    # subject is not within the sampled group, so skip that subject
                    continue

                for ptm_name in self.ptms_for_proteins[protein_id][subject]:
                    if ptm_name not in temp_dict:
                       temp_dict[ptm_name] = {
                           "all": {elem:0 for elem in self.all_groups},
                           "distinct": {elem:0 for elem in self.all_groups}
                       }
                    temp_dict[ptm_name]["all"][self.file_to_group[subject]] += 1

                for ptm_name in set(self.ptms_for_proteins[protein_id][subject]):
                    temp_dict[ptm_name]["distinct"][self.file_to_group[subject]] += 1
                    
            result_dict[protein_id] = temp_dict

        logger.info(f"Building PTM DFs, convert to CSV and saveing to {self.output_dir}")

        df = self.build_ptm_csv(result_dict)
        df.to_csv(os.path.join(self.output_dir, "ptms.csv"), index=False)
        


