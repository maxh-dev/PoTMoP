import json 
import os 
import shutil
import pandas as pd 
from src.logger import logger

def load_or_create_json(json_path):
    if os.path.exists(json_path):
        with open(json_path, "r") as file:
            return json.load(file)
    else:
        data = {}
        with open(json_path, "w") as file:
            json.dump(data, file, indent=4)
        return data
    

def get_transcript_id_from_header(header):
    return header.split("|")[1]

def safe_to_numeric(series):
    try:
        return pd.to_numeric(series)
    except ValueError:
        return series
    

def remove_duplicate_sequences(tuples_list):
    seen_sequences = set()
    unique_tuples = []
    duplicate_elems = []

    duplicate_count = 0

    for header, sequence in tuples_list:
        if sequence not in seen_sequences:
            unique_tuples.append((header, sequence))
            seen_sequences.add(sequence)
        else:
            logger.debug(f"Duplicate sequence found for {header[:-1]}")
            duplicate_count += 1
            duplicate_elems.append(header)
    logger.info(f"Duplicate count: {duplicate_count}. Duplicates were removed. {len(tuples_list) - duplicate_count} elements remain.")

    return unique_tuples, duplicate_elems


def copy_file_to_dir(source_file_path, destination_dir_path):
    if os.path.isfile(source_file_path):
        shutil.copy(source_file_path, destination_dir_path)
        logger.debug(f"Copied {source_file_path} to {destination_dir_path}")
    else:
        logger.info(f"{source_file_path} not found.")