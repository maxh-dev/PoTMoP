import os
import json
import toml


def init_base_folder_structure(config):
    os.makedirs(config["pipeline"]["base_dir"], exist_ok=True)
    os.makedirs(config["pipeline"]["base_external_dir"], exist_ok=True)
    os.makedirs(config["pipeline"]["base_af_pdb_dir"], exist_ok=True)
    os.makedirs(config["pipeline"]["base_af_custom_dir"], exist_ok=True)

def init_mutation_info_json(config):
    mutation_info_dict = {
        "mutation_count": 0,
        "mutation_type_to_id": {}, 
        "sequences_to_id": {},
        "mutation_information": {}
    }
    file_path = config["pipeline"]["base_mutation_json"]
    if os.path.exists(file_path):
        os.remove(file_path)

    with open(file_path, "w") as f:
        json.dump(mutation_info_dict, f)
    
def update_config_toml(config):
    with open("config.toml", "w") as file:
        toml.dump(config, file)

def get_run_id(config):
    last_id = config["pipeline"]["last_run_id"]
    id = last_id + 1
    config["pipeline"]["last_run_id"] = id
    config["pipeline"]["selected_run_id"] = id
    update_config_toml(config)
    return id

def init_run_folder_structure(run_data_dir, run_name):
    run_config = {}

    analysis_dir = os.path.join(run_data_dir, run_name)  

    protein_db_dir = os.path.join(analysis_dir, "protein_db")  
    protein_db = os.path.join(protein_db_dir, "protein_db.fas")  
    protein_reference = os.path.join(protein_db_dir, "protein_reference.json")  
    protein_db_filtered = os.path.join(protein_db_dir, "protein_db_filtered.fas")  
    muttaions_db_dir = os.path.join(protein_db_dir, "mutation_files") 

    ms_data_db_dir = os.path.join(analysis_dir, "ms_data")  
    ms_log2_fold_change = os.path.join(ms_data_db_dir, "log2_fold_change.csv")  
    ms_corrected_p_values = os.path.join(ms_data_db_dir, "corrected_p_values.csv")  
    experiment_data = os.path.join(ms_data_db_dir, "combined_modified_peptide.tsv")  
    experiment_groups = os.path.join(ms_data_db_dir, "groups.csv")  
    ptm_extractor_output = os.path.join(ms_data_db_dir, "ptms_for_proteins.json")  
    ptm_formatter_output = os.path.join(ms_data_db_dir, "ptms")  

    alphafold_dir = os.path.join(analysis_dir, "alphafold") 
    fasta_for_execution_dir = os.path.join(alphafold_dir, "fasta_for_execution")  
    predictions_dir = os.path.join(alphafold_dir, "predictions")   

    structure_reference = os.path.join(analysis_dir, "structure_reference.json")  

    protein_structure_ptms_csv = os.path.join(analysis_dir, "protein_structure_ptms.csv")  

    os.makedirs(analysis_dir, exist_ok=True)

    os.makedirs(protein_db_dir, exist_ok=True)
    os.makedirs(muttaions_db_dir, exist_ok=True)

    os.makedirs(ms_data_db_dir, exist_ok=True)
    os.makedirs(ptm_formatter_output, exist_ok=True)

    os.makedirs(alphafold_dir, exist_ok=True)
    os.makedirs(fasta_for_execution_dir, exist_ok=True)
    os.makedirs(predictions_dir, exist_ok=True)

    run_config["protein_db"] = {}
    run_config["protein_db"]["fasta_location"] = protein_db
    run_config["protein_db"]["protein_reference"] = protein_reference
    run_config["protein_db"]["mutation_directory"] = muttaions_db_dir

    run_config["protein_filter"] = {}
    run_config["protein_filter"]["fasta_location"] = protein_db
    run_config["protein_filter"]["filtered_fasta_location"] = protein_db_filtered
    run_config["protein_filter"]["log2_fold_change"] = ms_log2_fold_change
    run_config["protein_filter"]["corrected_p_values"] = ms_corrected_p_values
    run_config["protein_filter"]["p_value"] = 0.05

    run_config["structure_fetcher"] = {}
    run_config["structure_fetcher"]["fasta_for_execution"] = fasta_for_execution_dir
    run_config["structure_fetcher"]["structure_reference_json_location"] = structure_reference

    run_config["structure_integrator"] = {}
    run_config["structure_integrator"]["structure_reference_json_location"] = structure_reference
    run_config["structure_integrator"]["predictions_dir"] = predictions_dir

    run_config["ptm_extractor"] = {}
    run_config["ptm_extractor"]["experiment_data_path"] = experiment_data
    run_config["ptm_extractor"]["experiment_groups_path"] = experiment_groups
    run_config["ptm_extractor"]["fasta_path"] = protein_db_filtered
    run_config["ptm_extractor"]["interim_output_path"] = ptm_extractor_output
    run_config["ptm_extractor"]["output_path"] = ptm_formatter_output

    run_config["analysis"] = {}
    run_config["analysis"]["structure_reference_json_location"] = structure_reference
    run_config["analysis"]["ptm_dir"] = ptm_formatter_output
    run_config["analysis"]["protein_structure_ptms"] = protein_structure_ptms_csv

    run_config_file = os.path.join(analysis_dir, "run_config.json") 
    with open(run_config_file, 'w') as f:
        json.dump(run_config, f, indent=4)

    
