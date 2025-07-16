import os 
import toml
import json 

def load_toml_file(file_path):
    try:
        with open(file_path, 'r') as file:
            config = toml.load(file)
            return config
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
    except toml.TomlDecodeError as e:
        print(f"Error decoding TOML file: {e}")


def load_json_file(file_path):
    with open(file_path, 'r') as f:
        json_obj = json.load(f)

    return json_obj


def open_run_config(global_config):
    run_name = f"run_{global_config['pipeline']['selected_run_id']}"
    run_config_path = os.path.join(global_config["pipeline"]["run_dir"], run_name, "run_config.json")  

    return load_json_file(run_config_path)

    