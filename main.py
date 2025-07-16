from src.config import load_toml_file
from src.logger import logger

from pipeline import *
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Data Pipeline CLI")
    parser.add_argument(
        "step", 
        type=str, 
        nargs="?",  # Makes the argument optional
        default="all",  # Default step when run without arguments
        help="Pipeline step to execute (ingest, process, store, external, all)"
    )
    parser.add_argument("--config", type=str, default="config.yaml", help="Path to configuration file")
    
    args = parser.parse_args()
    config = load_toml_file("config.toml")

    run_pipeline(args.step, config)