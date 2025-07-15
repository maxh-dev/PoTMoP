from src.config import open_run_config
from src.logger import logger

from src.pipeline_initialzation import init_base_folder_structure, init_mutation_info_json, init_run_folder_structure, get_run_id
from src.protein_db_creation.transcript_fasta_writer import TranscriptFastaWriter
from src.protein_db_creation.mutation_processor import MutationProcessor
from src.protein_processing.ProteinFilter import ProteinFilter
from src.protein_processing.ProteinStructureFetcher import ProteinStructureFetcher
from src.alphafold_processing.BoltzStructureIntegrator import BoltzStructureIntegrator
from src.ptm_extraction.PtmExtractor import PtmExtractor
from src.ptm_formatting.PtmFormatter import PtmFormatter
from src.alphafold_processing.af_data_formatter import AlphaFoldDataFormatter


def prepare_pipeline(config):
    logger.info(f"Initialise the pipeline.")
    init_base_folder_structure(config)
    init_mutation_info_json(config)
    
def prepare_run(config):
    logger.info(f"Initialise a run of the pipeline.")
    
    run_id = get_run_id(config)
    init_run_folder_structure(config["pipeline"]["run_dir"], f"run_{run_id}")

def create_db(config):
    logger.info(f"Custom protein DB creation.")
    run_config = open_run_config(config)
    # Writes all Ensemble protein coding transcripts to a fasta file
    
    tfw = TranscriptFastaWriter(
        release_number = config["protein_db"]["ensemble_release_number"], 
        sequence_to_uniprot = config["protein_db"]["sequence_to_uniprot_location"],
        output_file = run_config["protein_db"]["fasta_location"],
        protein_reference_json = run_config["protein_db"]["protein_reference"]
    )
    tfw.write_transcripts_to_fasta()
    
    # Adds SNP-Mutations from all VCF and MAF files in the mutation_directory to the fasta file 
    mp = MutationProcessor(
        release_number = config["protein_db"]["ensemble_release_number"], 
        mutation_info_json = config["pipeline"]["base_mutation_json"],
        mutation_directory = run_config["protein_db"]["mutation_directory"],
        fasta_file = run_config["protein_db"]["fasta_location"],
        protein_reference_json = run_config["protein_db"]["protein_reference"]
    )
    mp.translate_mutations()
    mp.insert_new_mutations()
    mp.append_mutations_list_to_fasta()


def filter_significant_proteins(config):
    logger.info("Filtering for significant proteins.")
    run_config = open_run_config(config)
    spf = ProteinFilter(
        fasta_file = run_config["protein_filter"]["fasta_location"],
        filtered_fasta = run_config["protein_filter"]["filtered_fasta_location"],
        protein_reference_json = run_config["protein_db"]["protein_reference"],
        log2_fold_change = run_config["protein_filter"]["log2_fold_change"],
        corrected_p_values = run_config["protein_filter"]["corrected_p_values"],
        take_all_proteins = True,
        p_thresh = run_config["protein_filter"]["p_value"]
    )
    significant_protein_dict = spf.filter_significant_proteins()

    logger.info("Fetching existing proteins.")
    psf = ProteinStructureFetcher(
        protein_dict = significant_protein_dict,
        structure_prediction_length_limit = 2400,
        structure_output_dir = config["pipeline"]["base_af_pdb_dir"],
        structure_custom_output_dir = config["pipeline"]["base_af_custom_dir"],
        protein_output_path = run_config["structure_fetcher"]["structure_reference_json_location"],
        af_fasta_path = run_config["structure_fetcher"]["fasta_for_execution"]
    )
    psf.process_protein_dict()
    psf.fetch_structures()
    psf.save_to_json()
    psf.create_fasta_for_boltz_execution()


def process_new_alphafold2_runs(config):
    logger.info("Process alphafold runs.")
    run_config = open_run_config(config)
    si = BoltzStructureIntegrator(
        af_data_dir = config["pipeline"]["base_af_custom_dir"], 
        predictions_dir = run_config["structure_integrator"]["predictions_dir"],
        structure_reference_json_location = run_config["structure_integrator"]["structure_reference_json_location"]
    )
    si.load_boltz_data()


def preprocess_ms_data(config):
    logger.info("Preprocess MS data.")
    run_config = open_run_config(config)

    # Given experiment data and a protein fasta file, this call returns all PTMs for all given Proteins
    ppe = PtmExtractor(
        experiment_file = run_config["ptm_extractor"]["experiment_data_path"], 
        selected_ptms = config["msfragger_ptms"]["selected_ptms"],
        experiment_groups_file = run_config["ptm_extractor"]["experiment_groups_path"],
        protein_fasta = run_config["ptm_extractor"]["fasta_path"],
        output_path = run_config["ptm_extractor"]["interim_output_path"],
    )
    ppe.obtain_proteins()
    ppe.process_proteins(save_to_json=True)

    # Formats PTMs to CSV for StructureMap
    # only support Phospo and Acetyl
    pf = PtmFormatter(
        ptms_for_proteins_path=run_config["ptm_extractor"]["interim_output_path"],
        experiment_groups_file = run_config["ptm_extractor"]["experiment_groups_path"],
        output_dir=run_config["ptm_extractor"]["output_path"],
    )
    pf.process_ptm_dictionary()


def structure_map_preprocessing(config):
    logger.info("Run StructureMap preprocessing steps.")
    # This inlcudes converting the cifs to rows in a df, computing metrics (like pPSE) and merging the ptm df into the structural df. 
    run_config = open_run_config(config)

    afdf = AlphaFoldDataFormatter(
        custom_dir = config["pipeline"]["base_af_custom_dir"],
        pdb_dir = config["pipeline"]["base_af_pdb_dir"],
        ptm_dir = run_config["ptm_extractor"]["output_path"],
        structure_reference_json = run_config["structure_fetcher"]["structure_reference_json_location"],
        protein_structure_ptms_csv = run_config["analysis"]["protein_structure_ptms"],
        annotation_config = config["af_annotation_settings"],
    )
    afdf.format_af_data()
    afdf.annotate_af_accessibility()
    afdf.annoatte_idrs()
    afdf.annotate_ptms()
    afdf.save_alphafold_ptms()
    
    return



def run_pipeline(step, config):
    """Run a specific pipeline step based on user input."""
    steps = {
        "init": prepare_pipeline,
        "prepare": prepare_run,
        "create_db": create_db,
        "filter_significant_proteins": filter_significant_proteins,
        "process_new_alphafold2_runs": process_new_alphafold2_runs,
        "preprocess": preprocess_ms_data,
        "structure_map_preprocessing": structure_map_preprocessing
    }

    if step in steps:
        steps[step](config)
    else:
        print(f"Invalid step: {step}. Available steps: {list(steps.keys())}")
