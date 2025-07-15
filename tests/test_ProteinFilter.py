import os 
import pytest
import json 
import pandas as pd 

from unittest.mock import Mock, MagicMock, patch
from src.protein_processing.ProteinFilter import ProteinFilter
from tests.conftest import data_path_for_tests

mock_protein_reference_json = os.path.join(data_path_for_tests(), "mock_protein_reference.json")
mock_output_file = os.path.join(data_path_for_tests(), "mock_output.fasta")


@pytest.fixture
def corrected_p_values():
    corrected_p_values = (
        ["Protein1", 0.04],
        ["Protein2", 0.03],
        ["Protein3", 0.07],
    )

    corrected_p_values = pd.DataFrame(
        data=corrected_p_values,
        columns=["Protein ID", "corrected_p_value"],
    )

    return corrected_p_values

@pytest.fixture
def df_log2_fold_change():
    df_log2_fold_change = (
        ["Protein1", 0.023],
        ["Protein2", -0.016],
        ["Protein3", 0.003],
    )

    df_log2_fold_change = pd.DataFrame(
        data=df_log2_fold_change,
        columns=["Protein ID", "log2_fold_change"],
    )

    return df_log2_fold_change

@pytest.fixture
def protein_reference_dict():
    protein_reference_dict = {
        "Protein1": {"uniprot_id": None, "mutation": None},
        "Protein2": {"uniprot_id": "UP02", "mutation": None},
        "Protein3": {"uniprot_id": None, "mutation": None},
        "Protein4": {"uniprot_id": None, "mutation": None}
    }

    return protein_reference_dict

@pytest.fixture
def fasta_data():
    fasta_data = [
        ("sp|Protein1|GP1_HUMAN", "ACD"),
        ("sp|Protein2|GP2_HUMAN", "MPP"),
        ("sp|Protein3|GP3_HUMAN", "MSQ"),
        ("sp|Protein4|GP4_HUMAN", "RKG")
    ]

    return fasta_data

@pytest.fixture
def filter(protein_reference_dict, corrected_p_values, df_log2_fold_change, fasta_data):
    filter = ProteinFilter.__new__(ProteinFilter)

    filter.filtered_fasta = mock_output_file
    filter.protein_reference_json = mock_protein_reference_json
    filter.protein_reference_dict = protein_reference_dict
    filter.log2_fold_change = df_log2_fold_change
    filter.corrected_p_values = corrected_p_values
    filter.take_all_proteins = False
    filter.p_thresh = 0.05

    filter.fasta_data = fasta_data
    filter.significant_protein_dict = {}
    
    yield filter

    if os.path.exists(mock_protein_reference_json):
        os.remove(mock_protein_reference_json)
    if os.path.exists(mock_output_file):
        os.remove(mock_output_file)


def test_filter_significant_proteins(filter):
    significant_protein_dict = filter.filter_significant_proteins()

    assert significant_protein_dict["Protein1"]["structure_present"] == False 
    assert significant_protein_dict["Protein2"]["structure_present"] == True 
    assert significant_protein_dict["Protein2"]["pdb_id"] == 'UP02' 

    with open(mock_protein_reference_json) as f:
        mock_protein_reference = json.load(f)

    assert mock_protein_reference["Protein1"]["significant"]["p_val"] == 0.04
    assert mock_protein_reference["Protein1"]["significant"]["log2_fold_change"] == 0.023

    assert mock_protein_reference["Protein2"]["uniprot_id"] == "UP02"

    with open(mock_output_file, "r") as f:
        mock_output_lines = f.readlines()

    assert mock_output_lines[0] == '>sp|Protein1|GP1_HUMAN\n'
    assert mock_output_lines[1] == 'ACD\n'
    assert mock_output_lines[2] == '>sp|Protein2|GP2_HUMAN\n'
    assert mock_output_lines[3] == 'MPP\n'


def test_filter_all_proteins(filter):
    filter.take_all_proteins = True
    significant_protein_dict = filter.filter_significant_proteins()

    for protein_name in ["Protein1","Protein2","Protein3"]:
        assert protein_name in significant_protein_dict