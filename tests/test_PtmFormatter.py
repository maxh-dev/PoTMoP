import os
import pytest
import json
import pandas as pd
from unittest.mock import Mock, MagicMock, patch, mock_open
from src.ptm_formatting.PtmFormatter import PtmFormatter
from tests.conftest import data_path_for_tests

@pytest.fixture
def mock_groups_df():
    groups_df = [
        ["Sample1", "CAC"],
        ["Sample2", "CTR"]
    ]
    return pd.DataFrame(groups_df, columns=["Sample","groups"])

@pytest.fixture
def formatter(mock_groups_df):
    formatter = PtmFormatter.__new__(PtmFormatter)
    formatter.ptms_for_proteins = {}
    formatter.all_groups = mock_groups_df['groups'].unique()
    formatter.file_to_group = dict(zip(mock_groups_df['Sample'], mock_groups_df['groups']))
    formatter.output_dir = ""
    yield formatter


@pytest.fixture
def mock_results_dict():
    return {
        "Protein1": {"Phospho(S)@3": {"all":{"CAC": 1, "CTR": 0}, "distinct": {"CAC": 1, "CTR": 0}}},
        "Protein2": {"Acetyl(K)@5": {"all":{"CAC": 1, "CTR": 2}, "distinct": {"CAC": 1, "CTR": 1}}},
    }

@pytest.fixture
def mock_results_dict_same_aa_diff_mods():
    return {
        "Protein1": {"Phospho(S)@3": {"all":{"CAC": 1, "CTR": 0}, "distinct": {"CAC": 1, "CTR": 0}}, "Acetyl(S)@3": {"all":{"CAC": 1, "CTR": 1}, "distinct": {"CAC": 1, "CTR": 1}}},
    }


@pytest.fixture
def mock_ptms_for_proteins():
    return {
        "Protein1": {"Sample1": ["Phospho(S)@3"]},
        "Protein2": {"Sample1": ["Acetyl(K)@5"], "Sample2": ["Acetyl(K)@5"]}
    }


def test_extract_ptm_info(formatter):
    input_str = "Phospho(S)@3"
    expected_output = {"type": "Phospho", "aa": "S", "pos": "3"}
    assert formatter.extract_ptm_info(input_str) == expected_output

def test_compute_columns(formatter):
    supported_ptms = ["ac", "p"]
    columns, ptm_group_list = formatter.compute_columns(supported_ptms)

    assert columns == ['protein_id', 'AA', 'position', 'CAC_ac', 'CAC_ac_distinct', 'CAC_p', 'CAC_p_distinct', 'CTR_ac', 'CTR_ac_distinct', 'CTR_p', 'CTR_p_distinct']
    assert ptm_group_list == ['CAC_ac', 'CAC_ac_distinct', 'CAC_p', 'CAC_p_distinct', 'CTR_ac', 'CTR_ac_distinct', 'CTR_p', 'CTR_p_distinct']


def test_build_ptm_csv(formatter, mock_results_dict):
    df = formatter.build_ptm_csv(mock_results_dict)

    assert not df.empty
    assert list(df.columns) == ['protein_id', 'AA', 'position', 'CAC_ac', 'CAC_ac_distinct', 'CAC_p', 'CAC_p_distinct', 'CTR_ac', 'CTR_ac_distinct', 'CTR_p', 'CTR_p_distinct']
    assert df.shape == (2, 11) 

    assert df.loc[0, "protein_id"] == "Protein1"
    assert df.loc[0, "AA"] == "S"
    assert df.loc[0, "position"] == "3"
    assert df.loc[0, "CAC_ac"] == 0
    assert df.loc[0, "CAC_p"] == 1
    assert df.loc[0, "CTR_ac"] == 0
    assert df.loc[0, "CTR_p"] == 0
    assert df.loc[0, "CAC_ac_distinct"] == 0
    assert df.loc[0, "CAC_p_distinct"] == 1
    assert df.loc[0, "CTR_ac_distinct"] == 0
    assert df.loc[0, "CTR_p_distinct"] == 0
    
    assert df.loc[1, "protein_id"] == "Protein2"
    assert df.loc[1, "AA"] == "K"
    assert df.loc[1, "position"] == "5"
    assert df.loc[1, "CAC_ac"] == 1
    assert df.loc[1, "CAC_p"] == 0
    assert df.loc[1, "CTR_ac"] == 2
    assert df.loc[1, "CTR_p"] == 0
    assert df.loc[1, "CAC_ac_distinct"] == 1
    assert df.loc[1, "CAC_p_distinct"] == 0
    assert df.loc[1, "CTR_ac_distinct"] == 1
    assert df.loc[1, "CTR_p_distinct"] == 0


def test_build_ptm_csv_same_aa_diff_mods(formatter, mock_results_dict_same_aa_diff_mods):
    df = formatter.build_ptm_csv(mock_results_dict_same_aa_diff_mods)

    assert not df.empty
    assert list(df.columns) == ['protein_id', 'AA', 'position', 'CAC_ac', 'CAC_ac_distinct', 'CAC_p', 'CAC_p_distinct', 'CTR_ac', 'CTR_ac_distinct', 'CTR_p', 'CTR_p_distinct']

    assert df.shape == (1, 11) 

    assert df.loc[0, "protein_id"] == "Protein1"
    assert df.loc[0, "AA"] == "S"
    assert df.loc[0, "position"] == "3"
    assert df.loc[0, "CAC_ac"] == 1
    assert df.loc[0, "CAC_p"] == 1
    assert df.loc[0, "CTR_ac"] == 1
    assert df.loc[0, "CTR_p"] == 0
    assert df.loc[0, "CAC_ac_distinct"] == 1
    assert df.loc[0, "CAC_p_distinct"] == 1
    assert df.loc[0, "CTR_ac_distinct"] == 1
    assert df.loc[0, "CTR_p_distinct"] == 0

def test_process_ptm_dictionary(formatter, mock_ptms_for_proteins):
    formatter.ptms_for_proteins = mock_ptms_for_proteins
    
    with patch("pandas.DataFrame") as MockDataFrame:
        formatter.process_ptm_dictionary()

        assert "to_csv" in MockDataFrame.return_value.mock_calls[-1][0]
        assert MockDataFrame.return_value.mock_calls[-1][1] == ('ptms.csv',)
        #mock_to_csv.assert_called_once()

        assert MockDataFrame.call_args_list[0][1] == {'columns': ['protein_id', 'AA', 'position', 'CAC_ac', 'CAC_ac_distinct', 'CAC_p', 'CAC_p_distinct', 'CTR_ac', 'CTR_ac_distinct', 'CTR_p', 'CTR_p_distinct']}
        df_content = MockDataFrame.call_args_list[0][0][0] 
        assert df_content == [{'protein_id': 'Protein1', 'AA': 'S', 'position': '3', 'CAC_ac': 0, 'CAC_ac_distinct': 0, 'CAC_p': 1, 'CAC_p_distinct': 1, 'CTR_ac': 0, 'CTR_ac_distinct': 0, 'CTR_p': 0, 'CTR_p_distinct': 0}, {'protein_id': 'Protein2', 'AA': 'K', 'position': '5', 'CAC_ac': 1, 'CAC_ac_distinct': 1, 'CAC_p': 0, 'CAC_p_distinct': 0, 'CTR_ac': 1, 'CTR_ac_distinct': 1, 'CTR_p': 0, 'CTR_p_distinct': 0}]