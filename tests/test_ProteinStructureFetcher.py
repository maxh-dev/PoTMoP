
import os 
import pytest
import json 
import pandas as pd 

from unittest.mock import Mock, MagicMock, patch
import tempfile
from src.protein_processing.ProteinStructureFetcher import ProteinStructureFetcher
from tests.conftest import data_path_for_tests

mock_protein_reference_json = os.path.join(data_path_for_tests(), "mock_protein_reference.json")
mock_output_file = os.path.join(data_path_for_tests(), "mock_output.fasta")


@pytest.fixture
def significant_protein_dict():
    significant_protein_dict = {
        "Protein1": {'structure_present': False, 'pdb_id': None, 'sequence': 'ACDACD'},
        "Protein2": {'structure_present': True, 'pdb_id': 'UP02', 'sequence': 'MPP'},
        "Protein3": {'structure_present': True, 'pdb_id': 'UP03', 'sequence': 'MPX'}
    }

    return significant_protein_dict


@pytest.fixture
def fetcher(significant_protein_dict):
    fetcher = ProteinStructureFetcher.__new__(ProteinStructureFetcher)

    fetcher.protein_dict = significant_protein_dict

    fetcher.protein_output_dict = {}
    fetcher.protein_output_dict["present"] = {}
    fetcher.protein_output_dict["rerun"] = {}
    
    fetcher.pdb_present_protein_ids = []
    fetcher.rerun_alphafold_list = []

    fetcher.structure_prediction_length_limit = 5

    fetcher.pdb_cif_dir = "/cif"
    fetcher.pdb_pae_dir = "/pae"

    fetcher.custom_cif_dir = "/cust_cif"
    fetcher.custom_pae_dir = "/cust_pae"
    
    yield fetcher

    if os.path.exists(mock_protein_reference_json):
        os.remove(mock_protein_reference_json)
    if os.path.exists(mock_output_file):
        os.remove(mock_output_file)

def test_process_protein_dict(fetcher):
    fetcher.process_protein_dict()

    assert len(fetcher.pdb_present_protein_ids) == 2
    assert fetcher.pdb_present_protein_ids == ['UP02', 'UP03']

    assert len(fetcher.rerun_alphafold_list) == 1
    assert fetcher.rerun_alphafold_list == [('Protein1', 'ACDACD')]

    assert 'Protein1' in fetcher.protein_output_dict['rerun']
    assert fetcher.protein_output_dict['present'] == {}

def test_process_protein_dict_with_prev_structure_present(fetcher):
    with patch.object(fetcher, 'transcript_structure_available', side_effect=lambda transcript_id: transcript_id == 'Protein1'):
        fetcher.process_protein_dict()

        assert len(fetcher.pdb_present_protein_ids) == 2
        assert fetcher.pdb_present_protein_ids == ['UP02', 'UP03']

        assert len(fetcher.rerun_alphafold_list) == 0

        assert 'Protein1' in fetcher.protein_output_dict['present']
        assert fetcher.protein_output_dict['present']['Protein1']['path'] == ('/cust_cif/Protein1.cif', '/cust_pae/pae_Protein1.hdf')


def test_add_downloads_to_output_dict(fetcher):
    fetcher.protein_output_dict["rerun"]['Protein1'] = None
    fetcher.add_downloads_to_output_dict(
        valid_existing_proteins_cif=['UP02'],
        invalid_proteins_cif=['UP03']
    )
    assert 'Protein2' in fetcher.protein_output_dict['present']
    assert fetcher.protein_output_dict['present']['Protein2']['path'] == ('/cif/UP02.cif', '/pae/pae_UP02.hdf')

    assert 'Protein1' in fetcher.protein_output_dict['rerun']
    assert 'Protein3' in fetcher.protein_output_dict['rerun']

def test_fetch_structures(fetcher):
    with patch("src.protein_processing.ProteinStructureFetcher.download_alphafold_cif") as mock_cif, \
         patch("src.protein_processing.ProteinStructureFetcher.download_alphafold_pae") as mock_pae, \
         patch("os.makedirs"):
        
        mock_cif.return_value = (["UP02"], ["UP03"], [])
        mock_pae.return_value = (["UP02"], ["UP03"], [])

        fetcher.fetch_structures()

        assert 'Protein2' in fetcher.protein_output_dict['present']
        assert fetcher.protein_output_dict['present']['Protein2']['path'] == ('/cif/UP02.cif', '/pae/pae_UP02.hdf')

        assert 'Protein3' in fetcher.protein_output_dict['rerun']
        assert fetcher.protein_output_dict['rerun']['Protein3']['pdb_id'] == 'UP03'
        assert fetcher.protein_output_dict['rerun']['Protein3']['path'] is None


def test_create_fasta_for_boltz_execution(fetcher):
    fetcher.protein_output_dict["rerun"]["Protein1"] = {'sequence': 'ACDACD', 'pdb_id': None, 'path': None}
    fetcher.protein_output_dict["rerun"]["Protein3"] = {'sequence': 'MPX', 'pdb_id': 'UP03', 'path': None}

    with tempfile.TemporaryDirectory() as temp_dir, patch("src.protein_processing.ProteinStructureFetcher.FastaWriter") as mock_fasta_writer:
        fetcher.af_fasta_path = temp_dir      
        fetcher.create_fasta_for_boltz_execution()
        assert mock_fasta_writer.call_count == 1  

        mock_fasta_writer_call = mock_fasta_writer.call_args_list[0][0]
        assert mock_fasta_writer_call[0] == f"{temp_dir}/Protein3.fasta"
        assert mock_fasta_writer_call[1] == [("A|protein|", "MPX")]