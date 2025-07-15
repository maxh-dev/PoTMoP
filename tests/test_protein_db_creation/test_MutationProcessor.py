import os 
import pytest
import json
import time 

from unittest.mock import Mock, MagicMock, patch
from src.protein_db_creation.mutation_processor import MutationProcessor
from src.util.fasta_handler import FastaReader
from src.util.protein_nameing import create_fasta_header
from pyensembl import EnsemblRelease
from tests.conftest import data_path_for_tests

mock_protein_reference_json = os.path.join(data_path_for_tests(), "mock_protein_reference.json")
mock_output_file = os.path.join(data_path_for_tests(), "mock_output.fasta")

@pytest.fixture
def mutation_info_dict_empty():
    mutation_info_dict = {
        "mutation_count": 0,
        "mutation_type_to_id": {}, 
        "sequences_to_id": {},
        "mutation_information": {}
    }

    return mutation_info_dict

@pytest.fixture
def mutation_info_dict_with_entry():
    mutation_info_dict = {
        "mutation_count": 2,
        "mutation_type_to_id": {
            "1_420|ENST2_42_G/A|4_E/K": 1,
            "X_240|ENST4_24_G/A|2_E/K": 2
        }, 
        "sequences_to_id": {
            "AAAK*": [1],
            "CK*": [2]
        },
        "mutation_information": {
            "1": {
                "transcript_id": "ENST2",
                "mutation_db_id": "ENST2_MUT1",
                "mutation_str": "1_420|ENST2_42_G/A|4_E/K",
                "sequence": "AAAK*"
            },
            "2": {
                "transcript_id": "ENST4",
                "mutation_db_id": "ENST4_MUT2",
                "mutation_str": "X_240|ENST4_24_G/A|2_E/K",
                "sequence": "CK*"
            }
        }
    }
    return mutation_info_dict

@pytest.fixture
def mutation_list():
    mutation_list = [('ENST00000234038', '2_241182658|ENST00000234038_946_G/A|315_E/K', 'MAAERGAGQQQSQEMMEVDRRVESEESGDEEGKKHSSGIVADLSEQSLKDGEERGEEDPEEEHELPVDMETINLDRDAEDVDLNHYRIGKIEGFEVLKKVKTLCLRQNLIKCIENLEELQSLRELDLYDNQIKKIENLEALTELEILDISFNLLRNIEGVDKLTRLKKLFLVNNKISKIENLSNLHQLQMLELGSNRIRAIENIDTLTNLESLFLGKNKITKLQNLDALTNLTVLSMQSNRLTKIEGLQNLVNLRELYLSHNGIEVIEGLENNNKLTMLDIASNRIKKIENISHLTELQEFWMNDNLLESWSDLDKLKGARSLETVYLERNPLQKDPQYRRKVMLALPSVRQIDATFVRF*')]
    return mutation_list


@pytest.fixture
def mutation_list_mock_big():
    mutation_list = [
        ('ENST1', '1_420|ENST1_42_G/A|4_E/K', 'AAAK*'),
        ('ENST1', '1_420|ENST1_42_G/A|4_E/K', 'AAAK*'),
        ('ENST2', '1_420|ENST2_42_G/A|4_E/K', 'AAAK*'),
        ('ENST3', '2_4200|ENST3_402_C/A|6_P/H', 'QQQSAH*'),
    ]
    return mutation_list

@pytest.fixture
def mutation_info_dict_for_mutation_list_mock_big():
    mutation_info_dict = {
        'mutation_count': 3,
        'mutation_type_to_id': {'1_420|ENST1_42_G/A|4_E/K': 1, '1_420|ENST2_42_G/A|4_E/K': 2, '2_4200|ENST3_402_C/A|6_P/H': 3},
        'sequences_to_id': {'AAAK*': [1, 2], 'QQQSAH*': [3]},
        'mutation_information': {'1': {'transcript_id': 'ENST1', 'mutation_db_id': 'ENST1_MUT1', 'mutation_str': '1_420|ENST1_42_G/A|4_E/K', 'sequence': 'AAAK*'}, '2': {'transcript_id': 'ENST2', 'mutation_db_id': 'ENST1_MUT1', 'mutation_str': '1_420|ENST2_42_G/A|4_E/K', 'sequence': 'AAAK*'}, '3': {'transcript_id': 'ENST3', 'mutation_db_id': 'ENST3_MUT3', 'mutation_str': '2_4200|ENST3_402_C/A|6_P/H', 'sequence': 'QQQSAH*'}}
    }
    return mutation_info_dict


@pytest.fixture
def processor(mutation_info_dict_empty):
    processor = MutationProcessor.__new__(MutationProcessor)

    processor.release_number = 111
    processor.fasta_file = mock_output_file

    processor.mutation_info_json = "mock_mutation_info.json"
    processor.mutation_info_dict = mutation_info_dict_empty

    processor.protein_reference_json = mock_protein_reference_json
    processor.protein_reference_dict = {}
    processor.mutations_list = []

    processor.ensembl_data = EnsemblRelease(111)
    processor.mutations_dict = {
        'key1':{'Chromosome': '7', 'PositionGenomic': '4023366', 'REF': 'C', 'ALT_1': 'C', 'ALT_2': 'T'},
        'key2':{'Chromosome': '2', 'PositionGenomic': '241182658','REF': 'G', 'ALT_1': 'G', 'ALT_2': 'A'}
    }

    processor.mutations_counter = 1
    
    yield processor

    if os.path.exists(mock_protein_reference_json):
        os.remove(mock_protein_reference_json)
    if os.path.exists(mock_output_file):
        os.remove(mock_output_file)

def test_translate_mutations(processor, mutation_list):
    processor.translate_mutations()

    assert processor.mutations_list == mutation_list


def test_insert_new_mutations(processor, mutation_list_mock_big):
    with patch("json.dump") as mock_json_dump:
        processor.mutations_list = mutation_list_mock_big
        processor.insert_new_mutations()

        assert mock_json_dump.call_count == 1
        assert processor.mutation_info_dict["mutation_count"] == 3
        assert processor.mutation_info_dict["sequences_to_id"] == {'AAAK*': [1, 2], 'QQQSAH*': [3]}
        assert processor.mutation_info_dict["mutation_type_to_id"] == {'1_420|ENST1_42_G/A|4_E/K': 1, '1_420|ENST2_42_G/A|4_E/K': 2, '2_4200|ENST3_402_C/A|6_P/H': 3}
        assert processor.mutation_info_dict["mutation_information"]["1"]["mutation_db_id"] == "ENST1_MUT1"
        assert processor.mutation_info_dict["mutation_information"]["2"]["mutation_db_id"] == "ENST1_MUT1"
        assert processor.mutation_info_dict["mutation_information"]["3"]["mutation_db_id"] == "ENST3_MUT3"

def test_insert_new_mutations_non_empty_mutation_info_dict(processor, mutation_info_dict_with_entry, mutation_list_mock_big):
    with patch("json.dump") as mock_json_dump:
        processor.mutation_info_dict = mutation_info_dict_with_entry
        processor.mutations_list = mutation_list_mock_big
        processor.insert_new_mutations()

        assert mock_json_dump.call_count == 1
        assert processor.mutation_info_dict["mutation_count"] == 4
        assert processor.mutation_info_dict["mutation_information"]["1"]["mutation_db_id"] == "ENST2_MUT1"
        assert processor.mutation_info_dict["mutation_information"]["2"]["mutation_db_id"] == "ENST4_MUT2"
        assert processor.mutation_info_dict["mutation_information"]["3"]["mutation_db_id"] == "ENST2_MUT1"
        assert processor.mutation_info_dict["mutation_information"]["4"]["mutation_db_id"] == "ENST3_MUT4"


def test_append_mutations_list_to_fasta(processor, mutation_list_mock_big, mutation_info_dict_for_mutation_list_mock_big):
    processor.mutations_list = mutation_list_mock_big
    processor.mutation_info_dict = mutation_info_dict_for_mutation_list_mock_big

    mock_sequence_to_header = {
        "AAAK*": None
    }

    with patch.object(FastaReader, "read_fasta") as mock_read_fasta, \
         patch.object(FastaReader, "get_sequence_to_header_dict", return_value=mock_sequence_to_header) as mock_get_seq_dict, \
         patch.object(processor, "append_entry") as mock_append_entry, \
         patch.object(processor.ensembl_data, "transcript_by_id") as mock_transcript_by_id:

        mock_transcript_by_id.return_value.gene_name = "PPP1R7"
        processor.append_mutations_list_to_fasta()

        assert len(mock_append_entry.call_args_list) == 1
        assert mock_append_entry.call_args_list[0].args == ('>sp|ENST3_MUT3|PPP1R7_HUMAN\n', 'QQQSAH*')
        assert "ENST3_MUT3" in processor.protein_reference_dict