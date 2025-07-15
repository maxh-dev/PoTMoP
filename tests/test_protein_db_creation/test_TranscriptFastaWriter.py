import os 
import pytest
from unittest.mock import Mock, MagicMock, patch
import json
from src.protein_db_creation.transcript_fasta_writer import TranscriptFastaWriter
from tests.conftest import data_path_for_tests

mock_protein_reference_json = os.path.join(data_path_for_tests(), "mock_protein_reference.json")
mock_output_file = os.path.join(data_path_for_tests(), "mock_output.fasta")

@pytest.fixture
def writer():
    writer = TranscriptFastaWriter.__new__(TranscriptFastaWriter)
    
    
    writer.release_number = 100
    writer.sequence_to_uniprot = {"MS": {"id": "U00000", "name": "T0000_HUMAN"}}
    writer.protein_reference_dict = {}
    writer.protein_reference_json = mock_protein_reference_json
    writer.output_file = mock_output_file
    
    mock_ensembl_data = Mock()

    mock_transcript_1 = Mock()
    mock_transcript_1.transcript_id = "ENST000001"
    mock_transcript_1.transcript_name = "ENST000001"
    mock_transcript_1.gene_id = "ENSG000001"
    mock_transcript_1.gene_name = "BRCA1"
    mock_transcript_1.protein_id = "P38398"
    mock_transcript_1.is_protein_coding = True
    mock_transcript_1.biotype = "protein_coding"
    mock_transcript_1.complete = True
    mock_transcript_1.coding_sequence = "ATGAGTUAG"

    mock_transcript_2 = Mock()
    mock_transcript_2.transcript_id = "ENST000002"
    mock_transcript_2.transcript_name = "ENST000002"
    mock_transcript_2.gene_id = "ENSG000002"
    mock_transcript_2.gene_name = "TP53"
    mock_transcript_2.protein_id = "P04637"
    mock_transcript_2.is_protein_coding = True
    mock_transcript_2.biotype = "protein_coding"
    mock_transcript_2.complete = False
    mock_transcript_2.coding_sequence = "ATGCAG"

    mock_transcript_3 = Mock()
    mock_transcript_3.transcript_id = "ENST000003"
    mock_transcript_3.transcript_name = "ENST000003"
    mock_transcript_3.gene_id = "ENSG000003"
    mock_transcript_3.gene_name = "BTG2"
    mock_transcript_3.protein_id = "P17354"
    mock_transcript_3.is_protein_coding = True
    mock_transcript_3.biotype = "protein_coding"
    mock_transcript_3.complete = True
    mock_transcript_3.coding_sequence = "CAGAGTUAG"

    mock_ensembl_data.transcripts.return_value = [mock_transcript_1, mock_transcript_2, mock_transcript_3]

    writer.ensembl_data = mock_ensembl_data

    yield writer

    if os.path.exists(mock_protein_reference_json):
        os.remove(mock_protein_reference_json)
    if os.path.exists(mock_output_file):
        os.remove(mock_output_file)


def test_writer_initialization(writer):
    assert writer.release_number == 100
    assert writer.sequence_to_uniprot == {"MS": {"id": "U00000", "name": "T0000_HUMAN"}}
    assert writer.protein_reference_dict == {}

def test_writer_is_valid_aa_sequence(writer):
    seq1 = "ABC*"
    seq2 = "ABC"
    seq3 = "ACD*"

    assert writer.is_valid_aa_sequence(seq1) == False
    assert writer.is_valid_aa_sequence(seq2) == False
    assert writer.is_valid_aa_sequence(seq3) == True

def test_write_transcripts_to_fasta(writer):
    writer.write_transcripts_to_fasta()
    
    with open(mock_protein_reference_json) as f:
        mock_protein_reference = json.load(f)

    assert mock_protein_reference["ENST000001"] == {"uniprot_id": "U00000","mutation": None}
    assert mock_protein_reference["ENST000003"] == {"uniprot_id": None,"mutation": None}

    with open(mock_output_file, "r") as f:
        mock_output_lines = f.readlines()

    assert mock_output_lines[0] == '>sp|ENST000001|BRCA1_HUMAN\n'
    assert mock_output_lines[1] == 'MS\n'
    assert mock_output_lines[2] == '>sp|ENST000003|BTG2_HUMAN\n'
    assert mock_output_lines[3] == 'QS\n'