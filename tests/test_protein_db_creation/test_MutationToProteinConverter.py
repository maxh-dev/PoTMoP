import os 
import pytest
import json 

from unittest.mock import Mock, MagicMock, patch, PropertyMock
from src.protein_db_creation.mutation_processing.mutation_protein_converter import MutationToProteinConverter
from pyensembl import EnsemblRelease

@pytest.fixture
def mock_ensembl_data():
    mock_gene = MagicMock()
    mock_gene.transcripts = []

    mock_ensembl_data = MagicMock()
    mock_ensembl_data.genes_at_locus.return_value = [mock_gene]

    return mock_ensembl_data

@pytest.fixture
def converter(mock_ensembl_data):
    converter = MutationToProteinConverter.__new__(MutationToProteinConverter)

    converter.ensembl_data = mock_ensembl_data
    converter.chromosome = "1"
    converter.mutation_position_genomic = 125
    converter.ref_allele = "A"
    converter.alt_allele = "C"

    converter.matches = []
    converter.mutated_sequences = []

    return converter


@pytest.fixture
def mock_transcripts():
    #transcript_id
    mock_transcript = MagicMock()
    mock_transcript.start = 100
    mock_transcript.end = 200
    mock_transcript.transcript_id = "ENST001"
    mock_transcript.gene_name = "Gene1"
    mock_transcript.start_codon_complete = True
    mock_transcript.start_codon_spliced_offsets = [0, 1, 2]
    mock_transcript.stop_codon_spliced_offsets = [50, 51, 52]
    mock_transcript.spliced_offset.return_value = 25
    type(mock_transcript).coding_sequence = PropertyMock(return_value="ATG" * 10) 

    mock_transcript2 = MagicMock()
    mock_transcript2.start = 300
    mock_transcript2.end = 400
    mock_transcript2.transcript_id = "ENST002"
    mock_transcript2.start_codon_complete = False
    mock_transcript2.start_codon_spliced_offsets = [5, 6, 7]
    mock_transcript2.stop_codon_spliced_offsets = [60, 61, 62]
    mock_transcript2.spliced_offset.return_value = 30
    type(mock_transcript2).coding_sequence = PropertyMock(return_value="GTG" * 10)

    mock_transcript3 = MagicMock()
    mock_transcript3.start = 100
    mock_transcript3.end = 200
    mock_transcript3.transcript_id = "ENST003"
    mock_transcript3.start_codon_complete = True
    mock_transcript3.start_codon_spliced_offsets = [0, 1, 2]
    mock_transcript3.stop_codon_spliced_offsets = [50, 51, 52]
    mock_transcript3.spliced_offset.return_value = 25
    type(mock_transcript3).coding_sequence = PropertyMock(return_value="TGA" * 10) #REF does not match at position 25
    
    return [mock_transcript, mock_transcript2, mock_transcript3]

    

def test_find_matching_transcripts(converter, mock_transcripts):
    converter.ensembl_data.genes_at_locus.return_value[0].transcripts = mock_transcripts
    converter.find_matching_transcripts()


    assert len(converter.matches) == 1
    assert converter.matches[0][0].transcript_id == "ENST001"


def test_get_coding_sequences(converter, mock_transcripts):
    converter.matches = [(mock_transcripts[0], 25)]
    converter.get_coding_sequences()

    assert len(converter.mutated_sequences) == 1
    assert converter.mutated_sequences[0][1] == 'ATGATGATGATGATGATGATGATGCTGATG'

def test_get_coding_sequences_no_mutations(converter):
    converter.matches = []
    converter.get_coding_sequences()

    assert len(converter.mutated_sequences) == 0

def test_compute_mutated_sequences(converter, mock_transcripts):
    converter.ensembl_data.genes_at_locus.return_value[0].transcripts = mock_transcripts
    converter.compute_mutated_sequences()

    assert len(converter.mutated_sequences) == 1
    assert converter.mutated_sequences[0][1] == 'ATGATGATGATGATGATGATGATGCTGATG'

def test_compute_mutated_sequences_no_matching_transcripts(converter, mock_transcripts):
    converter.ensembl_data.genes_at_locus.return_value[0].transcripts = [mock_transcripts[1], mock_transcripts[2]]
    converter.compute_mutated_sequences()

    assert len(converter.mutated_sequences) == 0

def test_get_mutated_proteins_no_mutated_sequences(converter):
    converter.mutated_sequences = []

    result_list = converter.get_mutated_proteins()

    assert len(result_list) == 0

def test_get_mutated_proteins(converter, mock_transcripts):
    converter.mutated_sequences = [(mock_transcripts[0], 'ATGATGATGATGATGATGATGATGCTGATG', 25)]

    result_list = converter.get_mutated_proteins()

    assert len(result_list) == 1
    assert result_list == [('ENST001', '1_125|ENST001_25_A/C|8_M/L', 'MMMMMMMMLM*')]