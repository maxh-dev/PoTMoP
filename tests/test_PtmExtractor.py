import os
import pytest
import json
import pandas as pd
from unittest.mock import Mock, MagicMock, patch, mock_open
from src.ptm_extraction.PtmExtractor import PtmExtractor
from tests.conftest import data_path_for_tests


@pytest.fixture
def mock_groups_df():
    groups_df = [
        ["Sample1", "CAC"],
        ["Sample2", "CTR"]
    ]
    return pd.DataFrame(groups_df, columns=["Sample","groups"])

@pytest.fixture
def extractor(mock_groups_df):
    extractor = PtmExtractor.__new__(PtmExtractor)
    extractor.experiment_file = ""
    extractor.selected_ptms = ["Phospho", "Acetyl"]
    extractor.groups_df = mock_groups_df

    extractor.protein_fasta = None
    extractor.output_path = ""
    extractor.transcript_mods = {}

    yield extractor

@pytest.fixture
def mock_fasta_data():
    return [
        ("sp|Protein1|GP1_HUMAN", "ACD"),
        ("sp|Protein2|GP2_HUMAN", "MPP"),
        ("sp|Protein3|GP3_HUMAN", "MSQ")
    ]

@pytest.fixture
def mock_full_experiment_data():
    full_experiment_data = "Peptide Sequence\tModified Sequence\tPrev AA\tNext AA\tStart\tEnd\tPeptide Length\tCharges\tAssigned Modifications\tProtein\tProtein ID\tEntry Name\tGene\tProtein Description\tMapped Genes\tMapped Proteins\tSample1 Spectral Count\tSample2 Spectral Count\tSample1 Intensity\tSample2 Intensity\tSample1 MaxLFQ Intensity\tSample2 MaxLFQ Intensity\tSample1 Match Type\tSample2 Match Type\tSample1 Localization\tSample2 Localization\n" \
                            "ACS\tACS[79.9663]\t\t\t\t\t\t\t3S(79.9663)\tsp|Protein1|Protein1_HUMAN\tProtein1\t\t\t\t\tsp|Protein2|Protein2_HUMAN, sp|Protein3|Protein3_HUMAN\t\t\t0.5\t0.0\t\t\t\t\t\t\n" \
                            "MKD\tMK[42.0106]D\t\t\t\t\t\t\t2K(42.0106)\tsp|Protein2|Protein2_HUMAN\tProtein2\t\t\t\t\t\t\t\t0.1\t0.8\t\t\t\t\n"
    return full_experiment_data

@pytest.fixture
def mock_experiment_data():
    return [
        "ACS\tACS[79.9663]\t\t\t\t\t\t\t3S(79.9663)\tsp|Protein1|Protein1_HUMAN\tProtein1\t\t\t\t\tsp|Protein2|Protein2_HUMAN, sp|Protein3|Protein3_HUMAN\t\t\t0.5\t0.0\t\t\t\t\t\t",
        "MKD\tMK[42.0106]D\t\t\t\t\t\t\t2K(42.0106)\tsp|Protein2|Protein2_HUMAN\tProtein2\t\t\t\t\t\t\t\t0.1\t0.8\t\t\t\t"
    ]

def test_obtain_proteins(extractor, mock_fasta_data):
    with patch("src.ptm_extraction.PtmExtractor.FastaReader") as MockFastaReader:
        mock_fasta_reader = MockFastaReader.return_value
        mock_fasta_reader.get_data.return_value = mock_fasta_data
        
        extractor.obtain_proteins()

        proteins = ["Protein1", "Protein2", "Protein3"]
        for protein in proteins:
            assert protein in extractor.transcript_mods

        for protein in proteins:
            assert "Sample1" in extractor.transcript_mods[protein]
            assert "Sample2" in extractor.transcript_mods[protein]

        for protein in proteins:
            assert protein in extractor.transcript_sequence

def test_process_header(extractor):
    line = "Peptide Sequence\tModified Sequence\tPrev AA\tNext AA\tStart\tEnd\tPeptide Length\tCharges\tAssigned Modifications\tProtein\tProtein ID\tEntry Name\tGene\tProtein Description\tMapped Genes\tMapped Proteins\tSample1 Spectral Count\tSample2 Spectral Count\tSample1 Intensity\tSample2 Intensity\tSample1 MaxLFQ Intensity\tSample2 MaxLFQ Intensity\tSample1 Match Type\tSample2 Match Type\tSample1 Localization\tSample2 Localization"
    prot_accession_idx, pep_seq_idx, pep_mod_seq_idx, assigned_mod_idx, protein_idx, mapped_proteins_idx, exp_idx, exp_names = extractor.process_header(line)

    assert prot_accession_idx == 10 # same as Protein ID
    assert pep_seq_idx == 0
    assert pep_mod_seq_idx == 1
    assert assigned_mod_idx == 8
    assert protein_idx == 9
    assert mapped_proteins_idx == 15
    assert exp_idx == [18, 19]
    assert exp_names == ["Sample1", "Sample2"]


def test_get_exact_indexes(extractor):
    mod_sequence = "[79.9663]ACD[79.9663]M"
    indexes = extractor.get_exact_indexes(mod_sequence)

    assert indexes == [0,3]

def test_get_offset(extractor):
    sequence = "MSQACDMPPA"
    peptide = "MPP"
    offset = extractor.get_offset(peptide, sequence)

    assert offset == 6


def test_process_modifications(extractor):
    mod_sequence = "[42.0106]ACS[79.9663]MPP"
    peptide_offset = 3
    sequence = "MSKACSMPPA"
    modifications = extractor.process_modifications(mod_sequence, peptide_offset, sequence)

    assert 'Phospho(S)@6' in modifications


def test_process_row(extractor, mock_experiment_data):
    line = mock_experiment_data[0]
    extractor.transcript_sequence = {
        "Protein1": "ACSDCA",
        "Protein2": "AMDDCONRADACSCOOL",
    }
    extractor.transcript_mods["Protein1"] = {"Sample1":[], "Sample2":[]}
    extractor.transcript_mods["Protein2"] = {"Sample1":[], "Sample2":[]}
    extractor.process_row(line, 0, 1, 8, 9, 15, [18, 19], ["Sample1", "Sample2"])

    assert extractor.transcript_mods["Protein1"]["Sample1"] == ["Phospho(S)@3"]
    assert extractor.transcript_mods["Protein2"]["Sample1"] == ["Phospho(S)@13"]

    assert extractor.transcript_mods["Protein1"]["Sample2"] == []
    assert extractor.transcript_mods["Protein2"]["Sample2"] == []


def test_process_proteins(extractor, mock_full_experiment_data):
    extractor.transcript_sequence = {
        "Protein1": "ACSDCA",
        "Protein2": "AMKDCONRADACSCOOL",
    }
    extractor.transcript_mods["Protein1"] = {"Sample1":[], "Sample2":[]}
    extractor.transcript_mods["Protein2"] = {"Sample1":[], "Sample2":[]}

    with patch("builtins.open", mock_open(read_data=mock_full_experiment_data)) as mock_file:
        with patch("json.dump") as mock_json_dump:
            with patch("tqdm.tqdm", side_effect=lambda x, **kwargs: x):
                extractor.experiment_file = "mock_experiment_file.tsv"
                extractor.output_path = "mock_output.json"
                
                result = extractor.process_proteins(save_to_json=True)
                
                assert mock_json_dump.call_count == 1 
                assert isinstance(result, dict)
                assert "Sample1" in result["Protein2"]
                assert "Sample2" in result["Protein2"]
                assert result["Protein2"]["Sample1"] == ['Phospho(S)@13', 'Acetyl(K)@3']
                assert result["Protein2"]["Sample2"] == ['Acetyl(K)@3']
