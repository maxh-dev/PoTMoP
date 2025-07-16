import os 
import pytest
import json 

from unittest.mock import mock_open, patch
from src.protein_db_creation.mutation_processing.maf_parser import MAFParser


@pytest.fixture
def parser():
    parser = MAFParser.__new__(MAFParser)
    parser.variants = {}
    return parser


@pytest.fixture
def mock_maf_content():
    return """Chromosome\tStart_position\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tVariant_Type\tVariant_Classification\n""" \
           "1\t12345\tA\tT\tT\tSNP\tMissense_Mutation\n" \
           "2\t67890\tG\tC\tC\tSNP\tSilent\n" \
           "3\t13579\tT\tA\tA\tINDEL\tFrameshift_Mutation\n"


def test_parse(parser, mock_maf_content):    
    with patch("builtins.open", mock_open(read_data=mock_maf_content)):
        parser.parse("dummy_path.maf")
    
    expected_variant_key = "1:12345_A/T/T"
    expected_variant_data = {
        'Chromosome': '1',
        'PositionGenomic': '12345',
        'REF': 'A',
        'ALT_1': 'T',
        'ALT_2': 'T',
    }
    
    assert expected_variant_key in parser.variants
    assert parser.variants[expected_variant_key] == expected_variant_data
    
    assert "2:67890_G/C/C" not in parser.variants
    assert "3:13579_T/A/A" not in parser.variants
