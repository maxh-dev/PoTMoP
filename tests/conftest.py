import os 

import pytest
from unittest import mock
import json

def data_path_for_tests():
    test_data_path = os.path.join(os.getcwd(), "tests", "data")
    os.makedirs(test_data_path, exist_ok=True)
    return test_data_path
