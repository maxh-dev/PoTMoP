import os
import json
import pytest
import numpy as np
import h5py
from unittest.mock import MagicMock, patch, mock_open
from src.alphafold_processing.BoltzStructureIntegrator import BoltzStructureIntegrator


@pytest.fixture
def mock_structure_reference():
    return {
        "present": {},
        "rerun": {
            "Protein1": {"sequence": "ACD", "pdb_id": None, "path": None},
            "Protein2": {"sequence": "MPX", "pdb_id": "UP02", "path": None}
        }
    }


@pytest.fixture
def integrator(mock_structure_reference):
    with patch("builtins.open", mock_open(read_data=json.dumps(mock_structure_reference))):
        return BoltzStructureIntegrator(
            af_data_dir="/mock/af_data",
            predictions_dir="/mock/predictions",
            structure_reference_json_location="/mock/structure_reference.json"
        )


@patch("os.makedirs")
def test_create_folders(mock_makedirs, integrator):
    with patch("os.path.join", side_effect=lambda *args: "/".join(args)):
        cif_dir, pae_dir = integrator.create_folders()

    assert cif_dir == "/mock/af_data/cif"
    assert pae_dir == "/mock/af_data/pae"
    mock_makedirs.assert_any_call(cif_dir, exist_ok=True)
    mock_makedirs.assert_any_call(pae_dir, exist_ok=True)

@patch("numpy.load")
@patch("h5py.File")
def test_boltz_pae_npz_to_hdf(mock_h5py_file, mock_np_load, integrator):
    mock_np_load.return_value = {"pae": np.array([[1.0, 2.0], [3.0, 4.0]])}

    mock_hdf = MagicMock()
    mock_h5py_file.return_value.__enter__.return_value = mock_hdf

    integrator.boltz_pae_npz_to_hdf("/mock/input.npz", "/mock/output.hdf")

    mock_np_load.assert_called_once_with("/mock/input.npz")
    mock_h5py_file.assert_called_once_with("/mock/output.hdf", "w")
    mock_hdf.create_dataset.assert_called_once()
    
    _, kwargs = mock_hdf.create_dataset.call_args
    np.testing.assert_array_equal(kwargs["data"], np.array([1.0, 2.0, 3.0, 4.0]))


@patch("os.makedirs")  
@patch("shutil.copy")  
@patch("glob.glob") 
@patch("os.listdir") 
@patch("builtins.open", new_callable=mock_open)
def test_load_boltz_data(mock_open_file, mock_listdir, mock_glob, mock_shutil_copy, mock_makedirs, integrator):
    """Test the load_boltz_data method with mocked filesystem operations."""
    
    mock_listdir.return_value = ["protein_1", "protein_2"]
    
    def glob_side_effect(path):
        if "*.cif" in path:
            return [path.replace("*", "model")]
        elif "*.npz" in path:
            return [path.replace("*", "model")]
        return []
    
    mock_glob.side_effect = glob_side_effect

    integrator.structure_reference = {"present": {}, "rerun": {"protein_1": {}, "protein_2": {}}}

    with patch.object(integrator, "boltz_pae_npz_to_hdf") as mock_boltz_pae:
        integrator.load_boltz_data()
    
        mock_shutil_copy.call_args_list[0][0] == ("mock/predictions/protein_1/model.cif", "/mock/cif/protein_1.cif")
        mock_shutil_copy.call_args_list[1][0] == ("mock/predictions/protein_2/model.cif", "/mock/cif/protein_2.cif")
        
        mock_boltz_pae.call_args_list[0][0] == ("mock/predictions/protein_1/model.npz", "/mock/pae/pae_protein_1.hdf")
        mock_boltz_pae.call_args_list[1][0] == ("mock/predictions/protein_2/model.npz", "/mock/pae/pae_protein_2.hdf")
