import os
import json 
import pandas as pd
import re
import numpy as np
import tqdm
import Bio.PDB.MMCIF2Dict
import structuremap.utils
import h5py
from Bio.SeqUtils import seq1

from src.util.general import safe_to_numeric
from structuremap.processing import partition_df_by_prots, get_neighbors
from structuremap.processing import  get_smooth_score, annotate_proteins_with_idr_pattern, get_extended_flexible_pattern, get_proximity_pvals, perform_enrichment_analysis, perform_enrichment_analysis_per_protein, evaluate_ptm_colocalization, extract_motifs_in_proteome
from structuremap.plotting import plot_enrichment, plot_ptm_colocalization

from src.logger import logger

def annotate_accessibility(
    df: pd.DataFrame,
    max_dist: float,
    max_angle: float,
    structure_reference: dict,
) -> pd.DataFrame:
    """
    Half sphere exposure as calculated in
    https://onlinelibrary.wiley.com/doi/10.1002/prot.20379
    but with paired aligned error metric included.

    Parameters
    ----------
    df : pd.DataFrame
        pd.DataFrame of formatted AlphaFold data across various proteins.
        Such a dataframe is gerated by format_alphafold_data.
    max_dist : float
        Float specifying the maximum distance between two amino acids.
    max_angle : float
        Float specifying the maximum angle (in degrees) between two
        amino acids.
    structure_reference: : dict
        Path to the dict that maps protein ids to their cif and pae files.
    Returns
    -------
    : pd.DataFrame
        Dataframe repporting the number of neighboring amino acids at the
        specified maximum distance and angle per protein, amino acid and
        position.
    """
    # Adapted from structuremap 0.0.10, original function: structuremap.processing.annotate_accessibility
    # Modifications: custom pae path
    proteins = list()
    AA = list()
    AA_p = list()
    a_AA = list()
    for df_prot in partition_df_by_prots(df):
        protein_accession = df_prot.protein_id.values[0]
        _, pae_path = structure_reference[protein_accession]["path"]
        if pae_path is not None:
            with h5py.File(pae_path) as hdf_root:
                error_dist = hdf_root['dist'][...]
            size = int(np.sqrt(len(error_dist)))
            error_dist = error_dist.reshape(size, size)
            use_pae = 'pae'
        else:
            error_dist = np.zeros((df_prot.shape[0], df_prot.shape[0]))
            use_pae = 'nopae'
        idx_list = np.arange(0, df_prot.shape[0])
        res_a = get_neighbors(
            idx_list=idx_list,
            coord_a=np.vstack([df_prot.x_coord_ca.values,
                              df_prot.y_coord_ca.values,
                              df_prot.z_coord_ca.values]).T,
            coord_b=np.vstack([df_prot.x_coord_cb.values,
                              df_prot.y_coord_cb.values,
                              df_prot.z_coord_cb.values]).T,
            coord_c=np.vstack([df_prot.x_coord_c.values,
                              df_prot.y_coord_c.values,
                              df_prot.z_coord_c.values]).T,
            coord_n=np.vstack([df_prot.x_coord_n.values,
                              df_prot.y_coord_n.values,
                              df_prot.z_coord_n.values]).T,
            # If this step is slow, consider avoiding the vstack to create new arrays
            # Alternatively, it might be faster to use e.g. df[["x", "y", "z"]].values
            # as pandas might force this into a view rather than a new array
            position=df_prot.position.values,
            error_dist=error_dist,
            max_dist=max_dist,
            max_angle=max_angle)
        proteins.append(df_prot.protein_id.values)
        # using numeracal prot_numbers might be better.
        # In general it is good practice to reduce strings/objects in arrays/dfs
        # as much possible. Especially try to avoid repetetion of such types and
        # just use indices and a reference array. Rarely do you need this actual
        # values anyways.
        AA.append(df_prot.AA.values)
        AA_p.append(df_prot.position.values)
        a_AA.append(res_a)
    proteins = np.concatenate(proteins)
    AA = np.concatenate(AA)
    AA_p = np.concatenate(AA_p)
    a_AA = np.concatenate(a_AA)
    accessibility_df = pd.DataFrame({'protein_id': proteins,
                                     'AA': AA, 'position': AA_p})
    attribute_name = f'nAA_{max_dist}_{max_angle}_{use_pae}'
    accessibility_df[attribute_name] = a_AA
    return accessibility_df, attribute_name 

class AlphaFoldDataFormatter:
    def __init__(
        self,
        custom_dir: str,
        pdb_dir: str,
        ptm_dir: str,
        structure_reference_json: str,
        protein_structure_ptms_csv: str,
        annotation_config: dict[str, any]
    ) -> None:
        """
        Initialize the AlphaFoldDataFormatter class.

        Args:
            custom_dir (str): Directory for custom data.
            pdb_dir (str): Directory containing structual files.
            ptm_dir (str): Directory containing PTM data.
            structure_reference_json (str): Path to structure reference JSON file.
            protein_structure_ptms_csv (str): Path to save annotated PTM data.
            annotation_config (Dict[str, Any]): Configuration for annotations.
        """
        self.custom_dir = custom_dir
        self.pdb_dir = pdb_dir
        self.ptm_dir = ptm_dir
        self.structure_reference_json = structure_reference_json
        self.protein_structure_ptms_csv = protein_structure_ptms_csv
        self.annotation_config = annotation_config

        with open(structure_reference_json, 'r') as f:
            self.structure_reference = json.load(f)

        self.protein_ids = []
        self.alphafold_annotation = None
        self.alphafold_accessibility = None
        self.full_sphere_exposure_name = None
        self.part_sphere_exposure_name = None
        self.alphafold_accessibility_smooth_pattern_ext = None


    def format_af_data(self):
        # Adapted from structuremap 0.0.10, original function: structuremap.processing.format_alphafold_data
        # Modifications: custom cif path and protein nameing
        protein_number = 0
        alphafold_annotation_l = []
        for transcript_id in tqdm.tqdm(self.structure_reference["present"], desc="Converting structural data: "):
            cif_path, pae_path = self.structure_reference["present"][transcript_id]["path"]   
            protein_number += 1
            structure = Bio.PDB.MMCIF2Dict.MMCIF2Dict(cif_path)
            # boltz-1 cif files lack _atom_site.pdbx_sifts_xref_db_acc, so need to check before accessing it 
            if '_atom_site.pdbx_sifts_xref_db_acc' in structure:
                protein_id = [transcript_id] * len(structure['_atom_site.pdbx_sifts_xref_db_acc'])
                df = pd.DataFrame({'protein_id': protein_id,
                                    'protein_number': protein_number,
                                    'AA': structure['_atom_site.pdbx_sifts_xref_db_res'],
                                    'position': structure['_atom_site.label_seq_id'],
                                    'quality': structure['_atom_site.B_iso_or_equiv'],
                                    'atom_id': structure['_atom_site.label_atom_id'],
                                    'x_coord': structure['_atom_site.Cartn_x'],
                                    'y_coord': structure['_atom_site.Cartn_y'],
                                    'z_coord': structure['_atom_site.Cartn_z']})
            else:
                protein_id = [transcript_id] * len(structure['_atom_site.auth_comp_id'])
                one_letter_aa_codes = [seq1(aa) for aa in structure['_atom_site.auth_comp_id']]
                df = pd.DataFrame({'protein_id': protein_id,
                                    'protein_number': protein_number,
                                    'AA': one_letter_aa_codes,
                                    'position': structure['_atom_site.label_seq_id'],
                                    'quality': structure['_atom_site.B_iso_or_equiv'],
                                    'atom_id': structure['_atom_site.label_atom_id'],
                                    'x_coord': structure['_atom_site.Cartn_x'],
                                    'y_coord': structure['_atom_site.Cartn_y'],
                                    'z_coord': structure['_atom_site.Cartn_z']})
            
            df = df[df.atom_id.isin(['CA', 'CB', 'C', 'N'])].reset_index(drop=True)
            df = df.pivot(index=['protein_id',
                                    'protein_number',
                                    'AA', 'position',
                                    'quality'],
                            columns="atom_id")
            df = pd.DataFrame(df.to_records())

            df = df.rename(columns={"('x_coord', 'CA')": "x_coord_ca",
                                    "('y_coord', 'CA')": "y_coord_ca",
                                    "('z_coord', 'CA')": "z_coord_ca",
                                    "('x_coord', 'CB')": "x_coord_cb",
                                    "('y_coord', 'CB')": "y_coord_cb",
                                    "('z_coord', 'CB')": "z_coord_cb",
                                    "('x_coord', 'C')": "x_coord_c",
                                    "('y_coord', 'C')": "y_coord_c",
                                    "('z_coord', 'C')": "z_coord_c",
                                    "('x_coord', 'N')": "x_coord_n",
                                    "('y_coord', 'N')": "y_coord_n",
                                    "('z_coord', 'N')": "z_coord_n"})

            df = df.apply(safe_to_numeric)
           
            df['secondary_structure'] = 'unstructured'

            if '_struct_conf.conf_type_id' in structure.keys():
                start_idx = [int(i) for i in structure['_struct_conf.beg_label_seq_id']]
                end_idx = [int(i) for i in structure['_struct_conf.end_label_seq_id']]
                note = structure['_struct_conf.conf_type_id']

                for i in np.arange(0, len(start_idx)):
                    df['secondary_structure'] = np.where(
                        df['position'].between(
                            start_idx[i],
                            end_idx[i]),
                        note[i],
                        df['secondary_structure'])

            alphafold_annotation_l.append(df)

        alphafold_annotation = pd.concat(alphafold_annotation_l)
        alphafold_annotation = alphafold_annotation.sort_values(
            by=['protein_number', 'position']).reset_index(drop=True)

        alphafold_annotation['structure_group'] = [re.sub('_.*', '', i)
                                                for i in alphafold_annotation[
                                                'secondary_structure']]
        str_oh = pd.get_dummies(alphafold_annotation['structure_group'],
                                dtype='int64')
        alphafold_annotation = alphafold_annotation.join(str_oh)

        self.alphafold_annotation = alphafold_annotation

        return alphafold_annotation

    
    def annotate_af_accessibility(self) -> pd.DataFrame:
        """
        Annotate accessibility of AAs using pPSE values.

        Returns:
            pd.DataFrame: DataFrame containing accessibility annotations.
        """
        if self.alphafold_annotation is None:
            raise ValueError("AlphaFold data must be formatted before annotation.")
        
        logger.info("Annotate Accessibility")
        full_sphere_exposure, self.full_sphere_exposure_name = annotate_accessibility(
            df=self.alphafold_annotation,
            max_dist=self.annotation_config["full_sphere_exposure_max_dist"],
            max_angle=self.annotation_config["full_sphere_exposure_max_angle"],
            structure_reference=self.structure_reference["present"]
        )
        
        self.alphafold_accessibility = self.alphafold_annotation.merge(
            full_sphere_exposure, how='left', on=['protein_id', 'AA', 'position']
        )

        part_sphere_exposure, self.part_sphere_exposure_name = annotate_accessibility(
            df=self.alphafold_annotation, 
            max_dist=self.annotation_config["half_sphere_exposure_max_dist"], 
            max_angle=self.annotation_config["half_sphere_exposure_max_angle"], 
            structure_reference=self.structure_reference["present"]
        )

        self.alphafold_accessibility = self.alphafold_accessibility.merge(
            part_sphere_exposure, how='left', on=['protein_id','AA','position']
        )

        low_high_hse_thresh = self.annotation_config["low_high_hse_thresh"]
        self.alphafold_accessibility[f'high_acc_{low_high_hse_thresh}'] = np.where(
            self.alphafold_accessibility[self.part_sphere_exposure_name] <= low_high_hse_thresh, 1, 0)
        self.alphafold_accessibility[f'low_acc_{low_high_hse_thresh}'] = np.where(
            self.alphafold_accessibility[self.part_sphere_exposure_name] > low_high_hse_thresh, 1, 0)
        
        return self.alphafold_accessibility
    

    def annoatte_idrs(self) -> pd.DataFrame:
        """
        Annotate Intrinsically Disordered Regions (IDRs).

        Returns:
            pd.DataFrame: DataFrame with IDR annotations.
        """
        if self.alphafold_accessibility is None:
            raise ValueError("AlphaFold Accessibility must be created before annotation.")
        logger.info("Annotate IDRs")
        alphafold_accessibility_smooth = get_smooth_score(
            self.alphafold_accessibility, 
            np.array([self.full_sphere_exposure_name]), 
            [self.annotation_config["smooth_window_size"]]
        )
        attribute_name = f"{self.full_sphere_exposure_name}_smooth{self.annotation_config['smooth_window_size']}"

        alphafold_accessibility_smooth['IDR'] = np.where(
            alphafold_accessibility_smooth[attribute_name]<=self.annotation_config["idr_thresh"], 1, 0
        )

        alphafold_accessibility_smooth_pattern = annotate_proteins_with_idr_pattern(
            alphafold_accessibility_smooth,
            min_structured_length = self.annotation_config["min_structured_length"], 
            max_unstructured_length = self.annotation_config["max_unstructured_length"]
        )

        alphafold_accessibility_smooth_pattern_ext = get_extended_flexible_pattern(
            alphafold_accessibility_smooth_pattern, 
            ['flexible_pattern'], [5]
        )

        self.alphafold_accessibility_smooth_pattern_ext = alphafold_accessibility_smooth_pattern_ext
        
        return alphafold_accessibility_smooth_pattern_ext
    

    def annotate_ptms(self) -> pd.DataFrame:
        """
        Annotate Post-Translational Modifications (PTMs).

        Returns:
            pd.DataFrame: DataFrame with PTM annotations.
        """
        ptm_df = pd.read_csv(os.path.join(self.ptm_dir, "ptms.csv"))
        alphafold_ptms = self.alphafold_accessibility_smooth_pattern_ext.merge(ptm_df, how='left', on=['protein_id', 'AA', 'position'])

        successful_merges = alphafold_ptms[ptm_df.columns.difference(['protein_id', 'AA', 'position'])].notna().any(axis=1).sum()

        alphafold_ptms = alphafold_ptms.fillna(0)

        non_successful_merges = ptm_df.merge(
            self.alphafold_accessibility_smooth_pattern_ext, 
            how='left', on=['protein_id', 'AA', 'position'], 
            indicator=True
        ).query('_merge == "left_only"').drop(columns=['_merge'])[['protein_id', 'AA', 'position']]
        logger.info("Non successful merges:")
        for index, row in non_successful_merges.iterrows():
            logger.info(f"      {row['protein_id']} {row['AA']} {row['position']}")

        logger.info(f"Number of successful merges: {successful_merges} out of {len(ptm_df)}")

        ptm_site_dict = {'p':['S','T','Y'],
                        'p_reg':['S','T','Y']}
        
        self.alphafold_ptms = alphafold_ptms
        return alphafold_ptms
    
    def save_alphafold_ptms(self) -> None:
        """
        Save AlphaFold PTM annotations to a CSV file.
        """
        if isinstance(self.alphafold_ptms, pd.DataFrame):
            logger.info(f"Saveing alphafold_ptms dataframe to {self.protein_structure_ptms_csv}")
            self.alphafold_ptms.to_csv(self.protein_structure_ptms_csv, sep=';', encoding='utf-8')
        else:
            raise ValueError("alphafold_ptms must be created before saving to file.")


