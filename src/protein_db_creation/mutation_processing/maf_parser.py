import os
import glob
from tqdm import tqdm

from src.logger import logger

class MAFParser:
    def __init__(self, dirpath: str):
        """
        Initialize the MAFParser with required parameter dirpath. 
        
        :param dirpath: The path to the directory that contain the maf files.
        """
        self.dirpath = dirpath
        self.variants = {}
        
        if not os.path.exists(self.dirpath):
            try:
                os.makedirs(self.dirpath)
                logger.info(f"Directory {self.dirpath} created. Copy MAF files into dictionary.")
            except Exception as e:
                logger.error(f"Error creating directory: {e}")
        
        self.parse_all_maf_files()

    def parse_all_maf_files(self) -> None:
        """
        Reads all maf files in the specified directory.
        """
        maf_txt_files = glob.glob(os.path.join(self.dirpath, "*.maf.txt"))
        maf_files = glob.glob(os.path.join(self.dirpath, "*.maf"))

        maf_files.extend(maf_txt_files)
        
        if not maf_files:
            logger.warning(f"No MAF files found in directory: {self.dirpath}. No mutations will be translated.")
            return
        
        for maf_file in tqdm(maf_files, desc="Processing MAF-Files"):
            logger.debug(f"Processing file: {maf_file}")
            self.parse(maf_file)
    
    def parse(self, filepath: str) -> None:
        """
        Parses a single MAF file and stores some of its content in the variants dictionary. The dictionary maps a unique variant key to 
        inforamtion about the variant.
        
        :param filepath: Path to the MAF file to be parsed
        """
        try:
            headers = None
            with open(filepath, 'r') as file:
                for line in file:
                    fields = line.strip().split('\t')

                    if not headers:
                        headers = fields
                        continue

                    if len(fields) != len(headers):
                        raise ValueError(f"Invalid MAF format on line: {line}")

                    variant_data = dict(zip(headers, fields))

                    if variant_data['Variant_Type'] != "SNP" or variant_data['Variant_Classification'] == "Silent":
                        continue
                    
                    variant_key = f"{variant_data['Chromosome']}:{variant_data['Start_position']}_{variant_data['Reference_Allele']}/{variant_data['Tumor_Seq_Allele1']}/{variant_data['Tumor_Seq_Allele2']}"
                    self.variants[variant_key] = {
                        'Chromosome': variant_data['Chromosome'],
                        'PositionGenomic': variant_data['Start_position'],
                        'REF': variant_data['Reference_Allele'],
                        'ALT_1': variant_data['Tumor_Seq_Allele1'], 
                        'ALT_2': variant_data['Tumor_Seq_Allele2'],
                    }

        except FileNotFoundError:
            logger.error(f"Error: File not found at {filepath}")
        except Exception as e:
            logger.error(f"Error parsing VCF file {filepath}: {e}")

    def get_variants(self) -> dict:
        """
        Returns the dictionary containing parsed MAF data.

        :return: Dictionary of variants
        """
        return self.variants