import os
import glob
from tqdm import tqdm

from src.logger import logger

class VCFParser:
    def __init__(self, dirpath):
        """
        Initializes the VCFParser with the given directory path. If the directory doesn't
        exist, it creates it. It also collects and processes all VCF files in the directory.

        :param dirpath: Path to the directory containing VCF files
        """
        self.dirpath = dirpath
        self.variants = {}  # Dictionary to store parsed data
        
        # Check if the directory exists, create it if not
        if not os.path.exists(self.dirpath):
            try:
                os.makedirs(self.dirpath)
                logger.info(f"Directory {self.dirpath} created. Copy VCF files into dictionary.")
            except Exception as e:
                logger.error(f"Error creating directory: {e}")
        
        # Process all VCF files in the directory
        self.parse_all_vcf_files()

    def parse_all_vcf_files(self):
        """
        Processes all VCF files in the specified directory.
        """
        vcf_files = glob.glob(os.path.join(self.dirpath, "*.vcf"))
        
        if not vcf_files:
            logger.warning(f"No VCF files found in directory: {self.dirpath}. No mutations will be translated.")
            return

        for vcf_file in tqdm(vcf_files, desc="Processing VCF-Files"):
            logger.debug(f"Processing file: {vcf_file}")
            self.parse(vcf_file)
    
    def parse(self, filepath):
        """
        Parses a single VCF file and stores its content in the variants dictionary.
        
        :param filepath: Path to the VCF file to be parsed
        """
        try:
            with open(filepath, 'r') as file:
                for line in file:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')

                    if len(fields) < 8:
                        raise ValueError(f"Invalid VCF format on line: {line}")

                    chrom = fields[0] 
                    pos = fields[1]  
                    id_ = fields[2] 
                    ref = fields[3] 
                    alt = fields[4]  
                    qual = fields[5]  
                    filter_ = fields[6]  
                    info = fields[7]

                    if "," in alt:
                        alt1, alt2 = alt.split(",")
                    else: 
                        alt1 = ref # dirty fix to make it comparable to maf files
                        alt2 = alt 

                    variant_key = f"{chrom}:{pos}_{ref}/{alt1}/{alt2}"
                    self.variants[variant_key] = {
                        'ID': id_,
                        'REF': ref,
                        'ALT_1': alt1,
                        'ALT_2': alt2, 
                        'QUAL': qual,
                        'FILTER': filter_,
                        'INFO': info
                    }
        except FileNotFoundError:
            logger.error(f"Error: File not found at {filepath}")
        except Exception as e:
            logger.error(f"Error parsing VCF file {filepath}: {e}")

    def get_variants(self):
        """
        Returns the dictionary containing parsed VCF data.

        :return: Dictionary of variants
        """
        return self.variants
    