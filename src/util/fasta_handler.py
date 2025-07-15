from src.logger import logger

class FastaReader:
    def __init__(self, file_path):
        """
        Initializes the FastaReader with the given file path.
        :param file_path: Path to the FASTA file
        """
        self.file_path = file_path
        self.data = []  # List to store (header, sequence) tuples

    def read_fasta(self):
        """
        Reads the FASTA file and stores headers and sequences in a list.
        """
        try:
            with open(self.file_path, 'r') as file:
                header = None
                sequence_parts = []

                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save the previous sequence before starting a new one
                        if header is not None:
                            self.data.append((header, ''.join(sequence_parts)))
                        header = line[1:]  # Remove the '>'
                        sequence_parts = []  # Reset sequence parts for the new entry
                    else:
                        sequence_parts.append(line)

                # Add the last sequence after the loop
                if header is not None:
                    self.data.append((header, ''.join(sequence_parts)))
        except FileNotFoundError:
            logger.error(f"Error: File '{self.file_path}' not found.")
        except Exception as e:
            logger.error(f"An error occurred: {e}")

    def get_data(self):
        """
        Returns the list of (header, sequence) tuples.
        :return: List of tuples [(header, sequence), ...]
        """
        return self.data
    
    def get_sequence_to_header_dict(self):
        """
        Returns the list of (header, sequence) tuples.
        :return: List of tuples [(header, sequence), ...]
        """
        sequence_to_header_dict = {}
        for header, sequence in self.data:
            sequence_to_header_dict[sequence] = header
        return sequence_to_header_dict


class FastaWriter:
    def __init__(self, file_path, data):
        """
        Initializes the FastaWriter with the given file path and data.
        :param file_path: Path to the FASTA file to write to
        :param data: List of (header, sequence) tuples
        """
        self.file_path = file_path
        self.data = data  # List to store (header, sequence) tuples

    def write_fasta(self):
        """
        Writes the stored (header, sequence) tuples to a FASTA file.
        """
        try:
            with open(self.file_path, 'w') as file:
                for header, sequence in self.data:
                    file.write(f'>{header}\n')
                    file.write(f'{sequence}\n')
        except Exception as e:
            print(f"An error occurred while writing the file: {e}")

