import os
import subprocess
from typing import Optional
from absl import logging


class SignalP6:
    """Python wrapper for the SignalP 6.0 binary."""

    def __init__(self,
                 *,
                 binary_path: str,
                 organism: str = "other",
                 mode: str = "slow",
                 format: str = "none"):
        """Initializes the Python SignalP6 wrapper.

        Args:
          binary_path: The path to the SignalP 6.0 executable.
          organism: The organism type for prediction ('euk', 'gram+', 'gram-', 'arch', 'other').
          mode: Prediction mode ('fast' or 'slow').
          format: Output format ('short', 'long', 'none').

        Raises:
          ValueError: If the binary path is invalid or parameters are invalid.
        """
        if not os.path.isfile(binary_path):
            logging.error('SignalP6 binary not found at %s', binary_path)
            raise ValueError(f'SignalP6 binary not found at {binary_path}')

        valid_organisms = {"euk", "gram+", "gram-", "arch", "other"}
        if organism not in valid_organisms:
            raise ValueError(f"Invalid organism type: {organism}. Must be one of {valid_organisms}.")

        valid_modes = {"fast", "slow"}
        if mode not in valid_modes:
            raise ValueError(f"Invalid mode: {mode}. Must be 'fast' or 'slow'.")

        valid_formats = {"short", "long", "none"}
        if format not in valid_formats:
            raise ValueError(f"Invalid format: {format}. Must be 'short', 'long', or 'none'.")

        self.binary_path = binary_path
        self.organism = organism
        self.mode = mode
        self.format = format

    def predict(self, input_fasta_path: str, output_dir: str) -> str:
        """Executes the SignalP 6.0 prediction.

        Args:
          input_fasta_path: Path to the input FASTA file.
          output_dir: Directory where the output files will be saved.

        Return: 
            trimmed_output_path: Path to the trimmed output fasta file.
        Raises:
          RuntimeError: If SignalP 6.0 execution fails.
        """
        processed_fasta = os.path.join(output_dir, 'processed_entries.fasta')
        if not os.path.isfile(processed_fasta):

            os.makedirs(output_dir, exist_ok=True)
            cmd = [
                self.binary_path,
                '--fastafile', input_fasta_path,
                '--organism', self.organism,
                '--output_dir', output_dir,
                '--format', self.format,
                '--mode', self.mode
            ]

            logging.info('Launching SignalP6 with command: %s', ' '.join(cmd))
            process = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            stdout, stderr = process.communicate()
            retcode = process.wait()

            if retcode:
                logging.error('SignalP6 failed. SignalP6 stderr begin:')
                for error_line in stderr.decode('utf-8').splitlines():
                    if error_line.strip():
                        logging.error(error_line.strip())
                logging.error('SignalP6 stderr end')
                raise RuntimeError('SignalP6 failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                    stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

            logging.info('SignalP6 successfully completed.')

        # Check if processed entries file is not empty
        return self.trim_decision(
            input_fasta_path=input_fasta_path,
            processed_fasta=processed_fasta
        )
        
        
    def trim_decision(self, input_fasta_path: str, processed_fasta: str) -> str:
        if os.path.getsize(processed_fasta) > 0:
            return processed_fasta
        else:
            return input_fasta_path
        
