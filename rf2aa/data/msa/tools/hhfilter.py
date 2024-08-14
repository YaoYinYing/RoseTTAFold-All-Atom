import os
import subprocess
from typing import Optional
from absl import logging


class HHFilter:
    """Python wrapper for the HHfilter binary."""

    def __init__(self,
                 *,
                 binary_path: str,
                 maxseq: int = 100_000,
                 id_threshold: int = 90,
                 cov: int = 75):
        """Initializes the Python HHfilter wrapper.

        Args:
          binary_path: The path to the HHfilter executable.
          maxseq: Maximum number of sequences to retain after filtering.
          id_threshold: Sequence identity threshold for filtering (in percentage).
          cov: Minimum coverage percentage required for sequences.

        Raises:
          ValueError: If the binary path is invalid.
        """
        if not os.path.isfile(binary_path):
            logging.error('HHfilter binary not found at %s', binary_path)
            raise ValueError(f'HHfilter binary not found at {binary_path}')

        self.binary_path = binary_path
        self.maxseq = maxseq
        self.id_threshold = id_threshold
        self.cov = cov

    def filter(self, input_a3m_path: str, output_a3m_path: str) -> None:
        """Executes the HHfilter command.

        Args:
          input_a3m_path: Path to the input A3M file.
          output_a3m_path: Path where the output A3M file will be saved.

        Raises:
          RuntimeError: If HHfilter execution fails.
        """

        cmd = [
            self.binary_path,
            '-i', input_a3m_path,
            '-o', output_a3m_path,
            '-maxseq', str(self.maxseq),
            '-id', str(self.id_threshold),
            '-cov', str(self.cov)
        ]
        os.makedirs(os.path.dirname(output_a3m_path), exist_ok=True)

        logging.info('Launching HHfilter with command: %s', ' '.join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        stdout, stderr = process.communicate()
        retcode = process.wait()

        if retcode:
            logging.error('HHfilter failed. HHfilter stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('HHfilter stderr end')
            raise RuntimeError('HHfilter failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        logging.info('HHfilter successfully completed.')
        return output_a3m_path

