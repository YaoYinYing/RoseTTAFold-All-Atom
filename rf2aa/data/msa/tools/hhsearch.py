# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Library to run HHsearch from Python."""

import glob
import os
import subprocess
from typing import Sequence

from absl import logging

from rf2aa.data.msa import parsers
from rf2aa.data.msa.tools import utils
# Internal import (7716).


class HHSearch:
  """Python wrapper of the HHsearch binary."""

  def __init__(self,
               *,
               binary_path: str,
               databases: Sequence[str],
               maxseq: int = 1_000_000,
               mact: float= 0.35,
               min_align:int=50,
               max_align:int=500,
               min_hit:int=50,
               max_hit:int=500,
               n_cpu: int = 4,
               maxmem: int= 16,
               aliw: int=100_000,
               e_value: float = 100,
               p_value: float=5.0,
               verbose: int=0,
               ):
    """Initializes the Python HHsearch wrapper.

    Args:
      binary_path: The path to the HHsearch executable.
      databases: A sequence of HHsearch database paths. This should be the
        common prefix for the database files (i.e. up to but not including
        _hhm.ffindex etc.)
      maxseq: The maximum number of rows in an input alignment. Note that this
        parameter is only supported in HHBlits version 3.1 and higher.

    Raises:
      RuntimeError: If HHsearch binary not found within the path.
    """
    self.binary_path = binary_path
    self.databases = databases
    self.maxseq = maxseq
    self.mact=mact
    self.min_align=min_align
    self.max_align=max_align
    self.min_hit=min_hit
    self.max_hit=max_hit
    self.maxmem=maxmem
    self.aliw=aliw
    self.e_value=e_value
    self.p_value=p_value
    self.n_cpu=n_cpu
    self.verbose=verbose
    for database_path in self.databases:
      if not glob.glob(database_path + '_*'):
        logging.error('Could not find HHsearch database %s', database_path)
        raise ValueError(f'Could not find HHsearch database {database_path}')

  @property
  def output_format(self) -> str:
    return 'hhr'

  @property
  def input_format(self) -> str:
    return 'a3m'

  def query(self, input_a3m_path: str, output_hhr_path: str, output_atab_path: str):
    """Queries the database using HHsearch using a given a3m."""

    os.makedirs(os.path.dirname(output_hhr_path), exist_ok=True)
    os.makedirs(os.path.dirname(output_atab_path), exist_ok=True)

    db_cmd = []
    for db_path in self.databases:
      db_cmd.append('-d')
      db_cmd.append(db_path)
    cmd = [self.binary_path,
            '-i', input_a3m_path,
            '-o', output_hhr_path,
            '-atab', output_atab_path,
            '-maxseq', str(self.maxseq),
            '-cpu',str(self.n_cpu),
            '-b', str(self.min_align),
            '-B', str(self.max_align),
            '-z', str(self.min_hit),
            '-Z', str(self.max_hit),
            '-maxmem', str(self.maxmem),
            '-aliw', str(self.aliw),
            '-e', str(self.e_value),
            '-p', str(self.p_value),
            '-mact', str(self.mact),
            '-v', str(self.verbose),
            ] + db_cmd

    logging.info('Launching subprocess "%s"', ' '.join(cmd))
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  
    stdout, stderr = process.communicate()
    retcode = process.wait()

    if retcode:
      # Stderr is truncated to prevent proto size errors in Beam.
      raise RuntimeError(
          'HHSearch failed:\nstdout:\n%s\n\nstderr:\n%s\n' % (
              stdout.decode('utf-8'), stderr[:100_000].decode('utf-8')))


    return output_hhr_path


  def get_template_hits(self,
                        output_string: str,
                        input_sequence: str) -> Sequence[parsers.TemplateHit]:
    """Gets parsed template hits from the raw string output by the tool."""
    del input_sequence  # Used by hmmseach but not needed for hhsearch.
    return parsers.parse_hhr(output_string)
