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

"""Library to run HHblits from Python."""

import glob
import os
import subprocess
from typing import Any, List, Mapping, Optional, Sequence

from absl import logging
from rf2aa.data.msa.tools import utils
# Internal import (7716).


_HHBLITS_DEFAULT_P = 20
_HHBLITS_DEFAULT_Z = 500


class HHBlits:
  """Python wrapper of the HHblits binary."""

  def __init__(self,
               *,
               binary_path: str,
               databases: Sequence[str],
               n_cpu: int = 4,
               maxmem: int= 16,
               n_iter: int = 3,
               e_value: float = 0.001,
               mact: float= 0.35,
               neffmax: int= 20,
               cov: int= 25,
               maxseq: int = 1_000_000,
               realign_max: int = 100_000,
               maxfilt: int = 100_000,
               nodiff: bool=True,
               p: int = _HHBLITS_DEFAULT_P,
               z: int = _HHBLITS_DEFAULT_Z):
    """Initializes the Python HHblits wrapper.

    Args:
      binary_path: The path to the HHblits executable.
      databases: A sequence of HHblits database paths. This should be the
        common prefix for the database files (i.e. up to but not including
        _hhm.ffindex etc.)
      n_cpu: The number of CPUs to give HHblits.
      n_iter: The number of HHblits iterations.
      e_value: The E-value, see HHblits docs for more details.
      maxseq: The maximum number of rows in an input alignment. Note that this
        parameter is only supported in HHBlits version 3.1 and higher.
      realign_max: Max number of HMM-HMM hits to realign. HHblits default: 500.
      maxfilt: Max number of hits allowed to pass the 2nd prefilter.
        HHblits default: 20000.
      min_prefilter_hits: Min number of hits to pass prefilter.
        HHblits default: 100.
      all_seqs: Return all sequences in the MSA / Do not filter the result MSA.
        HHblits default: False.
      alt: Show up to this many alternative alignments.
      p: Minimum Prob for a hit to be included in the output hhr file.
        HHblits default: 20.
      z: Hard cap on number of hits reported in the hhr file.
        HHblits default: 500. NB: The relevant HHblits flag is -Z not -z.

    Raises:
      RuntimeError: If HHblits binary not found within the path.
    """
    self.binary_path = binary_path
    self.databases = databases

    for database_path in self.databases:
      if not glob.glob(database_path + '_*'):
        logging.error('Could not find HHBlits database %s', database_path)
        raise ValueError(f'Could not find HHBlits database {database_path}')

    self.n_cpu = n_cpu
    self.maxmem=maxmem
    self.n_iter = n_iter
    self.e_value = e_value
    self.maxseq = maxseq
    self.mact=mact
    self.neffmax=neffmax
    self.cov = cov
    self.nodiff = nodiff
    self.realign_max = realign_max
    self.maxfilt = maxfilt
    self.p = p
    self.z = z

  def query(self, input_fasta_path: str, save_dir: str='hhblits', output_prefix: str='output') -> str:
    """Queries the database using HHblits."""

    os.makedirs(save_dir, exist_ok=True)
    
    a3m_path = os.path.join(save_dir, f'{output_prefix}.a3m')
    if os.path.isfile(a3m_path):
      logging.warning('Found existing a3m file %s', a3m_path)
      return a3m_path

    db_cmd = []
    for db_path in self.databases:
      db_cmd.append('-d')
      db_cmd.append(db_path)
    cmd = [
        self.binary_path,
        '-i', input_fasta_path,
        '-cpu', str(self.n_cpu),
        '-maxmem', str(self.maxmem),
        '-mact', str(self.mact),
        '-neffmax', str(self.neffmax),
        '-cov', str(self.cov),
        '-oa3m', a3m_path,
        '-o', '/dev/null',
        '-n', str(self.n_iter),
        '-e', str(self.e_value),
        '-maxseq', str(self.maxseq),
        '-realign_max', str(self.realign_max),
        '-maxfilt', str(self.maxfilt),
    ]
    
    if self.nodiff:
      cmd += ['-nodiff']
    if self.p != _HHBLITS_DEFAULT_P:
      cmd += ['-p', str(self.p)]
    if self.z != _HHBLITS_DEFAULT_Z:
      cmd += ['-Z', str(self.z)]
    cmd += db_cmd

    logging.info('Launching subprocess "%s"', ' '.join(cmd))
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  
    stdout, stderr = process.communicate()
    retcode = process.wait()

    if retcode:
      # Logs have a 15k character limit, so log HHblits error line by line.
      logging.error('HHblits failed. HHblits stderr begin:')
      for error_line in stderr.decode('utf-8').splitlines():
        if error_line.strip():
          logging.error(error_line.strip())
      logging.error('HHblits stderr end')
      raise RuntimeError('HHblits failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
          stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))


  
    return a3m_path
