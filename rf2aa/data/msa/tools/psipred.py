import os
import subprocess
from typing import Optional, Union
from absl import logging

from dataclasses import dataclass

from rf2aa.data.msa.tools import utils


@dataclass(frozen=True)
class PsipredHformat:
    aa: str
    conf: str
    pred: str


def ParsePsipredHformat(content: str) -> PsipredHformat:
    c=content.splitlines()
    
    resn=[line.strip('  AA: ') for line in c if line.startswith('  AA:')]
    conf=[line.strip('Conf: ') for line in c if line.startswith('Conf:')]
    pred=[line.strip('Pred: ') for line in c if line.startswith('Pred:')]

    return PsipredHformat(''.join(resn), ''.join(conf), ''.join(pred))



class SecondaryStructure:
    """Python wrapper for secondary structure prediction using PSIPRED and CS-BLAST."""

    def __init__(self,
                 *,
                 csblast_dir: str,
                 psipred_data_dir: str,
                 makemat_binary: str,
                 psipred_binary: str,
                 psipass2_binary: str,
                 blast_path: str):
        """Initializes the SecondaryStructure wrapper.

        Args:
          csblast_dir: Directory where CS-BLAST is installed.
          psipred_data_dir: Directory where PSIPRED data files are located.
          makemat_binary: Path to the makemat executable.
          psipred_binary: Path to the psipred executable.
          psipass2_binary: Path to the psipass2 executable.

        Raises:
          ValueError: If any provided path does not exist.
        """
        if not os.path.isdir(csblast_dir):
            raise ValueError(f"CS-BLAST directory not found at {csblast_dir}")
        csbuild_bin=os.path.join(csblast_dir,'bin/csbuild')
        if not os.path.isfile(csbuild_bin):
            raise FileNotFoundError(f'{csbuild_bin} not found')
        if not os.path.isdir(psipred_data_dir):
            raise ValueError(f"PSIPRED data directory not found at {psipred_data_dir}")
        if not os.path.isfile(makemat_binary):
            raise ValueError(f"makemat binary not found at {makemat_binary}")
        if not os.path.isfile(psipred_binary):
            raise ValueError(f"psipred binary not found at {psipred_binary}")
        if not os.path.isfile(psipass2_binary):
            raise ValueError(f"psipass2 binary not found at {psipass2_binary}")

        self.csblast_dir = csblast_dir
        self.psipred_data_dir = psipred_data_dir
        self.makemat_binary = makemat_binary
        self.psipred_binary = psipred_binary
        self.psipass2_binary = psipass2_binary
        self.blast_path=blast_path


    def predict(self, input_a3m_path: str, output_ss_path: str) -> None:
        """Executes secondary structure prediction.

        Args:
          input_a3m_path: Path to the input A3M file.
          output_ss_path: Path where the output secondary structure file will be saved.

        Raises:
          RuntimeError: If any of the subprocesses fail.
        """
        # Create temporary identifiers and filenames
        with utils.tmpdir_manager() as tmp_dir:

            input_a3m_path=os.path.abspath(input_a3m_path)
            output_ss_path=os.path.abspath(output_ss_path)


            temp_id = os.path.basename(input_a3m_path).replace('.a3m', '.tmp')
            temp_fasta = f"{temp_id}.fasta"
            temp_chk = f"{temp_id}.chk"
            temp_pn = f"{temp_id}.pn"
            temp_sn = f"{temp_id}.sn"
            temp_mtx = f"{temp_id}.mtx"
            temp_ss = f"{temp_id}.ss"
            temp_horiz = f"{temp_id}.horiz"

            # Run csbuild to create a checkpoint file
            csbuild_cmd = [
                os.path.join(self.csblast_dir, "bin", "csbuild"),
                "-i", input_a3m_path,
                "-I", "a3m",
                "-D", os.path.join(self.csblast_dir, "data", "K4000.crf"),
                "-o", os.path.join(tmp_dir,temp_chk),
                "-O", "chk"
            ]
            

            
            logging.info('Running csbuild with command: %s', ' '.join(csbuild_cmd))
            # command list will fail on this
            self._run_command(' '.join(csbuild_cmd), "csbuild", env=os.environ.copy(), stdout_redirect=os.path.join(tmp_dir,'csbuild.log'))
            

            # Prepare FASTA file and input files for makemat
            with open(input_a3m_path, 'r') as infile:
                with open(os.path.join(tmp_dir,temp_fasta), 'w') as outfile:
                    outfile.writelines(infile.readlines()[:2])

            with open(os.path.join(tmp_dir,temp_pn), 'w') as pnfile:
                pnfile.write(f"{temp_chk}\n")

            with open(os.path.join(tmp_dir,temp_sn), 'w') as snfile:
                snfile.write(f"{temp_fasta}\n")

            # Run makemat to create a matrix file
            makemat_cmd = [self.makemat_binary, "-P", temp_id]
            logging.info('Running makemat with command: %s', ' '.join(makemat_cmd))
            # command list will fail on this
            self._run_command(' '.join(makemat_cmd), "makemat", cwd=tmp_dir)

            # Run psipred to predict secondary structure
            psipred_cmd = [
                self.psipred_binary,
                temp_mtx,
                os.path.join(self.psipred_data_dir, "weights.dat"),
                os.path.join(self.psipred_data_dir, "weights.dat2"),
                os.path.join(self.psipred_data_dir, "weights.dat3")
            ]

            logging.info('Running psipred with command: %s', ' '.join(psipred_cmd))
            # command list will fail on this
            self._run_command(' '.join(psipred_cmd), "psipred", cwd=tmp_dir, stdout_redirect=os.path.join(tmp_dir,temp_ss))

            # Run psipass2 for refinement
            psipass2_cmd = [
                self.psipass2_binary,
                os.path.join(self.psipred_data_dir, "weights_p2.dat"),
                "1", "1.0", "1.0",
                os.path.join(tmp_dir,f"{temp_id}.csb.hhblits.ss2"),
                os.path.join(tmp_dir,temp_ss),
            ]

            logging.info('Running psipass2 with command: %s', ' '.join(psipass2_cmd))
            # command list will fail on this
            self._run_command(' '.join(psipass2_cmd), "psipass2", cwd=tmp_dir, stdout_redirect=os.path.join(tmp_dir,temp_horiz))

            psipred_res=ParsePsipredHformat(open(os.path.join(tmp_dir,temp_horiz),'r').read())

            # Extract secondary structure prediction
            self._format_output(psipred_res, output_ss_path)

    def _run_command(self, cmd:Union[str, list,tuple], name, cwd:Optional[str]=None, env:Optional[dict[str,str]]=None, back_to_curdir: bool=True, stdout_redirect: Optional[str]=None):
        """Executes a shell command and handles errors.

        Args:
          cmd: Command to execute.
          name: Name of the process for logging.

        Raises:
          RuntimeError: If the command execution fails.
        """
        if cwd is not None:
            ori_wd = os.getcwd()
        process = subprocess.Popen(cmd, shell=True, encoding="utf-8",stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd, env=env)
        stdout, stderr = process.communicate()
        retcode = process.wait()
        if isinstance(stdout_redirect, str):
            open(stdout_redirect,'w').write(stdout)
        if cwd is not None and back_to_curdir:
            os.chdir(ori_wd)

        if retcode:
            logging.error('%s failed. stderr:', name)
            for error_line in stderr.splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            raise RuntimeError(f'{name} failed\nstdout:\n{stdout}\nstderr:\n{stderr[:500_000]}')

    def _format_output(self, psipred_res:PsipredHformat, output_file):
        """Formats the output secondary structure prediction.

        Args:
          psipred_res: PsipredHformat containing raw prediction data.
          output_file: Output file to save formatted secondary structure.

        """
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w') as out:
            out.write(f">ss_pred\n{psipred_res.pred}\n>ss_conf\n{psipred_res.conf}\n")
