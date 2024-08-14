import os
import shutil

from omegaconf import DictConfig
from absl import logging

from rf2aa.data.msa.pipeline import Pipeline


def make_msa(fasta_file: str, chain: str, model_runner: DictConfig):

    CONDA_PREFIX = os.environ.get("CONDA_PREFIX", None)
    if CONDA_PREFIX is None:
        logging.warning('You are runing RoseTTAFold-All-Atom WITHTOUT CONDA environment.')
        logging.warning('Make sure you have configured the correct paths of all binaries and data directories.')

    out_dir_base = os.path.abspath(model_runner.config.output_path)
    hash: str = model_runner.config.job_name
    out_dir = os.path.join(out_dir_base, hash, chain)
    os.makedirs(out_dir, exist_ok=True)

    # sequence databases
    DB_UR30 = model_runner.config.database_params.DB_UR30
    DB_BFD = model_runner.config.database_params.DB_BFD
    num_cpus = model_runner.config.database_params.num_cpus
    ram_gb = model_runner.config.database_params.mem
    template_database = model_runner.config.database_params.DB_PDB100

    # binaries 
    hhblits_binary= model_runner.config.binaries.hhblits
    hhfilter_binary= model_runner.config.binaries.hhfilter
    hhsearch_binary= model_runner.config.binaries.hhsearch
    signalp_binary= model_runner.config.binaries.signalp
    makemat_binary= model_runner.config.binaries.makemat
    psipred_binary= model_runner.config.binaries.psipred
    psipass2_binary= model_runner.config.binaries.psipass2

    # data dirs
    blast_dir= model_runner.config.data_dirs.blast
    psipred_dir= model_runner.config.data_dirs.psipred
    csblast_dir= model_runner.config.data_dirs.csblast

    msa_pipeline = Pipeline(
        hhblits_binary=hhblits_binary if hhblits_binary else shutil.which("hhblits"),
        hhfilter_binary=hhfilter_binary if hhfilter_binary else shutil.which("hhfilter"),
        hhsearch_binary=hhsearch_binary if hhsearch_binary else shutil.which("hhsearch"),
        signalp_binary=signalp_binary if signalp_binary else shutil.which("signalp"),
        makemat_binary=makemat_binary if makemat_binary else shutil.which("makemat"),
        blast_path=blast_dir if blast_dir else os.path.join(CONDA_PREFIX, "share/blast-2.2.26/blast-2.2.26"),
        psipred_binary=psipred_binary if psipred_binary else shutil.which("psipred"),
        psipred_data_dir=psipred_dir if psipred_dir else os.path.join(CONDA_PREFIX, "share", "psipred_4.01/data"),
        psipass2_binary=psipass2_binary if psipass2_binary else shutil.which("psipass2"),
        csblast_dir=csblast_dir if csblast_dir else os.path.join(CONDA_PREFIX, "share/csblast-2.2.3"),
        uniref30_database=DB_UR30,
        bfd_databse=DB_BFD,
        template_database=template_database,
        ncpu=num_cpus,
        max_mem=ram_gb,
        out_prefix="rf2aa_",
        save_dir=out_dir,
    )

    trim_fasta = msa_pipeline.run_signalp(fasta_path=os.path.abspath(fasta_file))

    out_a3m = msa_pipeline.run_msa_search(fasta_path=trim_fasta)
    ss2_file = msa_pipeline.run_psipred(msa_path=out_a3m)
    out_hhr, out_atab = msa_pipeline.run_hhsearch(msa_path=out_a3m, ss2_path=ss2_file)

    for f in [out_a3m, out_hhr, out_atab]:
        if not os.path.isfile(f):
            raise FileNotFoundError(f"{f} not found")

    return out_a3m, out_hhr, out_atab
