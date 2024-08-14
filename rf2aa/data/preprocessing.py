import os
import shutil

from omegaconf import DictConfig

from rf2aa.data.msa.pipeline import Pipeline


def make_msa(fasta_file: str, chain: str, model_runner: DictConfig):

    CONDA_PREFIX = os.environ.get("CONDA_PREFIX", None)
    if CONDA_PREFIX is None:
        raise NotImplementedError("Make sure conda is activated")

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

    msa_pipeline = Pipeline(
        hhblits_binary=shutil.which("hhblits"),
        hhfilter_binary=shutil.which("hhfilter"),
        hhsearch_binary=shutil.which("hhsearch"),
        signalp_binary=shutil.which("signalp"),
        makemat_binary=shutil.which("makemat"),
        blast_path=os.path.join(CONDA_PREFIX, "share/blast-2.2.26/blast-2.2.26"),
        psipred_binary=shutil.which("psipred"),
        psipred_data_dir=os.path.join(CONDA_PREFIX, "share", "psipred_4.01/data"),
        psipass2_binary=shutil.which("psipass2"),
        csblast_dir=os.path.join(CONDA_PREFIX, "share/csblast-2.2.3"),
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
