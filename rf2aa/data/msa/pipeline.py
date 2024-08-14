import os
from dataclasses import dataclass
from absl import logging


from rf2aa.data.msa.tools import (
    signalp,
    hhblits,
    hhfilter,
    hhsearch,
    psipred,
    utils,
)


@dataclass
class Pipeline:

    # binaries for sequence alignment
    hhblits_binary: str
    hhfilter_binary: str
    hhsearch_binary: str
    signalp_binary: str
    makemat_binary: str

    blast_path: str

    # binaries for secondary structure
    psipred_binary: str
    psipred_data_dir: str
    psipass2_binary: str
    csblast_dir: str

    # databases
    uniref30_database: str
    bfd_databse: str
    template_database: str

    # resources
    ncpu: int = 8
    max_mem: int = 32

    out_prefix: str = "t000_"

    save_dir: str = "rf2aa/features/"

    def __post_init__(self) -> None:
        for k, v in self.__dict__.items():
            if k.endswith("binary") and not os.path.exists(v):
                raise ValueError(f"{k} does not exist: {v}")

        os.environ["BLASTMAT"] = os.path.join(self.blast_path, "data")

    @staticmethod
    def check_a3m_seq_number(a3m_path: str, cutoff) -> bool:
        seq_num = len(
            [
                line
                for line in open(a3m_path, "r").read().split()
                if line.startswith(">")
            ]
        )
        return seq_num > cutoff

    def run_signalp(self, fasta_path):
        with utils.timing("SignalP"):
            save_dir = os.path.join(self.save_dir, "signalp")

            signalp_runner = signalp.SignalP6(
                binary_path=self.signalp_binary,
            )
            updated_fasta = signalp_runner.predict(
                input_fasta_path=fasta_path, output_dir=save_dir
            )

            return updated_fasta

    def run_hhblit(
        self,
        fasta_path: str,
        database: str,
        e_values: tuple[float] = (1e-10, 1e-6, 1e-3),
    ) -> tuple[str, bool]:

        db_alias = os.path.basename(database)

        save_dir = os.path.join(self.save_dir, "hhblits")
        tmp_dir = os.path.join(save_dir, db_alias)

        os.makedirs(tmp_dir, exist_ok=True)

        hhblits_runner = hhblits.HHBlits(
            binary_path=self.hhblits_binary,
            databases=[database],
            n_cpu=self.ncpu,
            maxmem=self.max_mem,
            n_iter=4,
            mact=0.35,
            neffmax=20,
            cov=25,
            maxseq=100_000_000,
            realign_max=100_000_000,
            maxfilt=1_000_000,
            nodiff=True,
        )

        for e_value in e_values:
            hhblits_res_path = os.path.join(tmp_dir, f"{self.out_prefix}.{e_value}.a3m")
            if not os.path.isfile(hhblits_res_path):
                logging.info(
                    f"Running HHblits against {db_alias} with E-value cutoff {e_value}"
                )
                hhblits_runner.e_value = e_value
                with utils.timing(f"hhblits vs {e_value=}"):
                    hhblits_res_path = hhblits_runner.query(
                        input_fasta_path=fasta_path,
                        save_dir=tmp_dir,
                        output_prefix=f"{self.out_prefix}.{e_value}",
                    )

            filtered_msa, passed = self.run_hhfilter(hhblits_res_path)

            if passed:
                return filtered_msa, True

        return filtered_msa, False

    def run_hhfilter(self, input_a3m_path) -> tuple[str, bool]:

        save_dir = os.path.join(self.save_dir, "hhfilter")

        hhfilter_runner = hhfilter.HHFilter(
            binary_path=self.hhfilter_binary,
            maxseq=100_000,
            id_threshold=90,
        )

        for cov_value, max_seq in zip((75, 50), (2_000, 4_000)):
            output_a3m_path = os.path.join(
                save_dir,
                f"{os.path.basename(input_a3m_path)[:-4]}.id90cov{cov_value}.a3m",
            )

            if not os.path.isfile(output_a3m_path):
                hhfilter_runner.cov = cov_value
                with utils.timing(f"hhfilter vs {cov_value=}, {max_seq=}"):
                    hhfilter_runner.filter(
                        input_a3m_path=input_a3m_path, output_a3m_path=output_a3m_path
                    )
            passed = self.check_a3m_seq_number(output_a3m_path, max_seq)
            logging.info(
                f"Coverage check ({os.path.basename(output_a3m_path)}): {passed}"
            )
            if passed:
                return output_a3m_path, True

        return output_a3m_path, False

    def run_msa_search(self, fasta_path: str) -> str:
        dbs = (self.uniref30_database, self.uniref30_database)
        e_value_groups = (
            (
                1e-10,
                1e-6,
                1e-3,
            ),
            1e-3,
        )

        for db, e_values in zip(dbs, e_value_groups):
            msa, finshed = self.run_hhblit(
                fasta_path=fasta_path, database=db, e_values=e_values
            )
            if finshed:
                # early return of MSA search
                return msa
            fasta_path = msa

        # fallback to final msa results
        return msa

    def run_psipred(
        self,
        msa_path: str,
    ) -> str:

        output_ss_path = os.path.join(
            self.save_dir, "psipred", f"{self.out_prefix}.ss2"
        )

        if os.path.isfile(output_ss_path):
            logging.info(f"Found existing secondary structure file {output_ss_path}")
            return output_ss_path

        psipred_runner = psipred.SecondaryStructure(
            csblast_dir=self.csblast_dir,
            psipred_binary=self.psipred_binary,
            psipass2_binary=self.psipass2_binary,
            makemat_binary=self.makemat_binary,
            psipred_data_dir=self.psipred_data_dir,
            blast_path=self.blast_path,
        )

        with utils.timing("psipred"):
            psipred_runner.predict(
                input_a3m_path=msa_path, output_ss_path=output_ss_path
            )

        return output_ss_path

    def run_hhsearch(self, msa_path: str, ss2_path: str) -> tuple[str, str]:
        hhsearch_runner = hhsearch.HHSearch(
            binary_path=self.hhsearch_binary,
            databases=[self.template_database],
            maxseq=100_000,
            n_cpu=self.ncpu,
            maxmem=self.max_mem,
            mact=0.05,
            min_align=50,
            max_align=500,
            min_hit=50,
            max_hit=500,
            e_value=100,
            p_value=5.0,
            verbose=0,
        )
        output_hhr_path = os.path.join(self.save_dir, "hhr", f"{self.out_prefix}.hhr")
        output_atab_path = os.path.join(
            self.save_dir, "atab", f"{self.out_prefix}.atab"
        )
        if os.path.isfile(output_hhr_path) and os.path.isfile(output_atab_path):
            return output_hhr_path, output_atab_path

        with utils.tmpdir_manager() as tmpdir:
            tmp_ss2_msa_path = os.path.join(tmpdir, f"{self.out_prefix}.msa0.ss2.a3m")

            with open(tmp_ss2_msa_path, "w") as ss2_msa_f:
                ss2_msa_f.write(open(ss2_path, "r").read().strip())
                ss2_msa_f.write("\n")
                ss2_msa_f.write(open(msa_path, "r").read().strip())
                ss2_msa_f.write("\n")
                with utils.timing("hhsearch"):
                    hhsearch_runner.query(
                        input_a3m_path=tmp_ss2_msa_path,
                        output_hhr_path=output_hhr_path,
                        output_atab_path=output_atab_path,
                    )
            return output_hhr_path, output_atab_path
