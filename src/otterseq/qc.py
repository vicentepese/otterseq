"""otterseq.qc is a module used to compute Quality Control on PLINK files.

Supported functionalities are:
- Inheritance By Descent (IBD)
- Computation of duplicated variants by coordinates and allele codes.
- Computation of duplicated variant rsIDs.
- Quality Control, encompassing MAF, missing genotyping and individual rates,
    exclusion of variants and individuals filtering.
"""

import os
import subprocess
from typing import ClassVar

import pandas as pd
from beartype import beartype

import ottersh


class OtterQC:
    """Class used to conduct Quality Control operations on PLINk files."""

    _OTTER_SH_PATH = ottersh.__path__[0]
    _IBD_SCRIPT = os.path.join(_OTTER_SH_PATH, "ibd.sh")
    _DUP_VARS = os.path.join(_OTTER_SH_PATH, "duplicate_vars.sh")
    _DUP_RSID = os.path.join(_OTTER_SH_PATH, "duplicate_rsid.sh")
    _QC_SCRIPT = os.path.join(_OTTER_SH_PATH, "qc.sh")
    _SUFFIXES: ClassVar[list[str]] = [".bim", ".bed", ".fam"]

    @beartype
    def ibd(
        self,
        filename: str,
        threshold: float | int = 0.25,
        pheno: str | None = None,
        exclude_indvs: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """Compute Identity By Descent (IBD) between individuals.

        This function computes the IBD between individuals using PLINK2. It
        creates a KING table and then applies a threshold to identify related
        individuals.

        Args:
            filename (str): Path to the file with prefix (e.g. `data/toy`,
                where the "data" folder contains a "toy.map", "toy.bed", and
                "toy.bim" file).
            threshold (float | int): Threshold for IBD. Must be between 0 and 1.
            pheno (str | None, optional): Path to the phenotype file (FID IID
                pheno, space-separated). Defaults to None.
            exclude_indvs (pd.DataFrame | None, optional): FID/IID of individuals
                to exclude before computing IBD. Defaults to None.

        Returns:
            pd.DataFrame: DataFrame with the FID and IID of individuals that
                exceed the threshold.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the threshold is not in the range (0, 1).
        """
        # Enforce logic
        if not 0 < threshold < 1:
            raise ValueError(f"threshold not in range (0,1). Got {threshold}")

        for suf in self._SUFFIXES:
            if not os.path.isfile(filename + suf):
                raise FileNotFoundError(f"{filename}{suf} not found")

        remove_path = None
        if exclude_indvs is not None:
            remove_path = filename + ".ibd_remove.tsv"
            exclude_indvs.to_csv(
                remove_path, sep="\t", header=False, index=False
            )

        command = [
            "bash",
            self._IBD_SCRIPT,
            "--bfile",
            filename,
            "--threshold",
            str(threshold),
        ]
        if pheno is not None:
            command.extend(["--pheno", pheno])
        if remove_path is not None:
            command.extend(["--remove", remove_path])
        subprocess.run(command, check=False)  # noqa: S603

        if remove_path is not None:
            os.remove(remove_path)

        # Return filtered out patients
        indv_out = pd.read_csv(
            filename + ".king.cutoff.out.id",
            sep=r"\s+",
            header=0,
            names=["FID", "IID"],
            dtype={"FID": str, "IID": str},
        )

        return indv_out

    @beartype
    def get_duplicate_vars(self, filename: str) -> list[str]:
        """Get the duplicated variants based on coordinates and allele codes.

        For instance, given a `.bim` file as follows:

        ```text
        1	rs0001	0	100100	G	A
        1	rs0001_dup1	0	100100	G	A
        1	rs0001_dup2	0	100100	G	A
        1	rs0002	0	100200	A	T
        1	rs0002_dup1	0	100200	A	T
        1	rs0003	0	100300	C	G
        1	rs0001	0	100300	T	G
        1	rs0004	0	100400	T	C
        1	rs0005	0	100500	T	A
        1	rs0002	0	100500	T	G
        ```

        This function creates the following `.dupvar` file:

        ```text
        CHR	POS	ALLELES	IDS
        1	100100	A,G	rs0001_dup1 rs0001_dup2
        1	100200	A,T	rs0002_dup1
        ```

        Indicating that `rs0001_dup1` and `rs0002_dup1` are duplicates. Then,
        it returns ["rs0001_dup1", "rs0002_dup1"]

        Args:
            filename (str): Path to the .bim file. Accepts the file without
                suffix (e.g, "toy"), or with suffix ("toy.bim")

        Returns:
            list[str]: List with duplicated variants rsIDs.
        """
        if not filename.endswith(".bim"):
            filename += ".bim"
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} not found.")

        basename = os.path.splitext(filename)[0]
        command = ["bash", self._DUP_VARS, "--bfile", basename]
        subprocess.run(  # noqa: S603
            command, capture_output=True, text=True, check=False
        )

        dupvar_path = basename + ".dupvar"
        if not os.path.isfile(dupvar_path):
            return []

        dup_vars_df = pd.read_csv(
            dupvar_path,
            usecols=[3],
            names=["rsid"],
            sep="\t",
            header=0,
        )
        dup_vars: list[str] = dup_vars_df.rsid.str.split().explode().to_list()
        return dup_vars

    @beartype
    def get_duplicate_rsids(self, filename: str) -> list[str]:
        """Get duplicated variants rsIDs.

        This function returns the variants that have duplicated rsIDs. Given the
        following `.bim` file:

        ```text
        1	rs0001	0	100100	G	A
        1	rs0001_dup1	0	100100	G	A
        1	rs0002	0	100200	A	T
        1	rs0002_dup1	0	100200	A	T
        1	rs0003	0	100300	C	G
        1	rs0001	0	100300	T	G
        1	rs0004	0	100400	T	C
        1	rs0005	0	100500	T	A
        1	rs0002	0	100500	T	G
        ```

        The function creates a `.rmdup.list` file with duplicated rsIDs:

        ```text
        rs0001
        rs0002
        ```

        Then, it returns a list ["rs0001", "rs0002"].

        Args:
            filename (str): Path to the .bim file Accepts path to the file with
                suffix ("toy.bim") or without suffix ("toy")

        Returns:
            list[str]: List of duplicated rsIDs.
        """
        if not filename.endswith(".bim"):
            filename += ".bim"
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} not found.")

        basename = os.path.splitext(filename)[0]
        command = ["bash", self._DUP_RSID, "--bfile", basename]
        subprocess.run(  # noqa: S603
            command, capture_output=True, text=True, check=False
        )

        rmdup_path = basename + ".rmdup.mismatch"
        if not os.path.isfile(rmdup_path):
            return []

        dup_rsid_df = pd.read_csv(rmdup_path, usecols=[0], names=["rsid"])
        dup_rsid: list[str] = dup_rsid_df.rsid.to_list()
        return dup_rsid

    @beartype
    def extract_duplicate_individuals(self, filename: str) -> pd.DataFrame:
        """Extract FID and IIDs from duplicated individuals from a .fam file.

        Args:
            filename (str): Path to the .fam file

        Returns:
            pd.DataFrame: DataFrame with the FID and IID of duplicated individuals
        """
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} was not found,")

        fam_data = pd.read_csv(
            filename,
            sep=r"\s+",
            header=None,
            usecols=[0, 1],
            names=["FID", "IID"],
        )
        dup_data = fam_data[fam_data.duplicated()]
        return dup_data

    @beartype
    def qc(  # noqa: C901
        self,
        filename: str,
        pheno: str | None = None,
        outpath: str | None = None,
        exclude_vars: list[str] | None = None,
        exclude_indvs: pd.DataFrame | None = None,
        maf: float | int | None = None,
        geno_miss: float | int | None = None,
        indv_miss: float | int | None = None,
    ) -> None:
        """Quality Control on PLINK files.

        This function runs a quality control on a PLINK file. It can be used to
        exclude variants and individuals based on missing genotyping rates,
        missing individual rates, and minor allele frequency.

        Args:
            filename (str): Path to the file with prefix (e.g. `data/toy`,
                where the "data" folder contains a "toy.map", "toy.bed", and
                "toy.bim" file).
            pheno (str | None, optional): Path to the phenotype file (FID IID
                pheno, space-separated). Defaults to None.
            outpath (str | None): Path to the directory where the output should
                be written. If None, uses `filename`.
            exclude_vars (list[str] | None): List of variants to exclude.
            exclude_indvs (pd.DataFrame | None): DataFrame with the FID and IID
                of individuals to exclude.
            maf (float | int | None): Minor allele frequency threshold.
            geno_miss (float | int | None): Missing genotyping rate threshold.
            indv_miss (float | int | None): Missing individual rate threshold.
        """
        # Enforce logic
        if maf is not None and not 0 < maf < 0.5:
            raise ValueError(f"maf not in range (0,1). Got {maf}")
        if geno_miss is not None and not 0 <= geno_miss <= 1:
            raise ValueError(f"geno_miss not in range (0,1). Got {geno_miss}")
        if indv_miss is not None and not 0 <= indv_miss <= 1:
            raise ValueError(f"pheno_miss not in range (0,1). Got {indv_miss}")

        for suf in self._SUFFIXES:
            if not os.path.isfile(filename + suf):
                raise FileNotFoundError(f"{filename}{suf} not found")

        outpath = outpath or os.path.dirname(filename) or "."

        # Extract basename from input file for exclude files
        basename = os.path.basename(filename)

        if exclude_vars is not None:
            rm_vars_path = os.path.join(outpath, basename + ".rmvars")
            with open(rm_vars_path, "w") as f:
                for var in exclude_vars:
                    f.write(f"{var}\n")

        if exclude_indvs is not None:
            rm_indv_path = os.path.join(outpath, basename + ".rmindv")
            exclude_indvs.to_csv(
                rm_indv_path, sep="\t", header=False, index=False
            )

        # Build command
        command = [
            "bash",
            self._QC_SCRIPT,
            "--bfile",
            filename,
            "--outpath",
            outpath,
        ]
        if pheno is not None:
            command.extend(["--pheno", pheno])

        # Add optional arguments
        optional_args = [
            ("--indv-miss", indv_miss),
            ("--geno-miss", geno_miss),
            ("--maf", maf),
            ("--rm-vars", rm_vars_path if exclude_vars is not None else None),
            ("--rm-indv", rm_indv_path if exclude_indvs is not None else None),
        ]

        # Flatten the optional arguments
        for flag, value in optional_args:
            if value is not None:
                command.extend([flag, str(value)])

        subprocess.run(  # noqa: S603
            command, capture_output=True, text=True, check=False
        )
