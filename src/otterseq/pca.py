"""Principal Component Analysis (PCA) module for genetic data.

This module provides the OtterPCA class for computing and visualizing PCA
results on PLINK binary files. It includes functionality for reading .fam
and .eigenvec files, performing PCA, plotting results, and matching
case-control samples based on PCA output.
"""

import logging
import os
import subprocess
from typing import ClassVar, Literal

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import polars as pl
from beartype import beartype
from scipy.spatial import distance_matrix

import ottersh
from otterseq.utils import read_fam_file, read_pca_file


class OtterPCA:
    """Class used to compute PCA and plot results."""

    _OTTER_SH_PATH = ottersh.__path__[0]
    _PCA_SCRIPT = os.path.join(_OTTER_SH_PATH, "pca.sh")
    _SUFFIXES: ClassVar[list[str]] = [".bim", ".bed", ".fam"]
    _PCA_SUFFIX: Literal[".eigenvec"] = ".eigenvec"

    @beartype
    def pca(
        self,
        filepath: str,
        pheno: str | None = None,
        outpath: str | None = None,
        exclude_hla: bool = True,
        n_pcs: int = 20,
        exclude_indvs: pd.DataFrame | None = None,
    ) -> None:
        """Compute Principal Component Analysis (PCA) on a PLINK binary file.

        Args:
            filepath (str): Path to the binary file, without suffix (e.g. data/toy)
            pheno (str | None, optional): Path to the phenotype file (FID IID
                pheno, space-separated). Defaults to None.
            outpath (str | None, optional): Path to the output file. If None,
                is the same as `filepath`. Defaults to None.
            exclude_hla (bool, optional): True to exclude HLA region from the PCA.
                Defaults to True.
            n_pcs (int, optional): Number of Principal Components. Defaults to 20.
            exclude_indvs (pd.DataFrame | None, optional): FID/IID of individuals
                to exclude from the PCA. Defaults to None.

        Raises:
            FileNotFoundError: If any of the PLINK binary files are missing.
            ValueError: If the number of PCs `n_pcs` is lesser than 0.
        """
        for suf in self._SUFFIXES:
            if not os.path.isfile(filepath + suf):
                raise FileNotFoundError(f"{filepath} not found.")
        if n_pcs <= 0:
            raise ValueError(
                f"Number of PCs cannot be lower or equal to 0. Got {n_pcs}"
            )

        outpath = outpath or filepath

        remove_path = None
        if exclude_indvs is not None:
            remove_path = filepath + ".pca_remove.tsv"
            exclude_indvs.to_csv(
                remove_path, sep="\t", header=False, index=False
            )

        command = [
            "bash",
            self._PCA_SCRIPT,
            "--bfile",
            filepath,
            "--outpath",
            outpath,
            "--exclude-hla",
            str(exclude_hla),
            "--pcs",
            str(n_pcs),
        ]
        if pheno is not None:
            command.extend(["--pheno", pheno])
        if remove_path is not None:
            command.extend(["--remove", remove_path])
        subprocess.run(  # noqa: S603
            command, text=True, capture_output=True, check=False
        )
        if remove_path is not None:
            os.remove(remove_path)

    @beartype
    def plot_pca(self, filename: str, plot: bool = True) -> go.Figure:
        """Plot PCA results.

        Args:
            filename (str | os.PathLike[str]): Path to the folder with files,
                without suffix (e.g, "data/toy")
            plot (bool, optional): If True, plots the PCA. Defaults to True.

        Raises:
            FileNotFoundError: If `filename`.eigenvec does not exist.

        Returns:
            go.Figure: plotly Figure which can be plotted with `.show()`
        """
        pca_filepath = filename + self._PCA_SUFFIX
        if not os.path.isfile(pca_filepath):
            raise FileNotFoundError(f"{pca_filepath} was not found.")

        fam_filepath = filename + ".fam"
        if not os.path.isfile(fam_filepath):
            logging.warning(
                FileNotFoundError(
                    f"{fam_filepath} file not found. Phenotypes will not be plotted"
                )
            )
            fam_df = None
        else:
            fam_df = read_fam_file(fam_filepath, vars_to_string=True)

        pca_df = read_pca_file(pca_filepath)

        if fam_df is not None:
            pca_df = pca_df.join(
                other=fam_df[["fid", "iid", "sex_str", "pheno_str"]],
                on=["fid", "iid"],
                how="left",
            )
            fig = px.scatter(
                data_frame=pca_df,
                x="pc1",
                y="pc2",
                color="pheno_str",
                symbol="sex_str",
                hover_data=["fid", "iid"],
                title=f"Principal Component Analysis of {pca_filepath}",
                labels={
                    "pc1": "PC1",
                    "pc2": "PC2",
                    "fid": "FID",
                    "iid": "IID",
                    "pheno_str": "Phenotype",
                    "sex_str": "Sex",
                },
            )
        else:
            fig = px.scatter(
                data_frame=pca_df,
                x="pc1",
                y="pc2",
                hover_data=["fid", "iid"],
                title=f"Principal Component Analysis of {pca_filepath}",
                labels={
                    "pc1": "PC1",
                    "pc2": "PC2",
                    "fid": "FID",
                    "iid": "IID",
                },
            )
        if plot:
            fig.show()
        return fig

    @beartype
    def match_case_controls(
        self, filename: str, n_controls: int, unique_controls: bool = False
    ) -> pl.DataFrame:
        """Match cases to controls based on PCA results.

        Args:
            filename (str): Path to the base filename (without suffix) containing the .fam and .eigenvec files.
            n_controls (int): Number of controls to match to each case.
            unique_controls (bool, optional): If True, ensures that controls are not repeated across cases. Defaults to False.

        Raises:
            FileExistsError: If either the .fam or .eigenvec file is not found.

        Returns:
            pl.DataFrame: DataFrame containing matched cases and controls.

        Example:
        ```python linenums="1"
            otter_pca = OtterPCA()
            matched_df = otter_pca.match_case_controls(filename="data/toy", n_controls=2, unique_controls=True)
            print(matched_df)
        ```
        ```text
            shape: (6, 10)
            ┌─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────────┬─────────┐
            │ fid │ iid │ pc0 │ pc1 │ pc2 │ pc3 │ pc4 │ pheno │ fid_iid │ pheno_str │
            ├─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────────┼─────────┤
            │ 1   │ 1   │ ... │ ... │ ... │ ... │ ... │ 2    │ 1-1     │ case    │
            │ 2   │ 2   │ ... │ ... │ ... │ ... │ ... │ 1    │ 2-2     │ control │
            │ 3   │ 3   │ ... │ ... │ ... │ ... │ ... │ 1    │ 3-3     │ control │
            │ 4   │ 4   │ ... │ ... │ ... │ ... │ ... │ 2    │ 4-4     │ case    │
            │ 5   │ 5   │ ... │ ... │ ... │ ... │ ... │ 1    │ 5-5     │ control │
            │ 6   │ 6   │ ... │ ... │ ... │ ... │ ... │ 1    │ 6-6     │ control │
            └─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────────┴─────────┘
        ```
        """
        fam_file = filename + ".fam"
        pca_file = filename + self._PCA_SUFFIX

        if not os.path.isfile(fam_file):
            raise FileExistsError(
                f"{fam_file} was not found. match_pca requires both .fam and .eigenvec files."
            )
        if not os.path.isfile(pca_file):
            raise FileExistsError(
                f"{pca_file} was not found. match_pca requires both .fam and .eigenvec files."
            )

        pca_df = read_pca_file(filename=filename)
        fam_df = read_fam_file(filename=filename, vars_to_string=False)
        merged_df = pca_df.join(
            other=fam_df[["fid", "iid", "pheno"]],
            on=["fid", "iid"],
            how="left",
        )
        merged_df = merged_df.with_columns(
            pl.concat_str(["fid", "iid"], separator="-").alias("fid_iid")
        )

        cases = merged_df.filter(pl.col("pheno") == 2)
        controls = merged_df.filter(pl.col("pheno") == 1)

        pc_cols = [col for col in merged_df.columns if col.startswith("pc")]
        cases_pcs = cases.select(pl.col(pc_cols)).to_numpy()
        controls_pcs = controls.select(pl.col(pc_cols)).to_numpy()

        cases_pcs_dist = distance_matrix(cases_pcs, controls_pcs)

        # Select unique controls
        matched_controls: set[str]
        if unique_controls:
            matched_controls = set()
            for case_dist in cases_pcs_dist:
                control_idx = np.argsort(case_dist).tolist()
                matched = 0
                for idx in control_idx:
                    control_id = controls["fid_iid"][idx]
                    if control_id not in matched_controls:
                        matched_controls.add(control_id)
                        matched += 1
                    if matched == n_controls:
                        break

        # Select n_controls closest controls (cases can share controls)
        else:
            pcs_idx = (
                np.argsort(cases_pcs_dist, axis=1)[:, :n_controls]
                .flatten()
                .tolist()
            )
            matched_controls = set(controls["fid_iid"][pcs_idx])

        matched_controls_list = list(matched_controls)
        merged_df_matched: pl.DataFrame = merged_df.filter(
            (pl.col("fid_iid").is_in(matched_controls_list))
            | (pl.col("pheno") == 2)
        )

        return merged_df_matched
