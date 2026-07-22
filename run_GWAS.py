#!/usr/bin/env python
"""Performs GWAS: Merge, Quality Control, PCA, Case matching, Logistic regression.

- Merge: from multiple binary .ped files, merge files using only common SNPs
- Quality Control (QC): filter based on minimum allele frequency, missing genotype
    and phenotype. Removes IBD-related individuals, duplicated variants and IIDs/FIDs.
- PCA: Principal Component Analysis - verifies that there is no batch biases.
- Case Matching: Matches cases with controls based on a predefined ratio,
    using PCs and Euclidean distance.
- Logistic Regression: Computes a logistic regression controlling for PCs.
"""

import glob
import logging
import os
import sys
from typing import Any

import pandas as pd
import plotly.express as px
import polars as pl
import yaml

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
)

from otterseq.operations import OtterOps
from otterseq.pca import OtterPCA
from otterseq.qc import OtterQC
from otterseq.snp import OtterSNP


def _prepare_plink_files(
    snp: OtterSNP,
    gwas_dir: str,
    gwas_bin_dir: str,
    gwas_out_dir: str,
    prefix: str,
) -> str:
    if not os.listdir(gwas_bin_dir):
        logging.info("Binarizing files")
        snp.binarize_files(gwas_dir, gwas_bin_dir)
    else:
        logging.info("Files already binarized")

    bed_files = glob.glob(os.path.join(gwas_bin_dir, "*.bed"))
    if len(bed_files) == 1:
        return os.path.join(
            gwas_bin_dir, os.path.splitext(os.path.basename(bed_files[0]))[0]
        )

    snp.merge_files(
        gwas_bin_dir, outpath=gwas_out_dir, prefix=prefix, only_common=True
    )
    return os.path.join(gwas_out_dir, prefix)


def _get_not_in_pheno(fam_path: str, pheno_path: str) -> pd.DataFrame:
    """Return FID/IID rows from .fam that are not present in the pheno file."""
    fam = pd.read_csv(
        fam_path,
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["FID", "IID"],
        dtype=str,
    )
    pheno = pd.read_csv(
        pheno_path,
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["FID", "IID"],
        dtype=str,
    )
    pheno_ids = set(zip(pheno["FID"], pheno["IID"], strict=False))
    mask = [
        tuple(row) not in pheno_ids
        for row in fam[["FID", "IID"]].itertuples(index=False)
    ]
    return fam[mask]


def _run_qc(
    qc_runner: OtterQC,
    plink_prefix: str,
    qc_outpath: str,
    settings: dict[str, Any],
) -> str:
    pheno = settings["file"]["pheno"]
    not_in_pheno = _get_not_in_pheno(plink_prefix + ".fam", pheno)
    logging.info(f"Individuals not in pheno file: {len(not_in_pheno)}")

    cutoff_file = plink_prefix + ".king.cutoff.out.id"
    if os.path.isfile(cutoff_file):
        logging.info(
            "IBD cutoff file already exists, skipping IBD computation"
        )
        ibd_excluded = pd.read_csv(
            cutoff_file, sep=r"\s+", header=0, names=["FID", "IID"], dtype=str
        )
    else:
        logging.info("Computing IBD")
        ibd_excluded = qc_runner.ibd(
            plink_prefix,
            pheno=pheno,
            threshold=settings["IBD_threshold"],
            exclude_indvs=not_in_pheno if not not_in_pheno.empty else None,
        )
    logging.info(f"Individuals excluded by IBD: {len(ibd_excluded)}")

    qc_prefix = os.path.join(qc_outpath, os.path.basename(plink_prefix))
    if os.path.isfile(qc_prefix + ".bed"):
        logging.info("QC output already exists, skipping QC")
        return qc_prefix

    dup_vars = qc_runner.get_duplicate_vars(plink_prefix)
    dup_rsids = qc_runner.get_duplicate_rsids(plink_prefix)
    exclude_vars = list(set(dup_vars + dup_rsids)) or None

    dup_indvs = qc_runner.extract_duplicate_individuals(plink_prefix + ".fam")
    all_excluded_indvs = pd.concat(
        [ibd_excluded, dup_indvs, not_in_pheno]
    ).drop_duplicates()

    qc_runner.qc(
        plink_prefix,
        pheno=pheno,
        outpath=qc_outpath,
        exclude_vars=exclude_vars,
        exclude_indvs=(
            all_excluded_indvs if not all_excluded_indvs.empty else None
        ),
        maf=settings["maf"],
        geno_miss=settings["genomiss"],
        indv_miss=settings["phenomiss"],
    )

    # qc.sh writes output to {outpath}/{basename(bfile)}
    return qc_prefix


def _run_pca_and_match(
    pca_runner: OtterPCA, qc_prefix: str, settings: dict[str, Any]
) -> None:
    pheno = settings["file"]["pheno"]
    # Output PCA to the same path as QC files so that .eigenvec and .fam
    # are co-located, which is required by plot_pca() and match_case_controls().
    if os.path.isfile(qc_prefix + ".eigenvec"):
        logging.info("PCA output already exists, skipping PCA computation")
    else:
        logging.info("Computing PCA")
        not_in_pheno = _get_not_in_pheno(qc_prefix + ".fam", pheno)
        pca_runner.pca(
            qc_prefix,
            pheno=pheno,
            outpath=qc_prefix,
            exclude_hla=True,
            exclude_indvs=not_in_pheno if not not_in_pheno.empty else None,
        )

    logging.info("Plotting PCA (batch)")
    fig_batch = pca_runner.plot_pca(qc_prefix, plot=True)
    fig_batch.write_html("GWAS_batch.html")

    logging.info("Matching cases and controls")
    matched_df = pca_runner.match_case_controls(
        qc_prefix,
        n_controls=settings["ControlCaseRatio"],
        unique_controls=False,
    )

    matched_df.select(["fid", "iid", "pheno"]).write_csv(
        settings["file"]["pheno_matched"], separator=" ", include_header=False
    )
    matched_df.filter(pl.col("pheno") == 1).select(["fid", "iid"]).write_csv(
        settings["file"]["matched_controls"], include_header=True
    )

    logging.info("Plotting PCA (matched)")
    n_cases = matched_df.filter(pl.col("pheno") == 2).shape[0]
    n_controls_matched = matched_df.filter(pl.col("pheno") == 1).shape[0]
    fig_matched = px.scatter(
        data_frame=matched_df,
        x="pc1",
        y="pc2",
        color=matched_df["pheno"].cast(pl.Utf8).to_list(),
        hover_data=["fid", "iid"],
        title=f"PCA (Matched) — Cases: {n_cases} // Controls: {n_controls_matched}",
        labels={"pc1": "PC1", "pc2": "PC2", "fid": "FID", "iid": "IID"},
    )
    fig_matched.show()
    fig_matched.write_html("GWAS_matched.html")


def main() -> None:
    """Run the full GWAS pipeline."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    with open("settings.yaml") as f:
        settings = yaml.safe_load(f)

    gwas_dir = settings["directory"]["GWAS"]
    gwas_bin_dir = settings["directory"]["GWAS_binaries"]
    gwas_out_dir = settings["directory"]["GWAS_out"]
    qc_outpath = settings["plinkFiles"]["GWASQC"]
    prefix = settings["plinkFiles"]["prefix"]

    snp = OtterSNP()
    qc_runner = OtterQC()
    pca_runner = OtterPCA()
    ops = OtterOps()

    plink_prefix = _prepare_plink_files(
        snp, gwas_dir, gwas_bin_dir, gwas_out_dir, prefix
    )
    qc_prefix = _run_qc(qc_runner, plink_prefix, qc_outpath, settings)
    _run_pca_and_match(pca_runner, qc_prefix, settings)

    logging.info("Computing logistic regression")
    ops.run_logistic_regression(qc_prefix, outpath=qc_prefix)
    logging.info("Logistic regression complete")


if __name__ == "__main__":
    main()
