"""OtterQC unit tests."""

import os
import tempfile
from typing import Any
from unittest.mock import patch

import pandas as pd
import pytest
from beartype.roar import BeartypeCallHintParamViolation

from otterseq.qc import OtterQC


def test_otter_qc_instantiation() -> None:
    """Test that OtterQC can be instantiated without __init__ method."""
    qc = OtterQC()
    assert isinstance(qc, OtterQC)
    # Test that class variables are accessible
    assert hasattr(qc, "_OTTER_SH_PATH")
    assert hasattr(qc, "_IBD_SCRIPT")
    assert hasattr(qc, "_DUP_VARS")
    assert hasattr(qc, "_DUP_RSID")
    assert hasattr(qc, "_QC_SCRIPT")
    assert hasattr(qc, "_SUFFIXES")


def test_qc_command_building() -> None:
    """Test that the improved command building works correctly."""
    qc = OtterQC()

    # Test with all optional arguments
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test files
        test_file = os.path.join(temp_dir, "test")
        for suffix in [".bed", ".bim", ".fam"]:
            with open(test_file + suffix, "w") as f:
                f.write("test")

        # Create exclude files
        exclude_vars = ["rs0001", "rs0002"]
        exclude_indvs = pd.DataFrame({"FID": ["FAM001"], "IID": ["1"]})

        # Mock subprocess.run to capture the command
        with patch("subprocess.run") as mock_run:
            qc.qc(
                filename=test_file,
                outpath=temp_dir,
                exclude_vars=exclude_vars,
                exclude_indvs=exclude_indvs,
                maf=0.01,
                geno_miss=0.1,
                indv_miss=0.05,
            )

            # Check that subprocess.run was called
            assert mock_run.called

            # Get the command that was passed to subprocess.run
            call_args = mock_run.call_args
            command = call_args[0][0]

            # Verify base command
            assert command[0] == "bash"
            assert command[1] == qc._QC_SCRIPT
            assert command[2] == "--bfile"
            assert command[3] == test_file
            assert command[4] == "--outpath"
            assert command[5] == temp_dir

            # Verify optional arguments were added
            assert "--indv-miss" in command
            assert "--geno-miss" in command
            assert "--maf" in command
            assert "--rm-vars" in command
            assert "--rm-indv" in command

            # Verify values
            indv_miss_idx = command.index("--indv-miss")
            assert command[indv_miss_idx + 1] == "0.05"

            geno_miss_idx = command.index("--geno-miss")
            assert command[geno_miss_idx + 1] == "0.1"

            maf_idx = command.index("--maf")
            assert command[maf_idx + 1] == "0.01"


def test_qc_command_building_minimal() -> None:
    """Test command building with minimal arguments."""
    qc = OtterQC()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test files
        test_file = os.path.join(temp_dir, "test")
        for suffix in [".bed", ".bim", ".fam"]:
            with open(test_file + suffix, "w") as f:
                f.write("test")

        with patch("subprocess.run") as mock_run:
            qc.qc(filename=test_file)

            call_args = mock_run.call_args
            command = call_args[0][0]

            # Should only have base command, no optional args
            assert (
                len(command) == 6
            )  # bash, script, --bfile, filename, --outpath, outpath
            assert "--indv-miss" not in command
            assert "--geno-miss" not in command
            assert "--maf" not in command
            assert "--rm-vars" not in command
            assert "--rm-indv" not in command


def test_qc_command_building_with_pheno() -> None:
    """Test that qc forwards the --pheno argument when provided."""
    qc = OtterQC()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test files
        test_file = os.path.join(temp_dir, "test")
        for suffix in [".bed", ".bim", ".fam"]:
            with open(test_file + suffix, "w") as f:
                f.write("test")

        with patch("subprocess.run") as mock_run:
            qc.qc(filename=test_file, pheno="pheno.txt", outpath=temp_dir)

            command = mock_run.call_args[0][0]

            assert "--pheno" in command
            assert command[command.index("--pheno") + 1] == "pheno.txt"


def test_qc_command_building_partial() -> None:
    """Test command building with some optional arguments."""
    qc = OtterQC()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test files
        test_file = os.path.join(temp_dir, "test")
        for suffix in [".bed", ".bim", ".fam"]:
            with open(test_file + suffix, "w") as f:
                f.write("test")

        with patch("subprocess.run") as mock_run:
            qc.qc(filename=test_file, maf=0.05, exclude_vars=["rs0001"])

            call_args = mock_run.call_args
            command = call_args[0][0]

            # Should have maf and rm-vars but not others
            assert "--maf" in command
            assert "--rm-vars" in command
            assert "--indv-miss" not in command
            assert "--geno-miss" not in command
            assert "--rm-indv" not in command


@pytest.mark.parametrize(
    argnames=("args", "error"),
    argvalues=[
        ({"filename": 1000}, BeartypeCallHintParamViolation),
        (
            {"filename": "tests", "threshold": "threshold"},
            BeartypeCallHintParamViolation,
        ),
        ({"filename": "tests/toy"}, FileNotFoundError),
        ({"filename": "tests/data/toy", "threshold": -10}, ValueError),
        ({"filename": "tests/data/toy", "threshold": 1}, ValueError),
        ({"filename": "tests/data/toy", "threshold": 10}, ValueError),
    ],
)
def test_ibd_errors(  # noqa: D103
    otter_qc: OtterQC, args: dict[str, Any], error: Exception
) -> None:
    with pytest.raises(error):
        otter_qc.ibd(**args)


@pytest.mark.parametrize(
    argnames=("filename", "threshold", "expected_iid"),
    argvalues=[
        ("tests/data/toy", 0.25, ["7"]),
        ("tests/data/toy", 0.2, ["1", "3", "7", "10"]),
    ],
)
def test_ibd(
    otter_qc: OtterQC, filename: str, threshold: float, expected_iid: list[str]
) -> None:
    """Test Inheritance By Descent (IBD)."""
    indv_out = otter_qc.ibd(filename, threshold)

    assert set(indv_out.IID) == set(
        expected_iid
    ), f"wrong IID for threshold {threshold}. \
        Expected {expected_iid}, got {indv_out.IID}"

    king_file = filename + ".kin0"
    cut_in_file = filename + ".king.cutoff.in.id"
    cut_out_file = filename + ".king.cutoff.out.id"

    assert os.path.isfile(king_file), f"File {king_file} was not created"
    assert os.path.isfile(cut_in_file), f"File {cut_in_file} was not created"
    assert os.path.isfile(cut_out_file), f"File {cut_out_file} was not created"

    os.remove(king_file)
    os.remove(cut_in_file)
    os.remove(cut_out_file)


def test_ibd_with_pheno_and_exclude_indvs() -> None:
    """Test that `ibd` forwards --pheno and --remove and cleans up temp file."""
    qc = OtterQC()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test files
        test_file = os.path.join(temp_dir, "test")
        for suffix in [".bed", ".bim", ".fam"]:
            with open(test_file + suffix, "w") as f:
                f.write("test")
        # ibd reads the KING cutoff output produced by PLINK
        with open(test_file + ".king.cutoff.out.id", "w") as f:
            f.write("FID\tIID\n")

        exclude_indvs = pd.DataFrame({"FID": ["FAM001"], "IID": ["1"]})
        with patch("subprocess.run") as mock_run:
            qc.ibd(
                filename=test_file,
                pheno="pheno.txt",
                threshold=0.2,
                exclude_indvs=exclude_indvs,
            )

            command = mock_run.call_args[0][0]

            assert "--pheno" in command
            assert command[command.index("--pheno") + 1] == "pheno.txt"
            assert "--remove" in command

        # The temporary remove file should be cleaned up after the call
        assert not os.path.isfile(test_file + ".ibd_remove.tsv")


@pytest.mark.parametrize(
    argnames=("filename", "error"),
    argvalues=[
        (1000, BeartypeCallHintParamViolation),
        ("tests/toy_qc", FileNotFoundError),
    ],
)
def test_duplicate_vars_errors(  # noqa: D103
    otter_qc: OtterQC, filename: int | str, error: Exception
) -> None:
    with pytest.raises(error):
        otter_qc.get_duplicate_vars(filename)


def test_duplicate_vars(
    otter_qc: OtterQC, filename_dup: str, duplicated_variants: list[str]
) -> None:
    """Test `OtterQC.get_duplicate_vars`."""
    dup_vars = otter_qc.get_duplicate_vars(filename_dup)

    basename = os.path.splitext(filename_dup)[0]
    path_dup_list = basename + ".dupvar"
    assert os.path.isfile(
        path_dup_list
    ), "Duplicated list of variants not created"
    assert set(dup_vars) == set(
        duplicated_variants
    ), "Duplicated variants do not match expected"

    os.remove(path_dup_list)
    os.remove(basename + ".log")


def test_duplicate_vars_no_output() -> None:
    """Return an empty list when PLINK produces no .dupvar file."""
    qc = OtterQC()

    with tempfile.TemporaryDirectory() as temp_dir:
        test_file = os.path.join(temp_dir, "test.bim")
        with open(test_file, "w") as f:
            f.write("test")

        # Mock subprocess so no .dupvar file is created
        with patch("subprocess.run"):
            assert qc.get_duplicate_vars(test_file) == []


@pytest.mark.parametrize(
    argnames=("filename", "error"),
    argvalues=[
        (1000, BeartypeCallHintParamViolation),
        ("tests/toy_qc", FileNotFoundError),
    ],
)
def test_duplicate_rsid_errors(  # noqa: D103
    otter_qc: OtterQC, filename: int | str, error: Exception
) -> None:
    with pytest.raises(error):
        otter_qc.get_duplicate_rsids(filename)


def test_duplicate_rsid(
    otter_qc: OtterQC, filename_dup: str, duplicated_rsid: list[str]
) -> None:
    """Test `OtterQC.get_duplicate_rsids`."""
    dup_rsids = otter_qc.get_duplicate_rsids(filename_dup)

    basename = os.path.splitext(filename_dup)[0]
    path_dup_list = basename + ".rmdup.mismatch"
    path_rm_dup_list = basename + ".rmdup.list"
    path_log = basename + ".log"
    assert os.path.isfile(
        path_dup_list
    ), "Duplicated list of rsids not created"
    assert set(dup_rsids) == set(
        duplicated_rsid
    ), "Duplicated variants do not match expected"

    os.remove(path_dup_list)
    os.remove(path_rm_dup_list)
    os.remove(path_log)


def test_duplicate_rsid_no_output() -> None:
    """Return an empty list when PLINK produces no .rmdup.mismatch file."""
    qc = OtterQC()

    with tempfile.TemporaryDirectory() as temp_dir:
        test_file = os.path.join(temp_dir, "test.bim")
        with open(test_file, "w") as f:
            f.write("test")

        # Mock subprocess so no .rmdup.mismatch file is created
        with patch("subprocess.run"):
            assert qc.get_duplicate_rsids(test_file) == []


def test_extract_duplicate_individuals_not_found() -> None:
    """Raise FileNotFoundError when the .fam file is missing."""
    qc = OtterQC()
    with pytest.raises(FileNotFoundError):
        qc.extract_duplicate_individuals("tests/data/no_file.fam")


def test_extract_duplicate_individuals() -> None:
    """Return duplicated FID/IID rows from a .fam file."""
    qc = OtterQC()

    with tempfile.TemporaryDirectory() as temp_dir:
        fam_file = os.path.join(temp_dir, "test.fam")
        with open(fam_file, "w") as f:
            f.write("FAM1 1 0 0 1 1\n")
            f.write("FAM1 1 0 0 1 1\n")  # duplicate of the previous row
            f.write("FAM2 2 0 0 1 1\n")

        dup = qc.extract_duplicate_individuals(fam_file)

        assert len(dup) == 1
        assert dup.iloc[0]["FID"] == "FAM1"


@pytest.mark.parametrize(
    argnames=("args", "error"),
    argvalues=[
        ({"filename": 1000}, BeartypeCallHintParamViolation),
        (
            {"filename": "filename", "outpath": 1000},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "exclude_vars": [1000]},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "exclude_indvs": [1000]},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "maf": "0.05"},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "geno_miss": "0.05"},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "indv_miss": "0.05"},
            BeartypeCallHintParamViolation,
        ),
        ({"filename": "tests/data/toy", "maf": 10}, ValueError),
        ({"filename": "tests/data/toy", "maf": 1}, ValueError),
        ({"filename": "tests/data/toy", "maf": 0}, ValueError),
        ({"filename": "tests/data/toy", "maf": -10}, ValueError),
        ({"filename": "tests/data/toy", "geno_miss": 10}, ValueError),
        ({"filename": "tests/data/toy", "geno_miss": -10}, ValueError),
        ({"filename": "tests/data/data", "indv_miss": 10}, ValueError),
        ({"filename": "tests/data/toy", "indv_miss": -10}, ValueError),
        ({"filename": "filename"}, FileNotFoundError),
    ],
)
def test_qc_error(  # noqa: D103
    otter_qc: OtterQC, args: dict[str, Any], error: Exception
) -> None:
    with pytest.raises(error):
        otter_qc.qc(**args)


@pytest.mark.parametrize(
    argnames=(
        "outpath",
        "exclude_vars",
        "exclude_indvs",
        "maf",
        "geno_miss",
        "indv_miss",
        "rm_vars",
        "rm_indv",
    ),
    argvalues=[
        ("test", None, None, None, None, None, [], []),
        ("test", None, None, 0.1, None, None, ["rs0002"], []),
        ("test", None, None, None, 0.1, None, ["rs0007"], []),
        ("test", None, None, None, None, 0.1, [], ["1", "2", "7"]),
        ("test", ["rs0002"], None, None, None, None, ["rs0002"], []),
        (
            "test",
            None,
            pd.DataFrame({"FID": ["FAM001"], "IID": ["1"]}),
            None,
            None,
            None,
            [],
            ["1"],
        ),
        ("test", None, None, 0.1, 0.1, None, ["rs0002", "rs0007"], []),
        ("test", None, None, 0.1, None, 0.1, ["rs0002"], ["1", "2", "7"]),
        ("test", None, None, None, 0.1, 0.1, [], ["1", "2", "7"]),
        ("test", None, None, 0.1, 0.1, 0.1, [], ["1", "2", "7"]),
    ],
)
def test_qc(  # noqa: D103
    otter_qc: OtterQC,
    filename_qc: str,
    outpath: str,
    exclude_vars: list[str] | None,
    exclude_indvs: pd.DataFrame | None,
    maf: float,
    geno_miss: float,
    indv_miss: float,
    rm_indv: list[str],
    rm_vars: list[str],
) -> None:
    """Test Quality Control."""
    os.makedirs(outpath, exist_ok=True)

    otter_qc.qc(
        filename_qc,
        outpath=outpath,
        exclude_vars=exclude_vars,
        exclude_indvs=exclude_indvs,
        maf=maf,
        geno_miss=geno_miss,
        indv_miss=indv_miss,
    )

    # Check that output files were created
    suffixes = [".bed", ".bim", ".fam", ".log"]
    for suffix in suffixes:
        outfile_path = os.path.join(outpath, "toy_qc" + suffix)
        assert os.path.isfile(
            outfile_path
        ), f"{outfile_path} was not correctly created"
        os.remove(outfile_path)

    # Clean up any additional PLINK output files
    additional_files = [".irem", ".nosex"]
    for suffix in additional_files:
        additional_file_path = os.path.join(outpath, "toy_qc" + suffix)
        if os.path.isfile(additional_file_path):
            os.remove(additional_file_path)

    # Check that exclude files were created if provided
    if exclude_vars is not None:
        rm_vars_path = os.path.join(outpath, "toy_qc.rmvars")
        assert os.path.isfile(
            rm_vars_path
        ), "Exclude variants file not created"
        with open(rm_vars_path) as f:
            file_vars = f.read().splitlines()
        assert set(file_vars) == set(
            exclude_vars
        ), "Exclude variants in file don't match expected"
        os.remove(rm_vars_path)

    if exclude_indvs is not None:
        rm_indv_path = os.path.join(outpath, "toy_qc.rmindv")
        assert os.path.isfile(
            rm_indv_path
        ), "Exclude individuals file not created"
        exclude_indvs_file = pd.read_csv(
            rm_indv_path, sep="\t", header=None, names=["FID", "IID"]
        )
        assert set(exclude_indvs_file.IID.astype(str)) == set(
            rm_indv
        ), "Exclude individuals in file don't match expected"
        os.remove(rm_indv_path)

    os.removedirs(outpath)
