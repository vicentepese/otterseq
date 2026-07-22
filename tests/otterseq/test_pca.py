"""Unit test for `OtterPCA`."""

import os
from typing import Any
from unittest.mock import patch

import pandas as pd
import plotly.graph_objects as go
import polars as pl
import pytest
from beartype.roar import BeartypeCallHintParamViolation
from pytest_mock import MockerFixture

from otterseq.pca import OtterPCA


@pytest.mark.parametrize(
    argnames=("args", "error"),
    argvalues=[
        ({"filepath": 100}, BeartypeCallHintParamViolation),
        (
            {"filepath": "filepath", "outpath": 100},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filepath": "filepath", "exclude_hla": "exclude_hla"},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filepath": "filepath", "n_pcs": "n_pcs"},
            BeartypeCallHintParamViolation,
        ),
        ({"filepath": "filepath"}, FileNotFoundError),
        ({"filepath": "tests/data_qc/toy", "n_pcs": -10}, ValueError),
    ],
)
def test_pca_errors(  # noqa: D103
    otter_pca: OtterPCA, args: dict[str, Any], error: Exception
) -> None:
    with pytest.raises(error):
        otter_pca.pca(**args)


@pytest.mark.parametrize(
    argnames=("n_pcs", "expected_n_pcs"),
    argvalues=[
        (2, 2),
        (5, 5),
        (10, 5),
    ],
)
def test_pca(
    otter_pca: OtterPCA, filename: str, n_pcs: int, expected_n_pcs: int
) -> None:
    """Test `OtterPCA.pca`."""
    with patch("subprocess.run") as mock_subprocess:
        # Mock the subprocess call to avoid running the actual PCA script
        mock_subprocess.return_value.returncode = 0

        otter_pca.pca(filepath=filename, exclude_hla=False, n_pcs=n_pcs)

        # Verify subprocess was called with correct arguments
        mock_subprocess.assert_called_once()
        call_args = mock_subprocess.call_args[0][0]
        assert call_args[0] == "bash"
        assert "--pcs" in call_args
        assert str(n_pcs) in call_args


def test_pca_with_pheno_and_exclude_indvs(
    otter_pca: OtterPCA, filename: str
) -> None:
    """Test that `pca` forwards --pheno and --remove and cleans up temp file."""
    exclude_indvs = pd.DataFrame({"FID": ["FAM001"], "IID": ["1"]})
    with patch("subprocess.run") as mock_subprocess:
        mock_subprocess.return_value.returncode = 0

        otter_pca.pca(
            filepath=filename,
            pheno="tests/data_pca/toy.pheno.tsv",
            exclude_hla=False,
            exclude_indvs=exclude_indvs,
        )

        call_args = mock_subprocess.call_args[0][0]
        assert "--pheno" in call_args
        assert (
            call_args[call_args.index("--pheno") + 1]
            == "tests/data_pca/toy.pheno.tsv"
        )
        assert "--remove" in call_args

    # The temporary remove file should be cleaned up after the call
    assert not os.path.isfile(filename + ".pca_remove.tsv")


def test_pca_plot(otter_pca: OtterPCA, filename_pca: str) -> None:
    """."""
    fig = otter_pca.plot_pca(filename=filename_pca, plot=False)

    assert isinstance(
        fig, go.Figure
    ), f"Output not of type go.Figure. Got{type(fig)}"


@pytest.mark.parametrize(
    argnames=("filename", "error"),
    argvalues=[
        ("tests/data_pca/nonexistent", FileExistsError),
        ("tests/data_pca/toy_no_fam", FileExistsError),
        ("tests/data_pca/toy_no_eigenvec", FileExistsError),
    ],
)
def test_match_case_controls_errors(
    otter_pca: OtterPCA,
    filename: str,
    error: FileNotFoundError,
) -> None:
    """Test error handling in match_case_controls."""
    with pytest.raises(error):
        otter_pca.match_case_controls(filename=filename, n_controls=1)


@pytest.mark.parametrize(
    argnames=("n_controls", "unique_controls", "expected_controls"),
    argvalues=[
        (1, False, 2),
        (2, False, 4),
        (1, True, 3),
        (2, True, 6),
    ],
)
def test_match_case_controls(
    otter_pca: OtterPCA,
    filename_pca: str,
    n_controls: int,
    unique_controls: bool,
    expected_controls: int,
) -> None:
    """Test case-control matching based on PCA results."""
    matched_df = otter_pca.match_case_controls(
        filename=filename_pca,
        n_controls=n_controls,
        unique_controls=unique_controls,
    )

    assert isinstance(
        matched_df, pl.DataFrame
    ), "Output should be a Polars DataFrame"

    # Check if the number of controls matches the expected value
    actual_controls = len(matched_df.filter(pl.col("pheno") == 1))
    assert (
        actual_controls == expected_controls
    ), f"Expected {expected_controls} controls, but got {actual_controls}"

    # Check if all cases are included
    fam_df = pl.read_csv(
        filename_pca + ".fam",
        has_header=False,
        separator=" ",
        new_columns=[
            "fid",
            "iid",
            "father_id",
            "mother_id",
            "sex",
            "pheno",
        ],
    )
    total_cases = len(fam_df.filter(pl.col("pheno") == 2))
    matched_cases = len(matched_df.filter(pl.col("pheno") == 2))
    assert (
        matched_cases == total_cases
    ), f"Expected {total_cases} cases, but got {matched_cases}"

    # Check if the correct columns are present
    expected_columns = [
        "fid",
        "iid",
        "pc1",
        "pc2",
        "pc3",
        "pc4",
        "pc5",
        "fid_iid",
        "pheno",
    ]
    assert all(
        col in matched_df.columns for col in expected_columns
    ), "Missing expected columns in the output DataFrame"


def test_match_case_controls_unique(
    otter_pca: OtterPCA, filename_pca: str
) -> None:
    """Test that unique_controls=True results in non-repeated controls."""
    matched_df = otter_pca.match_case_controls(
        filename=filename_pca, n_controls=2, unique_controls=True
    )

    controls = matched_df.filter(pl.col("pheno") == 1)
    unique_controls = controls.unique(subset=["fid_iid"])
    assert len(controls) == len(
        unique_controls
    ), "Controls should not be repeated when unique_controls=True"


@pytest.mark.parametrize(
    argnames=("plot", "expected_fig_show"),
    argvalues=[
        (True, True),
        (False, False),
    ],
)
def test_plot_pca(
    otter_pca: OtterPCA,
    filename_pca: str,
    plot: bool,
    expected_fig_show: bool,
    mocker: MockerFixture,
) -> None:
    """Test the plot_pca method with different plot options."""
    mock_show = mocker.patch("plotly.graph_objects.Figure.show")

    fig = otter_pca.plot_pca(filename=filename_pca, plot=plot)

    assert isinstance(fig, go.Figure)
    assert mock_show.called == expected_fig_show

    # Check if the figure contains the correct data
    all_points = set()  # type: ignore[var-annotated]
    for trace in fig.data:
        all_points.update(zip(trace.x, trace.y, strict=False))

    expected_points = {
        (0.490007, -0.678481),
        (-0.459736, -0.135041),
        (0.0437409, 0.525825),
        (0.158496, 0.0651892),
        (-0.130201, -0.099469),
        (-0.0493007, 0.0327187),
        (-0.459736, -0.135041),
        (0.285393, -0.0345281),
        (-0.257097, 0.000248211),
        (0.378435, 0.458578),
    }

    assert len(all_points) == len(
        expected_points
    ), "Number of points doesn't match"

    for point in all_points:
        assert any(
            pytest.approx(point, abs=1e-6) == expected_point
            for expected_point in expected_points
        ), f"Point {point} not found in expected data"

    # Check if the figure has the correct title
    assert (
        fig.layout.title.text
        == f"Principal Component Analysis of {filename_pca}.eigenvec"
    )

    # Check if the figure has the correct labels
    assert fig.layout.xaxis.title.text == "PC1"
    assert fig.layout.yaxis.title.text == "PC2"


def test_plot_pca_with_fam(otter_pca: OtterPCA, filename_pca: str) -> None:
    """Test plot_pca method for color and symbol information."""
    fig = otter_pca.plot_pca(filename=filename_pca)

    # Check if the figure contains color and symbol information
    assert hasattr(fig.data[0].marker, "color"), "Color information is missing"
    assert hasattr(
        fig.data[0].marker, "symbol"
    ), "Symbol information is missing"


def test_plot_pca_file_not_found(otter_pca: OtterPCA) -> None:  # noqa: D103
    with pytest.raises(FileNotFoundError):
        otter_pca.plot_pca(filename="nonexistent_file")
