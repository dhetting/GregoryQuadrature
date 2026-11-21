"""Behavior tests for the NumPy port, including an independent base-R oracle."""

import subprocess
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np
import pytest

from gregory_quadrature import Gregory_weights


def rational_weights(n_nodes, order):
    """Exact triangular recurrence, with both endpoint corrections accumulated."""
    coefficients = []
    for row in range(order - 1):
        coefficients.append(
            Fraction(1, row + 2)
            - sum(coefficients[col] / (row - col + 1) for col in range(row))
        )
    weights = [Fraction(1)] * n_nodes
    for col in range(order - 1):
        correction = sum(
            coefficients[row] * (-1) ** col * comb(row, col)
            for row in range(col, order - 1)
        )
        weights[col] -= correction
        weights[n_nodes - 1 - col] -= correction
    return np.array([float(value) for value in weights])


def test_published_eleven_node_order_eight_weights():
    expected = [
        0.304224537037037,
        1.4603835978836,
        0.453463955026455,
        1.47142857142857,
        0.728025793650793,
        1.16494708994709,
        0.728025793650793,
        1.47142857142857,
        0.453463955026455,
        1.4603835978836,
        0.304224537037037,
    ]
    result = Gregory_weights(11, 1, 8)
    assert result.shape == (11,)
    assert result.dtype == np.float64
    np.testing.assert_allclose(result, expected, rtol=0, atol=5e-15)


@pytest.mark.parametrize("n_nodes", [2, 3, 8, 31])
@pytest.mark.parametrize("h", [-0.25, 0, 0.25])
def test_order_two_is_trapezoid(n_nodes, h):
    expected = np.full(n_nodes, h, dtype=float)
    expected[[0, -1]] /= 2
    np.testing.assert_array_equal(Gregory_weights(n_nodes, h, 2), expected)


@pytest.mark.parametrize("order", range(2, 13))
def test_overlapping_endpoint_corrections_against_exact_rationals(order):
    for n_nodes in [order, order + 1, 2 * order + 3]:
        expected = rational_weights(n_nodes, order)
        actual = Gregory_weights(n_nodes, 1, order)
        np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(actual, actual[::-1], rtol=0, atol=2e-15)
        assert actual.sum() == pytest.approx(n_nodes - 1, abs=2e-13)


@pytest.mark.parametrize("order", range(2, 9))
def test_polynomial_exactness_and_first_nonexact_degree(order):
    # This endpoint formula has degree 2*floor(order/2)-1, including odd
    # accuracy orders: the actual R order-3 rule is not quadratic-exact.
    exact_degree = 2 * (order // 2) - 1
    for n_nodes in [order, 21]:
        x = np.linspace(0, 1, n_nodes)
        weights = Gregory_weights(n_nodes, 1 / (n_nodes - 1), order)
        for degree in range(exact_degree + 1):
            assert weights @ x**degree == pytest.approx(
                1 / (degree + 1), rel=0, abs=3e-15
            )
        # Also check the leading nonzero polynomial error against an
        # independently constructed exact-rational rule (not just a bound).
        degree = exact_degree + 1
        rational = rational_weights(n_nodes, order) / (n_nodes - 1)
        expected_error = rational @ x**degree - 1 / (degree + 1)
        assert abs(expected_error) > 1e-14
        assert weights @ x**degree - 1 / (degree + 1) == pytest.approx(
            expected_error, rel=0, abs=3e-15
        )


def test_smooth_integral_published_error_and_fourth_order_convergence():
    exact = np.exp(1) - np.exp(-1)
    x = np.linspace(-1, 1, 11)
    error = Gregory_weights(11, 0.2, 8) @ np.exp(x) - exact
    assert error == pytest.approx(2.3506e-8, abs=1e-12)
    errors = []
    for n_nodes in [17, 33, 65]:
        x = np.linspace(-1, 1, n_nodes)
        errors.append(
            abs(Gregory_weights(n_nodes, 2 / (n_nodes - 1), 4) @ np.exp(x) - exact)
        )
    assert errors[0] / errors[1] > 12
    assert errors[1] / errors[2] > 14


def test_signed_spacing_and_numpy_integer_scalars():
    positive = Gregory_weights(np.int64(11), np.float64(0.2), np.int32(8))
    np.testing.assert_array_equal(Gregory_weights(11, -0.2, 8), -positive)
    x = np.linspace(1, -1, 11)
    assert Gregory_weights(11, -0.2, 8) @ np.exp(x) == pytest.approx(
        np.exp(-1) - np.exp(1), abs=3e-8
    )


@pytest.mark.parametrize("n_nodes,order", [(1, 2), (0, 2), (4, 1), (4, 5)])
def test_invalid_node_order_range(n_nodes, order):
    with pytest.raises(ValueError, match="2 <= order <= n_nodes"):
        Gregory_weights(n_nodes, 1, order)


@pytest.mark.parametrize("bad", [True, 3.5, 3.0, "3", None, np.nan, np.inf, [3]])
@pytest.mark.parametrize("argument", ["n_nodes", "order"])
def test_counts_require_integer_scalars(bad, argument):
    arguments = {"n_nodes": 8, "h": 1, "order": 3}
    arguments[argument] = bad
    with pytest.raises(TypeError, match="integer"):
        Gregory_weights(**arguments)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_spacing_must_be_finite(bad):
    with pytest.raises(ValueError, match="finite"):
        Gregory_weights(8, bad, 3)


@pytest.mark.parametrize("bad", [True, "1", None, 1 + 0j, [1]])
def test_spacing_requires_real_scalar(bad):
    with pytest.raises(TypeError, match="real"):
        Gregory_weights(8, bad, 3)


def test_overflow_does_not_return_nonfinite_weights():
    with pytest.raises(ValueError, match="finite"):
        Gregory_weights(11, np.finfo(float).max, 8)


def test_actual_r_function_with_source_equivalent_base_r_primitives():
    """Not installed-pracma parity: pracma is unavailable and is not installed.

    Execute the unchanged package function body with only its two qualified
    primitive names substituted. Pascal k=1 follows the primary pracma source:
    https://raw.githubusercontent.com/cran/pracma/master/R/pascal.R
    """
    root = Path(__file__).resolve().parents[1]
    script = r"""
    Toeplitz <- function(column, row) {
      outer(seq_along(column), seq_along(row), Vectorize(function(i,j) {
        if (i >= j) column[i-j+1] else row[j-i+1]
      }))
    }
    pascal <- function(n, k) {
      stopifnot(k == 1)
      p <- matrix(0, n, n)
      for (i in seq_len(n)) for (j in seq_len(i))
        p[i,j] <- (-1)^(j-1) * choose(i-1,j-1)
      p
    }
    text <- readLines("R/gregory.R")
    text <- gsub("pracma::Toeplitz", "Toeplitz", text, fixed=TRUE)
    text <- gsub("pracma::pascal", "pascal", text, fixed=TRUE)
    eval(parse(text=text))
    for (order in 2:10) for (n in c(order,order+1,2*order+3)) {
      for (h in c(-0.2,0,1)) {
        cat(n, order, h, sprintf("%.17g", Gregory_weights(n,h,order)), "\n")
      }
    }
    """
    result = subprocess.run(
        ["Rscript", "--vanilla", "-e", script],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    rows = result.stdout.splitlines()
    assert len(rows) == 81
    for row in rows:
        values = row.split()
        n_nodes, order = map(int, values[:2])
        h = float(values[2])
        np.testing.assert_allclose(
            Gregory_weights(n_nodes, h, order),
            np.array(values[3:], dtype=float),
            rtol=2e-13,
            atol=2e-13,
        )
