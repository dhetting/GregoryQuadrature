"""Gregory weights for integration on an equally spaced grid.

NumPy translation of GregoryQuadrature's R/gregory.R, by Dylan Hettinger
and Bengt Fornberg. Copyright Dylan Hettinger. Licensed under GPL-3.0-only,
matching the original package's GPL-3 license.

Original algorithm: Bengt Fornberg, "Gregory's quadrature method",
https://www.colorado.edu/amath/sites/default/files/attached-files/gregory.pdf
Signed Pascal convention: pracma::pascal(n, 1), whose zero-based entry
(i, j) is (-1)**j * binomial(i, j) for j <= i:
https://raw.githubusercontent.com/cran/pracma/master/R/pascal.R
"""

from math import comb
from numbers import Real
from operator import index

import numpy as np


def Gregory_weights(n_nodes, h, order):
    """Return a one-dimensional float64 array of Gregory quadrature weights.

    ``n_nodes`` and ``order`` must be integer scalars satisfying
    ``2 <= order <= n_nodes``. ``h`` is a finite real scalar; negative spacing
    reverses the integration direction, and zero spacing gives zero weights.
    For function values on the grid, the integral is ``values @ weights``.

    Order 2 is the trapezoidal rule. Higher orders correct both endpoints;
    overlapping corrections are added independently, as in the R original.
    Large orders can suffer cancellation and are not guaranteed to improve
    accuracy. Nonfinite intermediate or output values raise ValueError.
    """
    for name, value in (("n_nodes", n_nodes), ("order", order)):
        if isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{name} must be an integer scalar")
        try:
            integer = index(value)
        except TypeError as error:
            raise TypeError(f"{name} must be an integer scalar") from error
        if name == "n_nodes":
            n_nodes = integer
        else:
            order = integer
    if not 2 <= order <= n_nodes:
        raise ValueError("require 2 <= order <= n_nodes")
    if isinstance(h, (bool, np.bool_)) or not isinstance(h, Real):
        raise TypeError("h must be a real scalar")
    try:
        h = float(h)
    except OverflowError as error:
        raise ValueError("h must be finite") from error
    if not np.isfinite(h):
        raise ValueError("h must be finite")

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            r = 1.0 / np.arange(1, order + 1, dtype=float)
            size = order - 1
            rows, columns = np.indices((size, size))
            # Lower triangular Toeplitz matrix with first column r[:-1].
            top = np.zeros((size, size))
            lower = rows >= columns
            top[lower] = r[(rows - columns)[lower]]
            gc = np.linalg.solve(top, r[1:])
            pascal = np.zeros((size, size))
            for row in range(size):
                for column in range(row + 1):
                    pascal[row, column] = (-1) ** column * comb(row, column)
            updates = np.sum(h * gc[:, None] * pascal, axis=0)
            weights = np.full(n_nodes, h)
            weights[:size] -= updates
            weights[-size:] -= updates[::-1]
    except (FloatingPointError, OverflowError) as error:
        raise ValueError(
            "Gregory coefficients and weights must remain finite"
        ) from error
    if not np.all(np.isfinite(weights)):
        raise ValueError("Gregory coefficients and weights must remain finite")
    return weights
