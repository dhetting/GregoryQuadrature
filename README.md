# GregoryQuadrature
R package for performing numerical integration using Gregory's Method. To install, run `install.packages("GregoryQuadrature")` from your R console.

Gregory's Formula expresses an integral as a sum of equally spaced integrand values: $\int_{x_0}^{x_0 + nh} f(y)dy = h(\frac{1}{2}f_0 + f_1 + \cdots + f_{n-1} + \frac{1}{2}f_n) + \frac{h}{12}\left(\Delta f_0 - \Delta f_{n-1}\right) - \frac{h}{24}\left(\Delta^2f_0 + \Delta^2f_{n-2}\right) + \cdots$ where $f_j$ is $f(x_o + jh)$,  $\Delta^i$ are function differences, and $h$ is the step size. See Ralston 1965 for derivation details.

In contrast to Newton-Cotes approximation formulas, the values are equally weighted through the middle of the integral interval with corrections applied only at the end values. In this way, Gregory's Formula is an enhancement to the well-known Trapezoidal Rule, preserving spectral accuracy through the middle and decreasing error at the ends relative to the Trapezoidal Rule. Gregory's Method also avoids the evaluation of high-order derivatives required by Newton-Cotes formulas.

Gregory's Formula is exact for polynomials of degree $n$ if the summation includes differences through order $n$. The error can be reduced by adding more terms, but rounding error can become problematic with higher order differences (say, greater than 10). A composite formulation can help overcome this or see Fornberg and Reeger 2019 for orders up to 20.

To use `GregoryQuadtrature` evaluate the integrand at equal intervals, apply the weights provided by the `Gregory_weights` method, and sum.

For example, to evaluate the integral of $e^x$ over the interval $[-1, 1]$:

```{r}
# define number of abscissas
n_nodes = 11

# define order of evaluation
order = 8

# calculate step-size
h = 2/(n_nodes-1)

# evaluate the integrand
x = pracma::linspace(-1, 1, n_nodes)
f = exp(x)

# determine weights
w = Gregory_weights(n_nodes, h, order)

# calculate the integral
int = f %*% w

# calculate the error for known integrals
exact = exp(1) - exp(-1)
error = int - exact

print(error)
```


Ralston, Anthony. A First Course in Numerical Analysis. New York: McGraw-Hill, 1965.

Fornberg, B., Reeger, J.A. An improved Gregory-like method for 1-D quadrature. Numer. Math. 141, 1–19 (2019). https://doi-org.mines.idm.oclc.org/10.1007/s00211-018-0992-0
