#' Calculate the Gregory quadrature weights for equispaced integration. If f is
#' a row vector containing the function values, the integral is approximated by 
#' the statement `f %*% t(w)` where w are the returned weights. Translated
#'  from https://www.colorado.edu/amath/sites/default/files/attached-files/gregory.pdf.
#' 
#' @param n_nodes Total number of nodes
#' @param h Step size
#' @param order Order of accuracy desired. 2, 3, 4, ... (with 2 giving the trapezoidal rule). The value must satisfy 2 <= order <= n_nodes
#'
#' @export Gregory_weights
#' @returns The weights to be used for the successive function values
#' @examples
#' n_nodes = 11
#' order = 8
#' h = 2/(n_nodes-1)
#' x = pracma::linspace(-1, 1, n_nodes)
#' f = exp(x)
#'
#' w = GregoryQuadrature::Gregory_weights(n_nodes, h, order)
#' int = f %*% w

#' # Exact value for integral
#' exact = exp(1) - exp(-1)
#'
#' error = int - exact
Gregory_weights <- function(n_nodes, h, order) {

  # Create the sequence of Gregory coefficients
  r <- 1. / seq(1, order)
  top_c <- r[1:order - 1]
  top_r <- c(r[1], replicate(order - 2, 0))
  top <- pracma::Toeplitz(top_c, top_r)
  gc <- solve(top, r[2:order])

  # Create the weights vector w and then update it at the two ends of the interval
  w <- replicate(n_nodes, 1) * h

  gc_rep <- replicate(order - 1, gc)

  p_chol <- pracma::pascal(order - 1, 1)
  w_updates <- colSums(h * gc_rep * p_chol)

  w[1:order - 1] <- w[1:order - 1] - w_updates[1:order - 1]

  w[(n_nodes - order + 2):n_nodes] <- w[(n_nodes - order + 2):n_nodes] - rev(w_updates[1:order - 1])

  return(w)
}

