% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ara_unconstrained_Linf.R
\name{ara_unconstrained_Linf}
\alias{ara_unconstrained_Linf}
\title{Unconstrained Adaptable Radial Axes (ARA) mappings using the L-infinity norm}
\usage{
ara_unconstrained_Linf(
  X,
  V,
  weights = rep(1, ncol(X)),
  solver = "glpkAPI",
  use_glpkAPI_simplex = TRUE,
  cluster = NULL
)
}
\arguments{
\item{X}{Numeric data matrix of dimensions N x n, where N is the number of observations, and n is the number of variables.}

\item{V}{Numeric matrix defining the axes or "axis vectors". Its dimensions are n x m, where 1<=m<=3 is the dimension
of the visualization space. Each row of \code{V} defines an axis vector.}

\item{weights}{Numeric array specifying optional non-negative weights associated with each variable. The function only considers
them if they do not share the same value. Default: array of n ones.}

\item{solver}{String indicating a package for solving the linear problem(s). It can be "clarabel" (default), "glpkAPI", "Rglpk", or "CVXR".}

\item{use_glpkAPI_simplex}{Boolean parameter that indicates whether to use the simplex algorithm (if \code{TRUE}) or an interior point method
(if \code{FALSE}), when using the \pkg{glpkAPI} solver. The default is \code{TRUE}.}

\item{cluster}{Optional cluster object related to the parallel package. If supplied, and \code{n_LP_problems} is N, the method
computes the mappings using parallel processing.}
}
\value{
A list with the three following entries:
\describe{
\item{\code{P}}{A numeric N x m matrix containing the mapped points. Each row is the low-dimensional representation
of a data observation in X.}
\item{\code{status}}{A vector of length N where the i-th element contains the status of the chosen solver when
calculating the mapping of the i-th data observation. The type of the elements depends on the particular chosen solver.}
\item{\code{objval}}{The numeric objective value associated with the solution to the optimization problem, considering
matrix norms.}
}
If the chosen solver fails to map one or more data observations (i.e., fails to solve the related optimization problems),
their rows in \code{P} will contain \code{NA} (not available) values. In that case, \code{objval} will also be \code{NA}.
}
\description{
\code{ara_unconstrained_Linf()} computes \strong{unconstrained} \strong{Adaptable Radial Axes} (ARA) mappings for
the \strong{L-infinity norm}
}
\details{
\code{ara_unconstrained_Linf()} computes low-dimensional point representations of high-dimensional
numerical data (\code{X}) according to the data visualization method "Adaptable Radial Axes" (Rubio-Sánchez, 2017),
which describes a collection of convex norm optimization problems aimed at minimizing estimates of original
values in \code{X} through dot products of the mapped points with the axis vectors (rows of \code{V}). This particular
function solves the unconstrained optimization problem in Eq. (10), for the L-infinity vector norm. Specifically, it solves
equivalent linear problems as described in (12). Optional non-negative weights (\code{weights}) associated with each data
variable can be supplied to solve the problem in Eq. (15).
}
\examples{
# Load data
library(ascentTraining) # contains the Auto MPG dataset
data("auto_mpg")

# Define subset of (numerical) variables
selected_variables <- c(1,4,5,6)   # 1:"mpg", 4:"horsepower", 5:"weight", 6:"acceleration")

# Retain only selected variables and rename dataset as X
X <- auto_mpg[, selected_variables] # Select a subset of variables
rm(auto_mpg)

# Remove rows with missing values from X
N <- nrow(X)
rows_to_delete <- NULL
for (i in 1:N) {
  if (sum(is.na(X[i, ])) > 0) {
    rows_to_delete <- c(rows_to_delete, -i)
  }
}
X <- X[rows_to_delete, ]

# Convert X to matrix
X <- apply(as.matrix.noquote(X), 2, as.numeric)

# Standardize data
Z <- scale(X)

# Define axis vectors (2-dimensional in this example)
library(geometry)
r <- c(0.8, 1, 1.2, 1)
theta <- c(225, 100, 315, 80) * 2 * pi / 360
V <- pol2cart(theta, r)

# Define weights
weights <- c(1, 0.75, 0.75, 1)

# Detect the number of available CPU cores
library(parallelly)
NCORES <- parallelly::availableCores(omit = 1)

# Create a cluster for parallel processing
cl <- parallel::makeCluster(NCORES)

# Compute the mapping
mapping <- ara_unconstrained_Linf(
  Z,
  V,
  weights = weights,
  solver = "glpkAPI",
  use_glpkAPI_simplex = TRUE,
  cluster = cl
)

# Stop cluster
parallel::stopCluster(cl)

# Select variables with labeled axis lines on ARA plot
axis_lines <- c(1,4)   # 1:"mpg", 4:"acceleration")

# Select variable used for coloring embedded points
color_variable <- 1    # "mpg"

# Draw the ARA plot
draw_ara_plot_2d_standardized(
  Z,
  X,
  V,
  mapping$P,
  weights = weights,
  axis_lines = axis_lines,
  color_variable = color_variable
)

}
\references{
M. Rubio-Sánchez, A. Sanchez, D. J. Lehmann: Adaptable radial axes plots for improved
multivariate data visualization. Computer Graphics Forum 36, 3 (2017), 389–399.
\href{https://onlinelibrary.wiley.com/doi/10.1111/cgf.13196}{doi:10.1111/cgf.13196}
}
