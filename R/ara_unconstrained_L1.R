#' Unconstrained Adaptable Radial Axes (ARA) mappings using the L1 norm
#'
#' @description
#' \code{ara_unconstrained_L1()} computes \strong{unconstrained} \strong{Adaptable Radial Axes} (ARA) mappings for the \strong{L1 norm}
#'
#' @details
#' \code{ara_unconstrained_L1()} computes low-dimensional point representations of high-dimensional
#' numerical data (\code{X}) according to the data visualization method "Adaptable Radial Axes" (Rubio-Sánchez, 2017),
#' which describes a collection of convex norm optimization problems aimed at minimizing estimates of original
#' values in \code{X} through dot products of the mapped points with the axis vectors (rows of \code{V}). This particular
#' function solves the unconstrained optimization problem in Eq. (10), for the L1 vector norm. Specifically, it solves
#' equivalent linear problems as described in (11). Optional non-negative weights (\code{weights}) associated with each data
#' variable can be supplied to solve the problem in Eq. (15).
#'
#' @param X
#' Numeric data matrix of dimensions N x n, where N is the number of observations, and n is the number of variables.
#' @param V
#' Numeric matrix defining the axes or "axis vectors". Its dimensions are n x m, where 1<=m<=3 is the dimension
#' of the visualization space. Each row of \code{V} defines an axis vector.
#' @param weights
#' Numeric array specifying optional non-negative weights associated with each variable. The function only considers
#' them if they do not share the same value. Default: array of n ones.
#' @param solver
#' String indicating a package for solving the linear problem(s). It can be "clarabel" (default), "glpkAPI", "Rglpk", or "CVXR".
#' @param use_glpkAPI_simplex
#' Boolean parameter that indicates whether to use the simplex algorithm (if \code{TRUE}) or an interior point method
#' (if \code{FALSE}), when using the \pkg{glpkAPI} solver. The default is \code{TRUE}.
#' @param cluster
#' Optional cluster object related to the parallel package. If supplied, and \code{n_LP_problems} is N, the method
#' computes the mappings using parallel processing.
#'
#' @returns
#' A list with the three following entries:
#' \describe{
#'   \item{\code{P}}{A numeric N x m matrix containing the mapped points. Each row is the low-dimensional representation
#'   of a data observation in X.}
#'   \item{\code{status}}{A vector of length N where the i-th element contains the status of the chosen solver when
#'   calculating the mapping of the i-th data observation. The type of the elements depends on the particular chosen solver.}
#'   \item{\code{objval}}{The numeric objective value associated with the solution to the optimization problem, considering
#'   matrix norms.}
#' }
#' If the chosen solver fails to map one or more data observations (i.e., fails to solve the related optimization problems),
#' their rows in \code{P} will contain \code{NA} (not available) values. In that case, \code{objval} will also be \code{NA}.
#'
#' @references
#' M. Rubio-Sánchez, A. Sanchez, D. J. Lehmann: Adaptable radial axes plots for improved
#' multivariate data visualization. Computer Graphics Forum 36, 3 (2017), 389–399.
#' [doi:10.1111/cgf.13196](https://onlinelibrary.wiley.com/doi/10.1111/cgf.13196)
#'
#' @export
#'
#' @examples
#' # Load data
#' library(ascentTraining) # contains the Auto MPG dataset
#' data("auto_mpg")
#'
#' # Define subset of (numerical) variables
#' selected_variables <- c(1,4,5,6)   # 1:"mpg", 4:"horsepower", 5:"weight", 6:"acceleration")
#'
#' # Retain only selected variables and rename dataset as X
#' X <- auto_mpg[, selected_variables] # Select a subset of variables
#' rm(auto_mpg)
#'
#' # Remove rows with missing values from X
#' N <- nrow(X)
#' rows_to_delete <- NULL
#' for (i in 1:N) {
#'   if (sum(is.na(X[i, ])) > 0) {
#'     rows_to_delete <- c(rows_to_delete, -i)
#'   }
#' }
#' X <- X[rows_to_delete, ]
#'
#' # Convert X to matrix
#' X <- apply(as.matrix.noquote(X), 2, as.numeric)
#'
#' # Standardize data
#' Z <- scale(X)
#'
#' # Define axis vectors (2-dimensional in this example)
#' library(geometry)
#' r <- c(0.8, 1, 1.2, 1)
#' theta <- c(225, 100, 315, 80) * 2 * pi / 360
#' V <- pol2cart(theta, r)
#'
#' # Define weights
#' weights <- c(1, 0.75, 0.75, 1)
#'
#' # Detect the number of available CPU cores
#' library(parallelly)
#' NCORES <- parallelly::availableCores(omit = 1)
#'
#' # Create a cluster for parallel processing
#' cl <- parallel::makeCluster(NCORES)
#'
#' # Compute the mapping
#' mapping <- ara_unconstrained_L1(
#'   Z,
#'   V,
#'   weights = weights,
#'   solver = "glpkAPI",
#'   use_glpkAPI_simplex = TRUE,
#'   cluster = cl
#' )
#'
#' # Stop cluster
#' parallel::stopCluster(cl)
#'
#' # Select variables with labeled axis lines on ARA plot
#' axis_lines <- c(1,4)   # 1:"mpg", 4:"acceleration")
#'
#' # Select variable used for coloring embedded points
#' color_variable <- 1    # "mpg"
#'
#' # Draw the ARA plot
#' draw_ara_plot_2d_standardized(
#'   Z,
#'   X,
#'   V,
#'   mapping$P,
#'   weights = weights,
#'   axis_lines = axis_lines,
#'   color_variable = color_variable
#' )
#'
ara_unconstrained_L1 <- function(
    X,
    V,
    weights = rep(1, ncol(X)),
    solver = "glpkAPI",
    use_glpkAPI_simplex = TRUE,
    cluster = NULL) {
  ###################   Check validity of input parameters   ###################

  # Check data types -----------------------------------------------------------

  if ((!inherits(X, "matrix") && !inherits(X, "array")) || !is.double(X)) {
    stop("Input error: X must be a numeric matrix")
  }

  if ((!inherits(V, "matrix") && !inherits(V, "array")) || !is.double(V)) {
    stop("Input error: V must be a numeric matrix")
  }

  if (!is.double(weights)) {
    stop("Input error: weights must be a numeric vector")
  }

  if (!is.character(solver)) {
    stop("Input error: solver must be a string")
  }

  if (!is.logical(use_glpkAPI_simplex)) {
    stop("Input error: use_glpkAPI_simplex must be logical (Boolean)")
  }

  if ((!is.null(cluster)) && !(inherits(cluster, "SOCKcluster") || inherits(cluster, "cluster"))) {
    stop("Input error: invalid cluster argument")
  }


  # Check dimensions of matrices -----------------------------------------------

  N <- nrow(X)
  nX <- ncol(X)
  nV <- nrow(V)
  m <- ncol(V)

  if ((m < 1) || (m > 3)) {
    stop("Input error: The dimensionality of the visualization space (columns of V) must be 1, 2, or 3")
  }

  if (nX != nV) {
    stop("Input error: The number of variables of X (columns) must match the number of variables of V (rows)")
  }

  n <- nX


  # Check that inputs have no missing values -----------------------------------

  if (any(is.na(X))) {
    stop("Input error: X cannot contain missing values)")
  }

  if (any(is.na(V))) {
    stop("Input error: V cannot contain missing values)")
  }

  if (any(is.na(weights))) {
    stop("Input error: weights cannot contain missing values)")
  }


  # Check additional preconditions on input parameters -------------------------

  if ((length(weights) != n) || (min(weights) < 0)) {
    stop("Input error: weights must be vector of length n with non-negative values")
  }

  if ((!pracma::strcmpi(solver, "clarabel")) && (!pracma::strcmpi(solver, "glpkAPI")) && (!pracma::strcmpi(solver, "Rglpk")) && (!pracma::strcmpi(solver, "CVXR"))) {
    stop('Input error: solver must be "clarabel", "glpkAPI", "Rglpk", or "CVXR"')
  }


  ############################   Compute mapping   #############################

  if (length(unique(weights)) > 1) {
    W <- diag(weights)
    X <- X %*% W
    V <- W %*% V
  }

  if (pracma::strcmpi(solver, "CVXR")) {
    outputs <- ara_unconstrained_L1_CVXR(
      X,
      V
    )
  } else {
    if (pracma::strcmpi(solver, "glpkAPI")) {
      outputs <- ara_unconstrained_L1_glpkAPI(
        X,
        V,
        use_glpkAPI_simplex,
        cluster
      )
    } else if (pracma::strcmpi(solver, "clarabel")) {
      outputs <- ara_unconstrained_L1_clarabel(
        X,
        V,
        cluster
      )
    } else {
      outputs <- ara_unconstrained_L1_Rglpk(
        X,
        V,
        cluster
      )
    }
    if (m == 1) {
      outputs$P <- t(outputs$P)
    }
  }

  if (!is.na(outputs$objval)) {
    rank_V <- Matrix::rankMatrix(V)[1]
    if (rank_V < m) {
      Q <- compute_orthogonal_projection_matrix(V, rank_V)
      outputs$P <- outputs$P %*% Q
    }
  }

  list(
    P = outputs$P,
    status = outputs$status,
    objval = outputs$objval
  )
}



#' @noRd
ara_unconstrained_L1_CVXR <- function(
    X,
    V) {
  N <- nrow(X)
  n <- ncol(X)
  m <- ncol(V)

  Pvar <- CVXR::Variable(N, m)
  Tvar <- CVXR::Variable(N, n)

  obj <- CVXR::Minimize(sum(Tvar))

  constraints <- list()
  constraints <- append(constraints, -Tvar <= Pvar %*% t(V) - X)
  constraints <- append(constraints, Pvar %*% t(V) - X <= Tvar)

  prob <- CVXR::Problem(obj, constraints)
  solution <- CVXR::solve(prob, solver = "ECOS")

  extract_CVXR_points_status_objval(
    solution,
    Pvar,
    V,
    N,
    m
  )
}



#' @noRd
ara_unconstrained_L1_glpkAPI <- function(
    X,
    V,
    use_glpkAPI_simplex,
    cluster) {
  N <- nrow(X)
  n <- ncol(X)
  m <- ncol(V)

  obj <- c(rep(1, n), rep(0, m))

  coo_lists <- ara_L1_norm_coo_lists(1, V, 0)

  ne <- 2 * n * (m + 1)

  nrows <- 2 * n
  ncols <- n + m

  kind <- rep(glpkAPI::GLP_CV, ncols)
  type_cols <- c(rep(glpkAPI::GLP_LO, n), rep(glpkAPI::GLP_FR, m))
  clower <- c(rep(0, n), rep(-Inf, m))
  cupper <- rep(Inf, n + m)

  type_rows <- rep(glpkAPI::GLP_UP, nrows)
  rlower <- rep(0, 2 * n)

  if (is.null(cluster)) {
    sol <- apply(X = X, MARGIN = 1, function(x) {
      min_unconstrained_glpkAPI(
        x,
        nrows,
        ncols,
        kind,
        clower,
        cupper,
        obj,
        type_cols,
        rlower,
        type_rows,
        ne,
        coo_lists$rows,
        coo_lists$cols,
        coo_lists$vals,
        use_glpkAPI_simplex,
        m
      )
    })
  } else {
    parallel::clusterEvalQ(cluster, library(glpkAPI))
    parallel::clusterExport(cluster, c("min_unconstrained_glpkAPI", "solve_glpkAPI_wrapper"), envir = environment())

    sol <- parallel::parApply(cluster, X = X, MARGIN = 1, function(x) {
      min_unconstrained_glpkAPI(
        x,
        nrows,
        ncols,
        kind,
        clower,
        cupper,
        obj,
        type_cols,
        rlower,
        type_rows,
        ne,
        coo_lists$rows,
        coo_lists$cols,
        coo_lists$vals,
        use_glpkAPI_simplex,
        m
      )
    })
  }

  sol_matrix <- pracma::Reshape(unlist(sol), m + 2, N)

  list(
    P = t(sol_matrix[1:m, ]),
    status = sol_matrix[m + 1, ],
    objval = sum(sol_matrix[m + 2, ])
  )
}


#' @noRd
ara_unconstrained_L1_clarabel <- function(
    X,
    V,
    cluster) {
  N <- nrow(X)
  n <- ncol(X)
  m <- ncol(V)

  obj <- c(rep(1, n), rep(0, m))
  # obj <- c(matrix(1, 1, n), matrix(0, 1, m))

  A <- rbind(cbind(diag(-1, n), -V), cbind(diag(-1, n), V)) # dense matrix

  cones <- list(l = (2 * n))

  if (is.null(cluster)) {
    sol <- apply(X = X, MARGIN = 1, function(x) {
      min_unconstrained_clarabel(
        x,
        A,
        obj,
        cones,
        m
      )
    })
  } else {
    parallel::clusterEvalQ(cluster, library(clarabel))
    parallel::clusterExport(cluster, c("min_unconstrained_clarabel", "solve_clarabel_wrapper"), envir = environment())

    sol <- parallel::parApply(cluster, X = X, MARGIN = 1, function(x) {
      min_unconstrained_clarabel(
        x,
        A,
        obj,
        cones,
        m
      )
    })
  }

  sol_matrix <- pracma::Reshape(unlist(sol), m + 2, N)

  list(
    P = t(sol_matrix[1:m, ]),
    status = sol_matrix[m + 1, ],
    objval = sum(sol_matrix[m + 2, ])
  )
}


#' @noRd
ara_unconstrained_L1_Rglpk <- function(
    X,
    V,
    cluster) {
  N <- nrow(X)
  n <- ncol(X)
  m <- ncol(V)

  obj <- c(rep(1, n), rep(0, m))

  coo_lists <- ara_L1_norm_coo_lists(1, V, 0)

  A <- slam::simple_triplet_matrix(
    coo_lists$rows,
    coo_lists$cols,
    coo_lists$vals,
    nrow = 2 * n,
    ncol = n + m
  )

  bounds <- list(lower = list(ind = (n + 1):(n + m), val = rep(-Inf, m)), upper = list())

  dirs <- rep("<=", 2 * n)

  if (is.null(cluster)) {
    sol <- apply(X = X, MARGIN = 1, function(x) {
      min_unconstrained_Rglpk(
        x,
        A,
        obj,
        bounds,
        dirs,
        m
      )
    })
  } else {
    parallel::clusterEvalQ(cluster, library(Rglpk))
    parallel::clusterExport(cluster, c("min_unconstrained_Rglpk", "solve_Rglpk_wrapper"), envir = environment())

    sol <- parallel::parApply(cluster, X = X, MARGIN = 1, function(x) {
      min_unconstrained_Rglpk(
        x,
        A,
        obj,
        bounds,
        dirs,
        m
      )
    })
  }

  sol_matrix <- pracma::Reshape(unlist(sol), m + 2, N)

  list(
    P = t(sol_matrix[1:m, ]),
    status = sol_matrix[m + 1, ],
    objval = sum(sol_matrix[m + 2, ])
  )
}
