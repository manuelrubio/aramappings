#' Title
#'
#' @param V
#' @param rank_V
#'
#' @returns
#' @export
#'
#' @examples
compute_orthogonal_projection_matrix <- function(
    V,
    rank_V) {
  m <- ncol(V)

  # Compute matrix U whose columns form a basis for the subspace spanned by the rows of V (of dimension rank_V)
  U <- zeros(m, rank_V)

  k <- 1
  while (nnz(V[k, ]) == 0) {
    k <- k + 1
  }

  U[, 1] <- t(V[k, ])

  if (rank_V > 1) {
    j <- k + 1
    U[, 2] <- t(V[j, ])

    while (rankMatrix(U)[1] < rank_V) {
      j <- j + 1
      U[, 2] <- t(V[j, ])
    }
  }

  # Orthogonal projection matrix
  U %*% pinv(U)
}


#' Title
#'
#' @param n_points
#' @param V
#' @param row_offset
#'
#' @returns
#' @export
#'
#' @examples
ara_L1_norm_coo_lists <- function(
    n_points,
    V,
    row_offset) {
  n <- nrow(V)
  m <- ncol(V)

  aux_vec <- repmat(pracma::Reshape((1 + row_offset):(2 * n * n_points + row_offset), 2 * n, n_points), m, 1)
  rows <- c(
    ((1 + row_offset):(2 * n * n_points + row_offset)),
    unlist(as.list(aux_vec))
  )

  aux_vec <- repmat(pracma::Reshape(1:(n * n_points), n, n_points), 2, 1)
  cols <- c(
    unlist(as.list(aux_vec)),
    rep((n_points * n + 1):(n_points * (m + n)), each = 2 * n)
  )

  vals <- c(
    rep(-1, 2 * n * n_points),
    repmat(unlist(as.list(rbind(-V, V))), 1, n_points)
  )

  list(rows = rows, cols = cols, vals = vals)
}


#' Title
#'
#' @param n_points
#' @param V
#' @param row_offset
#'
#' @returns
#' @export
#'
#' @examples
ara_Linf_norm_coo_lists <- function(
    n_points,
    V,
    row_offset) {
  n <- nrow(V)
  m <- ncol(V)

  aux_vec <- repmat(pracma::Reshape((1 + row_offset):(2 * n * n_points + row_offset), 2 * n, n_points), m, 1)
  rows <- c(
    (1 + row_offset):(2 * n * n_points + row_offset),
    unlist(as.list(aux_vec))
  )

  cols <- rep(1:(n_points * (m + 1)), each = 2 * n)

  vals <- c(
    rep(-1, 2 * n * n_points),
    repmat(unlist(as.list(rbind(-V, V))), 1, n_points)
  )

  list(rows = rows, cols = cols, vals = vals)
}



#' Title
#'
#' @param axis_vector
#' @param sort_indices
#' @param row_offset
#' @param col_offset
#'
#' @returns
#' @export
#'
#' @examples
ara_ordered_inequality_coo_lists <- function(
    axis_vector,
    sort_indices,
    row_offset,
    col_offset) {
  N <- length(sort_indices)
  m <- length(axis_vector)

  ranks <- order(sort_indices) # Data ranks for selected variable (in increasing order)

  rows <- rep(NA, 2 * m * (N - 1))
  cols <- rep(NA, 2 * m * (N - 1))
  vals <- rep(NA, 2 * m * (N - 1))

  v_block <- unlist(as.list(t(cbind(-axis_vector, axis_vector))))

  idx <- 1
  col <- col_offset + 1
  for (i in 1:N) {
    if (ranks[i] == 1) {
      rows[idx:(idx + m - 1)] <- rep(ranks[i], m) + row_offset
      cols[idx:(idx + m - 1)] <- col:(col + m - 1)
      vals[idx:(idx + m - 1)] <- axis_vector
      idx <- idx + m
      col <- col + m
    } else if (ranks[i] == N) {
      rows[idx:(idx + m - 1)] <- rep(N - 1, m) + row_offset
      cols[idx:(idx + m - 1)] <- col:(col + m - 1)
      vals[idx:(idx + m - 1)] <- -axis_vector
      idx <- idx + m
      col <- col + m
    } else {
      rows[idx:(idx + 2 * m - 1)] <- rep((ranks[i] - 1 + row_offset):(ranks[i] + row_offset), m)
      cols[idx:(idx + 2 * m - 1)] <- rep(col:(col + m - 1), each = 2)
      vals[idx:(idx + 2 * m - 1)] <- v_block
      idx <- idx + 2 * m
      col <- col + m
    }
  }

  list(
    rows = rows,
    cols = cols,
    vals = vals
  )
}


#' Title
#'
#' @param nrows
#' @param ncols
#' @param kind
#' @param clower
#' @param cupper
#' @param obj
#' @param type_cols
#' @param rlower
#' @param rupper
#' @param type_rows
#' @param ne
#' @param rows
#' @param cols
#' @param vals
#' @param use_interior_point
#' @param init_index
#' @param n_points
#' @param m
#'
#' @returns
#' @export
#'
#' @examples
solve_glpkAPI_wrapper <- function(
    nrows,
    ncols,
    kind,
    clower,
    cupper,
    obj,
    type_cols,
    rlower,
    rupper,
    type_rows,
    ne,
    rows,
    cols,
    vals,
    use_interior_point,
    init_index,
    n_points,
    m) {
  lp <- initProbGLPK()

  setObjDirGLPK(lp, GLP_MIN)

  addRowsGLPK(lp, nrows)
  addColsGLPK(lp, ncols)

  setColsKindGLPK(lp, j = c(1:ncols), kind)

  setColsBndsObjCoefsGLPK(lp = lp, j = c(1:ncols), lb = clower, ub = cupper, obj_coef = obj, type = type_cols)

  setRowsBndsGLPK(lp, i = c(1:nrows), lb = rlower, ub = rupper, type = type_rows)


  # load constraint matrix
  loadMatrixGLPK(lp = lp, ne = ne, ia = rows, ja = cols, ra = vals)


  # solve linear problem
  if (use_interior_point) {
    setInteriorParmGLPK(MSG_LEV, GLP_MSG_OFF)
    solveInteriorGLPK(lp)
    solution_status <- getSolStatIptGLPK(lp)

    if ((solution_status == 3) || (solution_status == 5)) {
      x <- getColsPrimIptGLPK(lp)

      P <- t(pracma::Reshape(x[init_index:(init_index + n_points * m - 1)], m, n_points))
      objval <- getObjValIptGLPK(lp)
    } else {
      print("Error: glpkAPI failed to compute an optimal solution")

      P <- matrix(NA, n_points, m)
      objval <- NA
    }
  } else {
    setSimplexParmGLPK(MSG_LEV, GLP_MSG_OFF)
    solveSimplexGLPK(lp)
    solution_status <- getSolStatGLPK(lp)

    if (solution_status == 5) {
      x <- getColsPrimGLPK(lp)

      P <- t(pracma::Reshape(x[init_index:(init_index + n_points * m - 1)], m, n_points))
      objval <- getObjValGLPK(lp)
    } else {
      print("Error: glpkAPI failed to compute an optimal solution")

      P <- matrix(NA, n_points, m)
      objval <- NA
    }
  }

  status <- rep(solution_status, n_points)

  # remove problem object
  delProbGLPK(lp)

  list(
    P = P,
    status = status,
    objval = objval
  )
}


#' Title
#'
#' @param A
#' @param b
#' @param obj
#' @param cones
#' @param P
#' @param init_index
#' @param n_points
#' @param m
#'
#' @returns
#' @export
#'
#' @examples
solve_clarabel_wrapper <- function(
    A,
    b,
    obj,
    cones,
    P = NULL,
    init_index,
    n_points,
    m) {
  clarabel_output <- clarabel(
    A = A,
    b = b,
    q = obj,
    P = P,
    cones = cones,
    control = list(verbose = FALSE),
    strict_cone_order = TRUE
  )

  if (clarabel_output$status == 2) {
    x <- clarabel_output$x

    P <- t(pracma::Reshape(x[init_index:(init_index + n_points * m - 1)], m, n_points))
    objval <- clarabel_output$obj_val
  } else {
    print("Error: clarabel failed to compute an optimal solution")

    P <- matrix(NA, n_points, m)
    objval <- NA
  }

  status <- rep(clarabel_output$status, n_points)

  list(
    P = P,
    status = status,
    objval = objval
  )
}


#' Title
#'
#' @param A
#' @param b
#' @param obj
#' @param bounds
#' @param dirs
#' @param init_index
#' @param n_points
#' @param m
#'
#' @returns
#' @export
#'
#' @examples
solve_Rglpk_wrapper <- function(
    A,
    b,
    obj,
    bounds,
    dirs,
    init_index,
    n_points,
    m) {
  rglpk_output <- Rglpk_solve_LP(
    obj = obj,
    mat = A,
    dir = dirs,
    rhs = b,
    bounds = bounds,
    types = NULL,
    max = FALSE,
    canonicalize_status = FALSE
  )

  if (rglpk_output$status == 5) {
    x <- rglpk_output$solution

    P <- t(pracma::Reshape(x[init_index:(init_index + n_points * m - 1)], m, n_points))
    objval <- rglpk_output$optimum
  } else {
    print("Error: Rglpk failed to compute an optimal solution")

    P <- matrix(NA, n_points, m)
    objval <- NA
  }

  status <- rep(rglpk_output$status, n_points)

  list(
    P = P,
    status = status,
    objval = objval
  )
}




#' Title
#'
#' @param x
#' @param nrows
#' @param ncols
#' @param kind
#' @param clower
#' @param cupper
#' @param obj
#' @param type_cols
#' @param rlower
#' @param type_rows
#' @param ne
#' @param rows
#' @param cols
#' @param vals
#' @param use_interior_point
#' @param m
#'
#' @returns
#' @export
#'
#' @examples
min_unconstrained_glpkAPI_point <- function(
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
    rows,
    cols,
    vals,
    use_interior_point,
    m) {
  rupper <- c(cbind(-x, x))

  glpkAPI_output <- solve_glpkAPI_wrapper(
    nrows,
    ncols,
    kind,
    clower,
    cupper,
    obj,
    type_cols,
    rlower,
    rupper,
    type_rows,
    ne,
    rows,
    cols,
    vals,
    use_interior_point,
    length(obj) - m + 1,
    1,
    m
  )
}


#' Title
#'
#' @param x
#' @param A
#' @param obj
#' @param cones
#' @param m
#'
#' @returns
#' @export
#'
#' @examples
min_unconstrained_clarabel_point <- function(
    x,
    A,
    obj,
    cones,
    m) {
  b <- c(cbind(-x, x))

  clarabel_output <- solve_clarabel_wrapper(
    A,
    b,
    obj,
    cones,
    NULL,
    length(obj) - m + 1,
    1,
    m
  )
}


#' Title
#'
#' @param x
#' @param A
#' @param obj
#' @param bounds
#' @param dirs
#' @param m
#'
#' @returns
#' @export
#'
#' @examples
min_unconstrained_Rglpk_point <- function(
    x,
    A,
    obj,
    bounds,
    dirs,
    m) {
  b <- c(cbind(-x, x))

  rglpk_output <- solve_Rglpk_wrapper(
    A,
    b,
    obj,
    bounds,
    dirs,
    length(obj) - m + 1,
    1,
    m
  )
}





#' Title
#'
#' @param x
#' @param nrows
#' @param ncols
#' @param kind
#' @param clower
#' @param cupper
#' @param obj
#' @param type_cols
#' @param type_rows
#' @param ne
#' @param rows
#' @param cols
#' @param vals
#' @param variable
#' @param use_interior_point
#' @param m
#'
#' @returns
#' @export
#'
#' @examples
min_exact_glpkAPI_point <- function(
    x,
    nrows,
    ncols,
    kind,
    clower,
    cupper,
    obj,
    type_cols,
    type_rows,
    ne,
    rows,
    cols,
    vals,
    variable,
    use_interior_point,
    m) {
  n <- length(x)

  # Equality constraint for x[variable]
  rlower <- c(x[variable], rep(-Inf, 2 * n))
  rupper <- c(x[variable], cbind(-x, x))

  solve_glpkAPI_wrapper(
    nrows,
    ncols,
    kind,
    clower,
    cupper,
    obj,
    type_cols,
    rlower,
    rupper,
    type_rows,
    ne,
    rows,
    cols,
    vals,
    use_interior_point,
    length(obj) - m + 1,
    1,
    m
  )
}



#' Title
#'
#' @param x
#' @param A
#' @param obj
#' @param cones
#' @param variable
#' @param m
#'
#' @returns
#' @export
#'
#' @examples
min_exact_clarabel_point <- function(
    x,
    A,
    obj,
    cones,
    variable,
    m) {
  b <- c(x[variable], cbind(-x, x))

  clarabel_output <- solve_clarabel_wrapper(
    A,
    b,
    obj,
    cones,
    NULL,
    length(obj) - m + 1,
    1,
    m
  )
}


#' Title
#'
#' @param x
#' @param A
#' @param obj
#' @param bounds
#' @param dirs
#' @param variable
#' @param m
#'
#' @returns
#' @export
#'
#' @examples
min_exact_Rglpk_point <- function(
    x,
    A,
    obj,
    bounds,
    dirs,
    variable,
    m) {
  b <- c(x[variable], cbind(-x, x))

  rglpk_output <- solve_Rglpk_wrapper(
    A,
    b,
    obj,
    bounds,
    dirs,
    length(obj) - m + 1,
    1,
    m
  )
}
