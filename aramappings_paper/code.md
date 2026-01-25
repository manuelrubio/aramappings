

``` r
############################### Clear environment ##############################

message("***********************************************************************************************\n")
```

```
## ***********************************************************************************************
```

``` r
message("Replication script for 'aramappings: An R Package for Computing Adaptable Radial Axes Mappings'\n")
```

```
## Replication script for 'aramappings: An R Package for Computing Adaptable Radial Axes Mappings'
```

``` r
message("***********************************************************************************************\n\n")
```

```
## ***********************************************************************************************
```

``` r
############################### Set working directory ##########################

if (Sys.getenv("RSTUDIO") == "1") {

  rm(list = ls())

  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory(function() {})[1])
}
```

```
## Error in setwd(dirname(rstudioapi::getActiveDocumentContext()$path)): no es posible cambiar el directorio de trabajo
```

``` r
message("Generating paper's radial axes plots in folder: ", getwd(), "\n")
```

```
## Generating paper's radial axes plots in folder: O:/GIT_working_directory_NPC/JSS/Version01/article/replication_script
```

``` r
############################### Load packages ##################################

if (!require(aramappings)) {
  message("Trying to install aramappings")
  install.packages("aramappings")
  if(!require(aramappings)){
    stop("could not install aramappings")
  }
}

if (!require(ggplot2)) {
  message("Trying to install package ggplot2")
  install.packages("ggplot2")
  if(!require(ggplot2)){
    stop("Could not install package ggplot2")
  }
}

if (!require(pracma)) {
  message("Trying to install package pracma")
  install.packages("pracma")
  if(!require(pracma)){
    stop("Could not install package pracma")
  }
}

if (!require(geometry)) {
  message("Trying to install package geometry")
  install.packages("geometry")
  if(!require(geometry)){
    stop("Could not install package geometry")
  }
}

if (!require(stats)) {
  message("Trying to install package stats")
  install.packages("stats")
  if(!require(stats)){
    stop("Could not install package stats")
  }
}

if (!require(grid)) {
  message("Trying to install package grid")
  install.packages("grid")
  if(!require(grid)){
    stop("Could not install package grid")
  }
}

if (!require(RSpectra)) {
  message("Trying to install package RSpectra")
  install.packages("RSpectra")
  if(!require(RSpectra)){
    stop("Could not install package RSpectra")
  }
}

if (!require(grDevices)) {
  message("Trying to install package grDevices")
  install.packages("grDevices")
  if(!require(grDevices)){
    stop("Could not install package grDevices")
  }
}

if (!require(viridis)) {
  message("Trying to install package viridis")
  install.packages("viridis")
  if(!require(viridis)){
    stop("Could not install package viridis")
  }
}

if (!require(liver)) {
  message("Trying to install package liver")
  install.packages("liver")
  if(!require(liver)){
    stop("Could not install package liver")
  }
}

if (!require(ascentTraining)) {
  message("Trying to install package ascentTraining")
  install.packages("ascentTraining")
  if(!require(ascentTraining)){
    stop("Could not install package ascentTraining")
  }
}


############################### Auxiliary functions ############################

radviz <- function(X, V, normalize = TRUE) {
  EPSILON <- 10^(-12)

  N <- nrow(X)
  n <- ncol(X)

  if (normalize) {
    # min-max normalization
    minimum <- apply(X, 2, min)
    maximum <- apply(X, 2, max)
    X <- (X - pracma::repmat(minimum, N, 1)) / (pracma::repmat(maximum, N, 1) - pracma::repmat(minimum, N, 1))

    # normalize so each row sums 1
    sum_row <- rowSums(X)
    for (i in 1:N)
    {
      s <- sum_row[i]
      if (s < EPSILON * n) {
        X[i, ] <- ones(1, n) / n
      } else {
        X[i, ] <- X[i, ] / s
      }
    }
  }

  return(X %*% V)
}

principal_component_biplot <- function(X, m = 2, c = 0, k = 1) {
  svd_res <- RSpectra::svds(X, m) # Compute singular value decomposition

  P <- k * svd_res$u %*% diag(svd_res$d^(1 - c)) # Optimal projected points
  V <- svd_res$v %*% diag(svd_res$d^c) / k # Optimal axis vectors

  return(list("P" = P, "V" = V))
}


draw_radial_axes_plot_2d_standardized <- function(
    Z,
    X,
    V,
    P,
    weights = rep(1, ncol(Z)),
    axis_lines = NULL,
    color_variable = NULL,
    projections_point = NULL,
    show_colorbar = FALSE,
    var_names_scale_factor = rep(1, ncol(Z)),
    var_names_horizontal_factor = 0.15,
    aspect_ratio = 1.618) {
  n <- nrow(V)

  # Declare local variables
  w <- x <- y <- P1 <- P2 <- V1 <- V2 <- V3 <- V4 <- NULL


  colnames(P) <- c("P1", "P2")
  variable_names <- colnames(X)
  if (!is.null(color_variable)) {
    color_variable_name <- variable_names[color_variable]
  } else {
    color_variable_name <- NULL
  }

  # Original data means and standard deviations
  meanX <- colMeans(X)
  stdX <- apply(X, 2, stats::sd)

  # Initialize plot
  my_plot <- ggplot2::ggplot()

  # Draw optional reference horizontal & vertical lines
  my_plot <- my_plot +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.5)

  # Draw optional reference circle
  center <- c(0, 0)
  npoints <- 100
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + cos(tt)
  yy <- center[2] + sin(tt)
  dfC <- data.frame(x = xx, y = yy)
  my_plot <- my_plot + ggplot2::geom_path(data = dfC, ggplot2::aes(x, y), linetype = "dotted", color = "gray80", linewidth = 0.5)

  # Set equal axes and empty background
  my_plot <- my_plot + ggplot2::coord_equal() + ggplot2::theme_void()


  # Draw perpendicular dashed segments to axis lines from the point with index "projections_point"
  if (!is.null(projections_point) && !is.null(axis_lines)) {
    p <- as.vector(P[projections_point, ])

    Y <- matrix(0, 0, 4)

    i <- 0
    for (sel_var in axis_lines) {
      i <- i + 1

      v <- V[sel_var, ]

      Y <- rbind(Y, c(p, v %*% ((v %*% p) / norm(v, type = "2")^2)))
    }

    dfY <- as.data.frame(Y)
    my_plot <- my_plot + ggplot2::geom_segment(
      data = dfY, mapping = ggplot2::aes(x = V1, y = V2, xend = V3, yend = V4),
      linetype = "dashed", alpha = 0.7, color = "darkorange", linewidth = 1
    )
  }


  # Draw mapped points on a Cartesian plane
  dfZ <- as.data.frame(Z)
  dfX <- as.data.frame(X)
  dfP <- as.data.frame(P)

  if (is.null(color_variable_name)) {
    dfPZ <- cbind(dfP, dfZ)

    my_plot <- my_plot + ggplot2::geom_point(data = dfPZ, mapping = ggplot2::aes(x = P1, y = P2), color = "#4488FF", size = 2, alpha = 1, show.legend = FALSE)
  } else {
    colorVar <- dfX[, match(color_variable_name, variable_names)]
    dfPZ <- cbind(dfP, dfZ, colorVar)

    my_plot <- my_plot + ggplot2::geom_point(data = dfPZ, mapping = ggplot2::aes(x = P1, y = P2, color = colorVar), size = 2, alpha = 1, show.legend = show_colorbar)

    # Color points using a viridis palette
    my_plot <- my_plot + ggplot2::scale_colour_viridis_c()

    if (show_colorbar) {
      my_plot <- my_plot + ggplot2::theme(legend.key.height = grid::unit(2, "cm")) + ggplot2::theme(
        legend.title = ggplot2::element_text(size = 16),
        legend.text = ggplot2::element_text(size = 16)
      )

      my_plot <- my_plot + ggplot2::labs(colour = color_variable_name)
    }
  }


  # Add a rectangular frame
  my_plot <- my_plot + ggplot2::theme(panel.border = ggplot2::element_rect(fill = "transparent", color = "gray70", linewidth = 2), legend.position = "right")

  # Set margins
  my_plot <- my_plot + ggplot2::theme(plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm"))



  #################### Add calibrated axes ####################

  if (!is.null(axis_lines)) {
    for (sel_var in axis_lines) {
      # Compute data estimates axis of selected variable (sel_var)
      estimates <- P %*% V[sel_var, ]

      # Draw axis line
      lo <- floor(min(estimates))
      up <- pracma::ceil(max(estimates))

      v <- matrix(unlist(V[sel_var, ]), nrow = 1)
      vs <- v / norm(v, type = "2")^2 # scaled axis vector
      E <- cbind(lo * vs, up * vs) # endpoints
      dfE <- as.data.frame(E)
      my_plot <- my_plot + ggplot2::geom_segment(
        data = dfE, mapping = ggplot2::aes(x = V1, y = V2, xend = V3, yend = V4),
        alpha = 0.7, color = "gray30", linewidth = 1
      )

      # Draw tick marks
      size_tick_mark <- 0.1
      v_perp <- v
      v_perp[1, 1] <- -v[1, 2]
      v_perp[1, 2] <- v[1, 1]
      v_perp <- v_perp / norm(v_perp, type = "2") * size_tick_mark / 2

      T <- matrix(0, up - lo + 1, 6)
      i <- 1
      for (k in lo:up) {
        T[i, ] <- cbind(k * vs + v_perp, k * vs - v_perp, k * vs)
        i <- i + 1
      }
      dfT <- as.data.frame(T)
      my_plot <- my_plot + ggplot2::geom_segment(
        data = dfT, mapping = ggplot2::aes(x = V1, y = V2, xend = V3, yend = V4),
        alpha = 0.7, color = "gray30", linewidth = 1
      )

      # Draw numeric values at tick marks
      TL <- cbind(rep(0, up - lo + 1, 1), T[, 5:6])
      colnames(TL) <- c("ticklabel", "x", "y")
      dfTL <- as.data.frame(TL)

      i <- 1
      for (k in lo:up) {
        dfTL[i, 1] <- sprintf("%.3f", stdX[sel_var] * k + meanX[sel_var])
        i <- i + 1
      }

      if (v_perp[1] < 0) {
        v_perp <- -v_perp
      }
      factor <- 3
      x_offset <- v_perp[1] * factor + 0.25
      y_offset <- v_perp[2] * factor
      my_plot <- my_plot + ggplot2::geom_text(
        data = dfTL, mapping = ggplot2::aes(x = x, y = y), label = dfTL$ticklabel, color = "gray30", fontface = "bold",
        nudge_x = x_offset, nudge_y = y_offset, check_overlap = T
      )
    }
  }


  # Draw axis vectors
  colnames(V) <- c("V1", "V2")
  dfV <- as.data.frame(V)

  max_weight <- max(weights)
  vector_colors <- rep(1, n)
  for (i in 1:n) {
    gray_value <- 1 - weights[i] / max_weight
    vector_colors[i] <- grDevices::rgb(gray_value, gray_value, gray_value)
  }

  dfV <- cbind(dfV, vector_colors)

  my_plot <- my_plot + ggplot2::geom_segment(
    data = dfV, mapping = ggplot2::aes(x = 0, y = 0, xend = V1, yend = V2),
    arrow = ggplot2::arrow(angle = 30, length = grid::unit(7, "pt"), type = "closed", ends = "last"),
    lineend = "round", linejoin = "mitre", alpha = 1, color = vector_colors, linewidth = 1
  )


  # Draw variable names
  vertical_offset <- 0.2

  V_names_pos <- diag(var_names_scale_factor) %*% V
  for (i in 1:n) {
    varname_length <- nchar(variable_names[i])

    if (V[i, 2] <= 0) {
      if (V[i, 1] <= 0) {
        V_names_pos[i, ] <- V_names_pos[i, ] + c(-varname_length * var_names_horizontal_factor, -vertical_offset)
      } else {
        V_names_pos[i, ] <- V_names_pos[i, ] + c(varname_length * var_names_horizontal_factor, -vertical_offset)
      }
    } else {
      if (V[i, 1] <= 0) {
        V_names_pos[i, ] <- V_names_pos[i, ] + c(-varname_length * var_names_horizontal_factor, vertical_offset + 0.1)
      } else {
        V_names_pos[i, ] <- V_names_pos[i, ] + c(varname_length * var_names_horizontal_factor, vertical_offset + 0.1)
      }
    }
  }

  dfVN <- data.frame(
    varnames = variable_names,
    x = V_names_pos[, 1],
    y = V_names_pos[, 2]
  )
  my_plot <- my_plot + ggplot2::geom_text(
    data = dfVN, mapping = ggplot2::aes(x = x, y = y),
    label = dfVN$varnames, family = "sans", fontface = "bold",
    color = "black", check_overlap = F, size = 8
  )


  # Set aspect ratio

  x_limits <- layer_scales(my_plot)$x$get_limits()
  y_limits <- layer_scales(my_plot)$y$get_limits()

  x_min <- x_limits[1]
  x_max <- x_limits[2]
  y_min <- y_limits[1]
  y_max <- y_limits[2]


  extend_x <- (aspect_ratio * (y_max - y_min) - (x_max - x_min)) / 2

  if (extend_x > 0) {
    my_plot <- my_plot + xlim(c(x_min - extend_x, x_max + extend_x))
  } else {
    extend_y <- ((x_max - x_min) / aspect_ratio - (y_max - y_min)) / 2

    my_plot <- my_plot + ylim(c(y_min - extend_y, y_max + extend_y))
  }



  # Draw plot
  my_plot
}


draw_radviz <- function(
    X,
    V,
    P,
    color_variable = NULL,
    show_colorbar = FALSE,
    aspect_ratio = 1.618) {
  colnames(P) <- c("P1", "P2")
  variable_names <- colnames(X)
  if (!is.null(color_variable)) {
    color_variable_name <- variable_names[color_variable]
  } else {
    color_variable_name <- NULL
  }

  # Initialize plot
  my_plot <- ggplot2::ggplot()

  # Draw optional reference circle
  center <- c(0, 0)
  npoints <- 100
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + cos(tt)
  yy <- center[2] + sin(tt)
  dfC <- data.frame(x = xx, y = yy)
  my_plot <- my_plot + ggplot2::geom_path(data = dfC, ggplot2::aes(x, y), linetype = "dotted", color = "gray80", linewidth = 0.5)

  # Set equal axes and empty background
  my_plot <- my_plot + ggplot2::coord_equal() + ggplot2::theme_void()


  # Draw anchor points
  dfV <- as.data.frame(V)
  my_plot <- my_plot + geom_point(data = dfV, mapping = aes(x = x, y = y), color = "gray70", size = 6, alpha = 0.5, show.legend = FALSE)

  dfV <- cbind(dfV, rbind(dfV[2:n, ], dfV[1, ]))
  colnames(dfV) <- c("x1", "y1", "x2", "y2")
  my_plot <- my_plot + geom_segment(data = dfV, mapping = aes(x = x1, y = y1, xend = x2, yend = y2), alpha = 1, color = "gray50", linetype = "dashed", linewidth = 1)



  # Draw mapped points on a Cartesian plane
  dfX <- as.data.frame(X)
  dfP <- as.data.frame(P)

  if (is.null(color_variable_name)) {
    dfPZ <- cbind(dfP, dfX)

    my_plot <- my_plot + ggplot2::geom_point(data = dfPZ, mapping = ggplot2::aes(x = P1, y = P2), color = "#4488FF", size = 2, alpha = 1, show.legend = FALSE)
  } else {
    colorVar <- dfX[, match(color_variable_name, variable_names)]
    dfPZ <- cbind(dfP, dfX, colorVar)

    my_plot <- my_plot + ggplot2::geom_point(data = dfPZ, mapping = ggplot2::aes(x = P1, y = P2, color = colorVar), size = 2, alpha = 1, show.legend = show_colorbar)

    # Color points using a viridis palette
    my_plot <- my_plot + ggplot2::scale_colour_viridis_c()

    if (show_colorbar) {
      my_plot <- my_plot + ggplot2::theme(legend.key.height = grid::unit(2, "cm")) + ggplot2::theme(
        legend.title = ggplot2::element_text(size = 16),
        legend.text = ggplot2::element_text(size = 16)
      )

      my_plot <- my_plot + ggplot2::labs(colour = color_variable_name)
    }
  }


  # Add a rectangular frame
  my_plot <- my_plot + ggplot2::theme(panel.border = ggplot2::element_rect(fill = "transparent", color = "gray70", linewidth = 2), legend.position = "right")

  # Set margins
  my_plot <- my_plot + ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 3, 1, "cm"))



  # Draw variable names
  vertical_offset <- 0.2
  var_names_horizontal_factor <- 0.05
  V_names_pos <- V
  for (i in 1:n) {
    varname_length <- nchar(variable_names[i])

    if (V[i, 2] <= 0) {
      if (V[i, 1] <= 0) {
        V_names_pos[i, ] <- V_names_pos[i, ] + c(-varname_length * var_names_horizontal_factor, -vertical_offset)
      } else {
        V_names_pos[i, ] <- V_names_pos[i, ] + c(varname_length * var_names_horizontal_factor, -vertical_offset)
      }
    } else {
      if (V[i, 1] <= 0) {
        V_names_pos[i, ] <- V_names_pos[i, ] + c(-varname_length * var_names_horizontal_factor, vertical_offset + 0.1)
      } else {
        V_names_pos[i, ] <- V_names_pos[i, ] + c(varname_length * var_names_horizontal_factor, vertical_offset + 0.1)
      }
    }
  }

  dfVN <- data.frame(
    varnames = variable_names,
    x = V_names_pos[, 1],
    y = V_names_pos[, 2]
  )
  my_plot <- my_plot + ggplot2::geom_text(
    data = dfVN, mapping = ggplot2::aes(x = x, y = y),
    label = dfVN$varnames, family = "sans", fontface = "bold",
    color = "black", check_overlap = F, size = 8
  )



  # Set aspect ratio

  x_limits <- layer_scales(my_plot)$x$get_limits()
  y_limits <- layer_scales(my_plot)$y$get_limits()

  x_min <- x_limits[1]
  x_max <- x_limits[2]
  y_min <- y_limits[1]
  y_max <- y_limits[2]


  extend_x <- (aspect_ratio * (y_max - y_min) - (x_max - x_min)) / 2

  if (extend_x > 0) {
    my_plot <- my_plot + xlim(c(x_min - extend_x, x_max + extend_x))
  } else {
    extend_y <- ((x_max - x_min) / aspect_ratio - (y_max - y_min)) / 2

    my_plot <- my_plot + ylim(c(y_min - extend_y, y_max + extend_y))
  }


  # Draw plot
  my_plot
}



################################################################################
############################### Figure 1 #######################################
################################################################################

message("\nFigure 1\n")
```

```
## 
## Figure 1
```

``` r
########################  Load data  ###########################################

data("cereal", package = "liver")


# Define subset of (numerical) variables
selected_variables <- c("calories", "protein", "sugars", "vitamins")
data_names <- cereal$name

X <- cereal[, selected_variables] # Select a subset of variables
rm(cereal)


# Remove rows with missing values from X
has_neg <- apply(X, 1, function(row) any(row < 0))
rows_to_delete <- which(has_neg) * (-1)
if (length(rows_to_delete) > 0) {
  X <- X[rows_to_delete, ]
  data_names <- data_names[rows_to_delete]
}


# Convert X to matrix
X <- apply(as.matrix.noquote(X), 2, as.numeric)

# Standardize the data
Z <- scale(X)


########################  Axis vectors for subfigures (a) and (d)  #############

r <- c(1, 0.75, 1, 1.2)
theta <- c(195, 60, 165, 300) * 2 * pi / 360

V_ad <- geometry::pol2cart(theta, r)
rownames(V_ad) <- selected_variables
colnames(V_ad) <- c("V1", "V2")


########################  Subfigure (a) - Star Coordinates  ####################

# Star Coordinates' linear mapping
P <- Z %*% V_ad

# Genearate the SC plot
filename <- "Figure_1a_Cereal_Star_Coordinates.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V_ad,
  P,
  color_variable = 1, # calories
  var_names_scale_factor = c(2.5, 2.2, 2.8, 1.2),
  var_names_horizontal_factor = 0.1,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_1a_Cereal_Star_Coordinates.pdf
```

``` r
########################  Subfigure (b) - RadViz  ##############################

n <- ncol(X)
variable_ordering <- c(2, 3, 4, 1) # "protein", "sugars", "vitamins", "calories"
X_radviz <- X[, variable_ordering]


r <- rep(1, n)
diff_angle <- 2 * pi / n
theta <- (1:n) * diff_angle - diff_angle / 2
V_radviz <- geometry::pol2cart(theta, r)
rownames(V_radviz) <- selected_variables[variable_ordering]

# RadViz (nonlinear) mapping
P_radviz <- radviz(X_radviz, V_radviz, normalize = TRUE)


# Genearate the RadViz plot
filename <- "Figure_1b_Cereal_RadViz.pdf"
pdf(file = filename)
print(draw_radviz(
  X_radviz,
  V_radviz,
  P_radviz,
  color_variable = 1, # calories
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_1b_Cereal_RadViz.pdf
```

``` r
########################  Subfigure (c) - Principal Component Biplot  ##########

PCB <- principal_component_biplot(Z)
P <- PCB$P
V <- PCB$V
colnames(P) <- c("P1", "P2")
colnames(V) <- c("V1", "V2")
rownames(V) <- selected_variables


# Genearate the Principal Component Biplot
filename <- "Figure_1c_Cereal_Biplot.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1:n,
  color_variable = 1, # calories
  projections_point = 4, # All-Bran with Extra Fiber
  show_colorbar = FALSE,
  var_names_scale_factor = c(7, 5, 5, 6),
  var_names_horizontal_factor = 0.1,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_1c_Cereal_Biplot.pdf
```

``` r
########################  Subfigure (d) - ARA Unconstrained L2  ################

mapping <- ara_unconstrained_l2(Z, V_ad)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_1d_Cereal_ARA_unconstrained_l2.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V_ad,
  P,
  axis_lines = 1:n,
  color_variable = 1, # calories
  projections_point = 4, # All-Bran with Extra Fiber
  show_colorbar = FALSE,
  var_names_scale_factor = c(2.2, 3.5, 2.2, 2.1),
  var_names_horizontal_factor = 0.1,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_1d_Cereal_ARA_unconstrained_l2.pdf
```

``` r
################################################################################
#######  Prepare data, axis vectors, and weights for Figures 3, 4, and 5  ######
################################################################################

data("auto_mpg", package = "ascentTraining")


# Define subset of (numerical) variables
selected_variables <- c("mpg", "horsepower", "weight", "acceleration")
data_names <- auto_mpg$car_name

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
if (!is.null(rows_to_delete)) {
  X <- X[rows_to_delete, ]
  data_names <- data_names[rows_to_delete]
}


# Convert X to matrix
X <- apply(as.matrix.noquote(X), 2, as.numeric)

# Standardize the data
Z <- scale(X)

# List of weights (list of ones if w is not used)
w <- c(1, 2 / 3, 2 / 3, 2 / 3)

# Matrix of axis vectors
r <- c(0.8, 1, 1.2, 1)
theta <- c(225, 100, 315, 80) * 2 * pi / 360
V <- geometry::pol2cart(theta, r)
rownames(V) <- selected_variables
colnames(V) <- c("V1", "V2")



################################################################################
############################### Figure 3 #######################################
################################################################################

message("\nFigure 3\n")
```

```
## 
## Figure 3
```

``` r
########################  Subfigure (a) - ARA Unconstrained l2  ################

mapping <- ara_unconstrained_l2(Z, V)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_3a_Auto-mpg_ARA_unconstrained_l2.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_3a_Auto-mpg_ARA_unconstrained_l2.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 976.08289793835
```

``` r
########################  Subfigure (b) - ARA Unconstrained l2 with weights  ###

mapping <- ara_unconstrained_l2(Z, V, weights = w)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_3b_Auto-mpg_ARA_unconstrained_l2_weights.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  weights = w,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_3b_Auto-mpg_ARA_unconstrained_l2_weights.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 999.076511507813
```

``` r
########################  Subfigure (c) - ARA Exact l2  ########################

mapping <- ara_exact_l2(Z, V, variable = 1)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_3c_Auto-mpg_ARA_exact_l2.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_3c_Auto-mpg_ARA_exact_l2.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 1131.23672744631
```

``` r
########################  Subfigure (d) - ARA Ordered l2  ########################

mapping <- ara_ordered_l2(Z, V, variable = 1)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_3d_Auto-mpg_ARA_ordered_l2.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_3d_Auto-mpg_ARA_ordered_l2.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 1045.98513347049
```

``` r
################################################################################
############################### Figure 4 #######################################
################################################################################

message("\nFigure 4\n")
```

```
## 
## Figure 4
```

``` r
########################  Subfigure (a) - ARA Unconstrained l1  ################

mapping <- ara_unconstrained_l1(Z, V)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_4a_Auto-mpg_ARA_unconstrained_l1.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_4a_Auto-mpg_ARA_unconstrained_l1.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 730.891216159581
```

``` r
########################  Subfigure (b) - ARA Unconstrained l1 with weights  ###

mapping <- ara_unconstrained_l1(Z, V, weights = w)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_4b_Auto-mpg_ARA_unconstrained_L1_weights.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  weights = w,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_4b_Auto-mpg_ARA_unconstrained_L1_weights.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 738.191082299998
```

``` r
########################  Subfigure (c) - ARA Exact l1  ########################

mapping <- ara_exact_l1(Z, V, variable = 1)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_4c_Auto-mpg_ARA_exact_L1.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_4c_Auto-mpg_ARA_exact_L1.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 755.976043492109
```

``` r
########################  Subfigure (d) - ARA Ordered l1  ########################

mapping <- ara_ordered_l1(Z, V, variable = 1)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_4d_Auto-mpg_ARA_ordered_L1.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_4d_Auto-mpg_ARA_ordered_L1.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 755.010073537356
```

``` r
################################################################################
############################### Figure 5 #######################################
################################################################################

message("\nFigure 5\n")
```

```
## 
## Figure 5
```

``` r
########################  Subfigure (a) - ARA Unconstrained l-inf  #############

mapping <- ara_unconstrained_linf(Z, V)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_5a_Auto-mpg_ARA_unconstrained_linf.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_5a_Auto-mpg_ARA_unconstrained_linf.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 2.60102201410351
```

``` r
########################  Subfigure (b) - ARA Unconstrained l-inf with weights

mapping <- ara_unconstrained_linf(Z, V, weights = w)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_5b_Auto-mpg_ARA_unconstrained_linf_weights.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  weights = w,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_5b_Auto-mpg_ARA_unconstrained_linf_weights.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 2.74022966229282
```

``` r
########################  Subfigure (c) - ARA Exact l-inf  #####################

mapping <- ara_exact_linf(Z, V, variable = 1)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_5c_Auto-mpg_ARA_exact_linf.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_5c_Auto-mpg_ARA_exact_linf.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 3.72333444576546
```

``` r
########################  Subfigure (d) - ARA Ordered l-inf  ###################

mapping <- ara_ordered_linf(Z, V, variable = 1)
P <- mapping$P
colnames(P) <- c("P1", "P2")

# Genearate the ARA plot
filename <- "Figure_5d_Auto-mpg_ARA_ordered_linf.pdf"
pdf(file = filename)
print(draw_radial_axes_plot_2d_standardized(
  Z,
  X,
  V,
  P,
  axis_lines = 1,
  color_variable = 1, # mpg
  show_colorbar = FALSE,
  var_names_scale_factor = c(2, 1.3, 1.2, 1.1),
  var_names_horizontal_factor = 0.12,
  aspect_ratio = 1.5
))
dev.off()
```

```
## png 
##   2
```

``` r
message("Created subfigure: ", filename, "\n")
```

```
## Created subfigure: Figure_5d_Auto-mpg_ARA_ordered_linf.pdf
```

``` r
message("  Objective value: ", mapping$objval, "\n")
```

```
##   Objective value: 2.63165292624249
```

``` r
sessionInfo()
```

```
## R version 4.5.2 (2025-10-31 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 10 x64 (build 19045)
## 
## Matrix products: default
##   LAPACK version 3.12.1
## 
## locale:
## [1] LC_COLLATE=Spanish_Spain.utf8  LC_CTYPE=Spanish_Spain.utf8    LC_MONETARY=Spanish_Spain.utf8 LC_NUMERIC=C                  
## [5] LC_TIME=Spanish_Spain.utf8    
## 
## time zone: Europe/Madrid
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] ascentTraining_1.0.0 liver_1.24           viridis_0.6.5        viridisLite_0.4.2    RSpectra_0.16-2      geometry_0.5.2      
##  [7] pracma_2.4.6         ggplot2_4.0.1        aramappings_0.1.1    devtools_2.4.6       usethis_3.2.1       
## 
## loaded via a namespace (and not attached):
##  [1] Matrix_1.7-4       gtable_0.3.6       compiler_4.5.2     Rcpp_1.1.0         gridExtra_2.3      scales_1.4.0       clarabel_0.11.1   
##  [8] fastmap_1.2.0      glpkAPI_1.3.4.1    lattice_0.22-7     R6_2.6.1           labeling_0.4.3     knitr_1.50         magic_1.6-1       
## [15] RColorBrewer_1.1-3 rlang_1.1.6        cachem_1.1.0       xfun_0.54          fs_1.6.6           S7_0.2.1           pkgload_1.4.1     
## [22] memoise_2.0.1      cli_3.6.5          withr_3.0.2        magrittr_2.0.4     rstudioapi_0.17.1  remotes_2.5.0      lifecycle_1.0.4   
## [29] vctrs_0.6.5        evaluate_1.0.5     glue_1.8.0         farver_2.1.2       sessioninfo_1.2.3  abind_1.4-8        pkgbuild_1.4.8    
## [36] purrr_1.2.0        tools_4.5.2        ellipsis_0.3.2
```

