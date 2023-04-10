
#

# Code written by Lloyd A. Courtenay #
# Lloyd A. Courtenay - ladc1995@gmail.com (Universidad de Salamanca [USAL]) #

#

# Load libraries and dependencies ----------------------------------------------------

# R libraries by L.A. Courtenay

library(GraphGMM) # for Geometric Morphometrics

# to install these libraries use the install_github function from the devtools library
# library(devtools)
# install_github("LACourtenay/GraphGMM")

# external dependencies

library(geomorph) # for basic geometric morphometric functions
library(shapes) # for basic geometric morphometric functions
library(Morpho) # for the implementation of CVA and bgPCA
library(ggplot2) # for plotting
library(gridExtra) # for plotting
library(abind) # for matrix and tensor operations

#

# Elliptic Fourier Analysis Functions -----------------------------------------------------------------

# function for the superimposition of outlines for EFA (if required)
# coordiante_array is the tensor containing all coordinates delimiting the outlines

superimpose_outlines <- function(coordinate_array) {
  
  superimposed_coords <- coordinate_array
  for (coord in 1:dim(coordinate_array)[3]) {
    
    target_coords <- coordinate_array[,,coord]
    
    # center
    
    target_coords <- apply(target_coords, 2, function(x) x - mean(x))
    
    # align
    
    target_coords <- target_coords %*% svd(var(target_coords))$u
    
    superimposed_coords[,,coord] <- target_coords
    
    
  }
  
  return(superimposed_coords)
  
}

# function for the calculation of EF coefficients from a single individual
# coordinates is the matrix containing the coordinates of the single outline
# n_harmonics are the number of harmonics the user wishes to calculate
# scale indicates whether the coefficients are to be normalised or not

elliptic_fourier_coefficients <- function(coordinates, n_harmonics, scale = FALSE) {
  
  n_points <- nrow(coordinates)
  
  dx <- coordinates[,1] - coordinates[,1][c(n_points, (1:(n_points - 1)))]
  dy <- coordinates[,2] - coordinates[,2][c(n_points, (1:(n_points - 1)))]
  dt <- sqrt(dx^2 + dy^2)
  dt[dt < 1e-10] <- 1e-10  # to avoid Nan - recomendation Momocs
  
  tp <- cumsum(dt) # in momocs tp is called t1
  tp_1 <- c(0, tp[-n_points]) # in momocs tp_1 is called t1m1
  perim_T <- sum(dt) # in momocs this is called T (changed for compatability reasons)
  
  an <- bn <- cn <- dn <- numeric(n_harmonics)
  
  for (harmonic in 1:n_harmonics) {
    Ti <- (perim_T / (2 * pi^2 * harmonic^2))
    numerator <- 2 * harmonic * pi # in momocs this is called r
    an[harmonic] <- Ti * sum(
      (dx / dt) * (cos((numerator * tp) / perim_T) - cos((numerator * tp_1) / perim_T))
    )
    bn[harmonic] <- Ti * sum(
      (dx / dt) * (sin((numerator * tp) / perim_T) - sin((numerator * tp_1) / perim_T))
    )
    cn[harmonic] <- Ti * sum(
      (dy / dt) * (cos((numerator * tp) / perim_T) - cos((numerator * tp_1) / perim_T))
    )
    dn[harmonic] <- Ti * sum(
      (dy / dt) * (sin((numerator * tp) / perim_T) - sin((numerator * tp_1) / perim_T))
    )
  }
  
  if (scale == TRUE) {
    
    a1 <- an[1]
    b1 <- bn[1]
    c1 <- cn[1]
    d1 <- dn[1]
    
    # calculate rotation angle
    
    psi <- 0.5 * atan(
      2 * (a1 * b1 + c1 * d1) /
        (a1^2 + c1^2 - b1^2 - d1^2)
    ) %% pi # originally theta - not sure why the modular division with pi
    
    phaseshift_matrix <- matrix(
      c(cos(psi), sin(psi), -sin(psi), cos(psi)), 2, 2
    ) # originally phaseshift
    
    transposed_matrix <- matrix(c(a1, c1, b1, d1), 2, 2) %*% phaseshift_matrix # originally M2
    
    v <- apply(transposed_matrix^2, 2, sum) # in momocs but i'm not sure why
    if (v[1] < v[2]) {
      psi <- psi + pi/2
    }
    
    psi <- ((psi + pi/2) %% pi) - (pi/2)
    
    a_prima <- a1 * cos(psi) + b1 * sin(psi) # originally Aa
    c_prima <- c1 * cos(psi) + d1 * sin(psi) # originally Cc
    
    Lambda <- sqrt(a_prima^2 + c_prima^2) # originally scale
    
    theta <- atan(c_prima / a_prima) %% pi # originally psi
    
    if (a_prima < 0) {
      theta <- theta + pi
    }
    
    size_coefficient <- 1 / Lambda
    
    rotation_matrix <- matrix(
      c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2
    )
    
    norm_an <- norm_bn <- norm_cn <- norm_dn <- numeric(n_harmonics)
    
    for (harmonic in 1:n_harmonics) {
      normalised_matrix <- size_coefficient * rotation_matrix %*%
        matrix(c(an[harmonic], cn[harmonic], bn[harmonic], dn[harmonic]),
               2, 2) %*%
        matrix(c(
          cos(harmonic * psi), sin(harmonic * psi),
          -sin(harmonic * psi), cos(harmonic * psi)
        ), 2, 2)
      norm_an[harmonic] <- normalised_matrix[1,1]
      norm_bn[harmonic] <- normalised_matrix[1,2]
      norm_cn[harmonic] <- normalised_matrix[2,1]
      norm_dn[harmonic] <- normalised_matrix[2,2]
    }
    
    coe_ind <- c(norm_an, norm_bn, norm_cn, norm_dn)
    
    if (coe_ind[1] < 0) {
      coe_ind <- coe_ind * -1 # recomendation of momocs
    }
    coe_ind[abs(coe_ind) < 1e-12] = 0 # recomendation of momocs
    
  } else {
    
    coe_ind <- c(an, bn, cn, dn)
    
  }
  
  return(coe_ind)
  
}

# function for the calculation of EF coefficients from an entire dataset.
# coordiante_array is the tensor containing all coordinates delimiting the outlines
# n_harmonics are the number of harmonics the user wishes to calculate
# scale indicates whether the coefficients are to be normalised or not

elliptic_fourier_analysis <- function(coordinate_array, n_harmonics, scale = FALSE) {
  
  coefficients <- array(numeric(), dim = c(0, n_harmonics * 4))
  
  for (individual in 1:dim(coordinate_array)[3]) {
    
    #coordinates <- coo[individual][[1]]
    coordinates <- coordinate_array[,,individual]
    
    coe_ind <- elliptic_fourier_coefficients(coordinates,
                                             n_harmonics = n_harmonics,
                                             scale = scale)
    
    coefficients <- abind::abind(coefficients, coe_ind, along = 1)
    
  }
  
  coefficients <- as.data.frame(coefficients)
  coef_names <- c()
  for(coef_val in c("A", "B", "C", "D")) {
    for(harmonic in 1:n_harmonics) {
      coef_names <- c(coef_names, paste0(coef_val, harmonic))
    }
  }
  colnames(coefficients) <- c(coef_names)
  
  return(coefficients)
  
}

# function for the calculation of harmonic power
# an, bn, cn and dn are the coefficients for a given harmonic

harmonic_power <- function(an, bn, cn, dn) {
  power <- (an^2 + bn^2 + cn^2 + dn^2)/2
  return(power)
}

# function to find and calculate the power of the optimal number of harmonics
# coordiante_array is the tensor containing all coordinates delimiting the outlines
# lower_limits and upper_limits mark the lower and upper acceptable threshold of harmonic power.
# i.e., if the user only wants harmonics representing more than 50% of harmonic power, then the
# lower limit should be set to 0.5

calculate_harmonic_power <- function(coordinate_array, lower_limits = 0.95, upper_limits = 0.99) {
  
  find_optimal_harmonic <- function(harmonic_results, target_threshold = 0.95) {
    
    median_power <- apply(harmonic_results, 2, median)
    
    return(which(median_power > target_threshold * 100)[1])
    
  }

  cat("\nCalculating optimal number of harmonics.")
  
  possible_n_harmonics <- floor(dim(coordinate_array)[1] / 2) - 1
  
  results <- array(numeric(), dim = c(0, (possible_n_harmonics - 1)))
  x <- 2:possible_n_harmonics
  
  for (individual in 1:dim(coordinate_array)[3]) {
    x_coefs <- elliptic_fourier_coefficients(coordinate_array[,,individual],
                                             n_harmonics = possible_n_harmonics,
                                             scale = TRUE)
    
    an <- x_coefs[1:possible_n_harmonics]
    bn <- x_coefs[(possible_n_harmonics + 1):(possible_n_harmonics * 2)]
    cn <- x_coefs[((possible_n_harmonics * 2) + 1):(possible_n_harmonics * 3)]
    dn <- x_coefs[((possible_n_harmonics * 3) + 1):(possible_n_harmonics * 4)]
    
    power <- harmonic_power(an, bn, cn, dn)[x]
    results <- abind::abind(results, power, along = 1)
    
  }
  
  results <- t(apply(results, 1, function(pow) {cumsum(pow) / sum(pow)})) * 100
  
  opt_rank <- ceiling(
    (find_optimal_harmonic(results, lower_limits) + find_optimal_harmonic(results, upper_limits)) / 2
  )
  
  decomp_c <- c(); decomp_samp <- c(); for(i in 1:ncol(results)) {
    decomp_c <- c(decomp_c, results[,i])
    decomp_samp <- c(decomp_samp, rep(i, nrow(results)))
  }
  
  power_plot_data <- data.frame(
    Harmonic_Power = decomp_c,
    Harmonic_Rank = as.factor(decomp_samp)
  )
  
  power_value <- median(power_plot_data$Harmonic_Power[power_plot_data$Harmonic_Rank == opt_rank])
  
  cat(paste0("\n", opt_rank, " EFA harmonics has a power of ", power_value, " %"))
  
  return(opt_rank)
  
}

#

# Main Functions - simulating handaxe datasets -----------------------------------------------------------------

# function to simulate a shape using the formulas from DOI: 10.1111/nyas.14680

# x is a single value that should lie between -L/2 and L/2
# L is the length of the geometric icon
# B is the maximum breadth
# w is the distance between two vertical lines corresponding to the maximum breath and y-axis
# D is the diameter of the icon at the distance L/4 from the pointed end

simulate_geometric_icon <- function(x, L, B, w, D) {

  # for computational ease the formula has been broken up into 4 subformulas that are plugged
  # into the final formula
 
  form_1 <- sqrt((5.5 * L^2) + (11 * L * w) + (4 * w^2))
  form_2 <- sqrt(L^2 + (2 * w * L) + (4 * w^2))
  form_3 <- sqrt(3) * B * L
  form_4 <- sqrt(L^2 + (8 * w * x) + (4 * w^2))
  
  value <- (B / 2) * sqrt(
    (L^2 - (4 * x^2)) / form_4
  ) * (
    1 - (
      (form_1 * (form_3 - (2 * D * form_2)) / (form_3 * (form_1 - (2 * form_2))))
    ) * (
      1 - sqrt(
        L * form_4 / (
          (2*(L - 2*w)*x^2 + ((L^2 + (8 * L * w) - (4 * w^2)) * x) + (2 * L * w^2) + (L^2 * w) + L^3)))
    )
  )
  
  return(value)
  
}

# function for the creation of a base geometric icon described using semilandmarks or points around
# the outline
# n_points is the parameter defining the number of semilandmarks or the number of points
# shape_params is a matrix containing the shape parameters that will be plugged into the mathematical
# funciton defining the base geometry
# param indicates the index of which of the shape_params will be used for the simulation

create_base_geometry <- function(n_points, shape_params, param = 1) {
  
  L_length <- shape_params[1,param]
  x <- seq(from = -L_length/2, to = L_length/2, length.out = 200)
  B_breadth <- shape_params[2,param]
  w_dist_max_breadth_and_y <- shape_params[3,param]
  D_diameter_L_4_from_point <- shape_params[4,param]
  
  shape <- simulate_geometric_icon(x, L_length, B_breadth,
                                 w_dist_max_breadth_and_y, D_diameter_L_4_from_point)
  
  x_axis_new <- c(shape, rev(shape * -1)) * 20
  y_axis_new <- c(x, rev(x)) * 20
  
  geometry_line <- sf::st_linestring(as.matrix(cbind(x_axis_new, y_axis_new)))
  
  cat("\nCalculating and extracting semilandmarks.")
  accepted_threshold = 0.01
  threshold = 0.5
  while(abs(threshold) > accepted_threshold) {
    trial_smlm <- sf::st_sample(geometry_line, n_points, type = "regular", exact = TRUE)
    #print(trial[[1]][1,1])
    threshold = trial_smlm[[1]][1,1]
  }
  
  base_geometry <- sf::st_coordinates(trial_smlm)[,1:2]
  base_geometry <- base_geometry[c(((n_points / 2) + 1):n_points, 1:(n_points / 2)),]
  
  return(base_geometry)
  
}

# function for the creation of a single sample of shapes based on the pertubation of landmarks for a single reference shape
# shape is a matrix of coordinate values defining the shape that will be used as reference to create the sample
# sigma is the deviation parameter of the gaussian distribution that is used for landmark deformation
# sample_size is the number of individuals to be created for this sample.
# NOTE: we recommend setting sigma to 1

# function from https://github.com/LACourtenay/GMM_ordination_classification_experimental_toolkit/blob/main/functions.R 

create_single_sample <- function(shape, sigma, sample_size) {
  
  p <- nrow(shape)
  k <- ncol(shape)
  
  single_sample <- array(dim = c(p, k, 0))
  
  for (individual in 1:sample_size) {
    
    x_var <- c()
    y_var <- c()
    
    for (lm in 1:p) {
      
      x_var <- c(x_var, rnorm(1, mean = shape[lm,1], sd = sigma))
      y_var <- c(y_var, rnorm(1, mean = shape[lm,2], sd = sigma))
      
    }
    
    ind_var <- as.matrix(cbind(x_var, y_var), ncol = 2)
    
    single_individual <- shape + ind_var
    
    single_sample <- abind::abind(single_sample, single_individual, along = 3)
    
  }
  
  return(single_sample)
  
}

# function to create the experimental dataset with the different geometric icons for morphological analyses
# n_points describes the number of landmarks or points used to define the oultine
# shape_params is a matrix containing the shape parameters that will be plugged into the mathematical
# funciton defining the base geometry
# sigma is the deviation parameter of the gaussian distribution that is used for landmark deformation
# sample_sizes is a vector of length g that defines the size of each of the samples to be created
# sample_names is a vector of length g that defines the names of the samples (if the user wishes to define sample names)

# function adapted from https://github.com/LACourtenay/GMM_ordination_classification_experimental_toolkit/blob/main/functions.R 

create_experimental_dataset <- function(
    n_points, shape_params,
    sigma = 1,
    sample_sizes = c(30, 30, 30, 30),
    sample_names = NULL
) {
  
  if (length(sample_sizes) != 4) {
    stop("the sample_size vector must be of length 4")
  }
  
  if (!is.null(sample_names)) {
    if (length(sample_names) != 4) {
      stop("the sample_names vector must be of length 4")
    }
  }
  
  shape1 <- create_base_geometry(n_points, shape_params, param = 1)
  shape2 <- create_base_geometry(n_points, shape_params, param = 2)
  shape3 <- create_base_geometry(n_points, shape_params, param = 3)
  shape4 <- create_base_geometry(n_points, shape_params, param = 4)
  
  group_1 <- create_single_sample(shape1, sigma, sample_sizes[1])
  group_2 <- create_single_sample(shape2, sigma, sample_sizes[2])
  group_3 <- create_single_sample(shape3, sigma, sample_sizes[3])
  group_4 <- create_single_sample(shape4, sigma, sample_sizes[4])
  
  experimental_sample <- abind::abind(group_1, group_2, group_3,
                                      group_4, along = 3)
  
  if (is.null(sample_names)) {
    experimental_labels <- as.factor(c(
      rep("S1", dim(group_1)[3]),
      rep("S2", dim(group_2)[3]),
      rep("S3", dim(group_3)[3]),
      rep("S4", dim(group_4)[3])
    ))
  } else {
    experimental_labels <- as.factor(c(
      rep(sample_names[1], dim(group_1)[3]),
      rep(sample_names[2], dim(group_2)[3]),
      rep(sample_names[3], dim(group_3)[3]),
      rep(sample_names[4], dim(group_4)[3])
    ))
  }
  
  return(list(
    coords = experimental_sample,
    labels = experimental_labels
  ))
  
}

# function for the visualisation of shape changes across the PC scores
# coordinates is the tensor of all landmark or outline coordinates used to predict
# shape changes
# pc_scores are the PC scores from the analsis
# gmm - set gmm to TRUE if the coordinates provided are landmark coordinates,
# and to FALSE is the coordinates are from an outline. Setting GMM to TRUE will visualise
# a deformation grid.

visualise_shape_change <- function(coordinates, pc_scores, gmm = TRUE) {
  
  par(mfrow = c(2,2))
  
  for (i in 1:2){
    
    negative_pc <- GraphGMM::morphological_predictor(
      coordinates,
      pc_scores[,i],
      min(pc_scores[,i])
    )
    
    if (gmm == TRUE) {
      
      geomorph::plotRefToTarget(
        GraphGMM::calc_central_morph(coordinates),
        negative_pc
      )
      
    } else {
      
      plot_outline(negative_pc)
      
    }
    
    title(main = paste0("-PC", i))
    
    positive_pc <- GraphGMM::morphological_predictor(
      coordinates,
      pc_scores[,i],
      max(pc_scores[,i])
    )
    
    if (gmm == TRUE) {
      
      geomorph::plotRefToTarget(
        GraphGMM::calc_central_morph(coordinates),
        positive_pc
      )
      
    } else {
      
      plot_outline(positive_pc)
      
    }
    
    title(main = paste0("+PC", i))
    
  }
  
  par(mfrow = c(1,1))
  
}

# simple funcion to plot and visualise the theoretical shapes being included within the study

plot_theoretical_shapes <- function(n_points, shape_params,
                                    sample_names = NULL) {
  
  shape1 <- create_base_geometry(n_points, shape_params, param = 1)
  shape2 <- create_base_geometry(n_points, shape_params, param = 2)
  shape3 <- create_base_geometry(n_points, shape_params, param = 3)
  shape4 <- create_base_geometry(n_points, shape_params, param = 4)
  
  par(mfrow = c(2,2))
  
  plot_outline(shape1)
  if (is.null(sample_names)) {
    title("S1")
  } else {
    title(sample_names[1])
  }
  
  plot_outline(shape2)
  if (is.null(sample_names)) {
    title("S2")
  } else {
    title(sample_names[2])
  }
  
  plot_outline(shape3)
  if (is.null(sample_names)) {
    title("S3")
  } else {
    title(sample_names[3])
  }
  
  plot_outline(shape4)
  if (is.null(sample_names)) {
    title("S4")
  } else {
    title(sample_names[4])
  }
  
  par(mfrow = c(1,1))
  
}

# function for the creation of a base geometric icon described using 28 fixed landmarks
# shape_params is a matrix containing the shape parameters that will be plugged into the mathematical
# funciton defining the base geometry
# param indicates the index of which of the shape_params will be used for the simulation

create_28_lm_shape <- function(shape_params, param = 1) {
  L_length <- shape_params[1,param]
  x <- 1 - c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.80, 0.85, 0.9, 0.95, 1)
  B_breadth <- shape_params[2,param]
  w_dist_max_breadth_and_y <- shape_params[3,param]
  D_diameter_L_4_from_point <- shape_params[4,param]
  
  bordes_y <- simulate_geometric_icon(x - 0.5, L_length, B_breadth,
                                    w_dist_max_breadth_and_y, D_diameter_L_4_from_point)
  bordes_y <- c(bordes_y, rev(bordes_y * -1))
  bordes_x <- c(x, rev(x))
  bordes_shape <- as.matrix(cbind(bordes_y, bordes_x))[-c(1,15),]
  
  return(bordes_shape * 20)
}

# function to create the experimental dataset with the different geometric icons for morphological analyses described
# using 28 fixed landmarks
# shape_params is a matrix containing the shape parameters that will be plugged into the mathematical
# funciton defining the base geometry
# sigma is the deviation parameter of the gaussian distribution that is used for landmark deformation
# sample_sizes is a vector of length g that defines the size of each of the samples to be created
# sample_names is a vector of length g that defines the names of the samples (if the user wishes to define sample names)

create_experimental_28_lm_dataset <- function(
    shape_params,
    sigma = 1,
    sample_sizes = c(30, 30, 30, 30),
    sample_names = NULL
) {
  
  if (length(sample_sizes) != 4) {
    stop("the sample_size vector must be of length 4")
  }
  
  if (!is.null(sample_names)) {
    if (length(sample_names) != 4) {
      stop("the sample_names vector must be of length 4")
    }
  }
  
  shape1 <- create_28_lm_shape(shape_params, param = 1)
  shape2 <- create_28_lm_shape(shape_params, param = 2)
  shape3 <- create_28_lm_shape(shape_params, param = 3)
  shape4 <- create_28_lm_shape(shape_params, param = 4)
  
  group_1 <- create_single_sample(shape1, sigma, sample_sizes[1])
  group_2 <- create_single_sample(shape2, sigma, sample_sizes[2])
  group_3 <- create_single_sample(shape3, sigma, sample_sizes[3])
  group_4 <- create_single_sample(shape4, sigma, sample_sizes[4])
  
  experimental_sample <- abind::abind(group_1, group_2, group_3,
                                      group_4, along = 3)
  
  if (is.null(sample_names)) {
    experimental_labels <- as.factor(c(
      rep("S1", dim(group_1)[3]),
      rep("S2", dim(group_2)[3]),
      rep("S3", dim(group_3)[3]),
      rep("S4", dim(group_4)[3])
    ))
  } else {
    experimental_labels <- as.factor(c(
      rep(sample_names[1], dim(group_1)[3]),
      rep(sample_names[2], dim(group_2)[3]),
      rep(sample_names[3], dim(group_3)[3]),
      rep(sample_names[4], dim(group_4)[3])
    ))
  }
  
  return(list(
    coords = experimental_sample,
    labels = experimental_labels
  ))
  
}

#

# Miscellaneous Functions -----------------------------------------------------------------

# Function for setting the alpha of colours

add_alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

# a function that can be used to create a plot of either a bgPCA, a CVA or a PCA
# data is the object produced by either the groupPCA, CVA or pca_plot functions.
# labels are the list of factor labels for the sample
# type is either "bgpca", "cva" or "pca"
# main is the title of the plot
# cv is a boolean value that indicates whether a Cv-bgPCA or Cv-CVA is to be visualised
# centroid_plot is a boolean parameter that defines whether the entire scatter plot is to be visualised
#       or only the centroid of distributions

# code originally from:
# https://github.com/LACourtenay/GMM_ordination_classification_experimental_toolkit/blob/main/functions.R

create_plot <- function(data, labels, type, main, cv, centroid_plot) {
  
  if (type == "bgpca") {
    
    var1 <- data$Variance[1,2][[1]] * 100
    var2 <- data$Variance[2,2][[1]] * 100
    
    if (cv == TRUE) {
      data <- data$CV
      xlabel <- ggplot2::xlab(paste("Cv bgPC1 (",
                                    round(var1,2),"%)", sep = ""))
      ylabel <- ggplot2::ylab(paste("Cv bgPC2 (",
                                    round(var2,2),"%)", sep = ""))
    } else {
      data <- data$Scores
      xlabel <- ggplot2::xlab(paste("bgPC1 (",
                                    round(var1,2),"%)", sep = ""))
      ylabel <- ggplot2::ylab(paste("bgPC2 (",
                                    round(var2,2),"%)", sep = ""))
    }
    
    
  } else if (type == "cva") {
    
    var1 <- data$Var[1,2][[1]]
    var2 <- data$Var[2,2][[1]]
    
    if (cv == TRUE) {
      data <- data$CVcv
      xlabel <- ggplot2::xlab(paste("Cv CV1 (",
                                    round(var1,2),"%)", sep = ""))
      ylabel <- ggplot2::ylab(paste("Cv CV2 (",
                                    round(var2,2),"%)", sep = ""))
    } else {
      data <- data$CVscores
      xlabel <- ggplot2::xlab(paste("CV1 (",
                                    round(var1,2),"%)", sep = ""))
      ylabel <- ggplot2::ylab(paste("CV2 (",
                                    round(var2,2),"%)", sep = ""))
    }
    
  } else {
    
    var1 <- data$variance[1] * 100
    var2 <- data$variance[2] * 100
    xlabel <- ggplot2::xlab(paste0("PC1 (", round(var1, 2), "%)"))
    ylabel <- ggplot2::ylab(paste0("PC2 (", round(var2, 2), "%)"))
    
    data <- data$pc_scores[,1:2]
    
  }
  
  data <- data.frame(
    x = data[,1],
    y = data[,2],
    Sample = as.factor(labels)
  )
  
  if (centroid_plot == FALSE) {
    base_plot <- ggplot2::ggplot(data = data,
                                 ggplot2::aes(x = x, y = y, colour = Sample))
    plot_colours <- ggplot2::scale_color_manual(values = c("black","red","blue",
                                                           "orange","darkgreen","violetred"))
    final_plot <- base_plot +
      ggplot2::geom_point(stat = "identity", size = 4) +
      plot_colours
    
    final_plot <- final_plot +
      ggplot2::stat_ellipse(size = 1) +
      xlabel + ylabel +
      ggplot2::ggtitle(main) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
        plot.title = ggplot2::element_text(face = "bold", size = 20),
        plot.subtitle = ggplot2::element_text(size = 15),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA),
        axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                             margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
        axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                             margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 13),
        legend.title = ggplot2::element_text(size = 18, face = "bold"),
        legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
        legend.box.background = ggplot2::element_rect(colour = "black"),
        legend.position = "bottom"
      ) +
      ggplot2::geom_vline(xintercept = 0,
                          colour = "black",
                          size = 0.5,
                          linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 0,
                          colour = "black",
                          linetype = "dashed",
                          size = 0.5)
  } else {
    
    centroids <- array(dim = c(3, 2)); target_index <- 1; for (target_sample in levels(data$Sample)) {
      centroids[target_index,] <- c(
        mean(data[data$Sample == target_sample, 1]), mean(data[data$Sample == target_sample, 2])
      )
      target_index <- target_index + 1
    }
    xlim <- c(min(data[,1]), max(data[,1]))
    ylim <- c(min(data[,1]), max(data[,1]))
    
    final_plot <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::geom_point(ggplot2::aes(x = centroids[1,1], y = centroids[1,2]), size = 7,
                          color = "black") +
      ggplot2::geom_point(ggplot2::aes(x = centroids[2,1], y = centroids[2,2]), size = 7,
                          color = "red") +
      ggplot2::geom_point(ggplot2::aes(x = centroids[3,1], y = centroids[3,2]), size = 7,
                          color = "blue") +
      xlabel + ylabel + ggplot2::ggtitle(main) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
        plot.title = ggplot2::element_text(face = "bold", size = 20),
        plot.subtitle = ggplot2::element_text(size = 15),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA),
        axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                             margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
        axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                             margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 13),
        legend.title = ggplot2::element_text(size = 18, face = "bold"),
        legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
        legend.box.background = ggplot2::element_rect(colour = "black"),
        legend.position = "bottom"
      ) +
      ggplot2::geom_vline(xintercept = 0,
                          colour = "black",
                          size = 0.5,
                          linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 0,
                          colour = "black",
                          linetype = "dashed",
                          size = 0.5)
    
  }
  
  
  return(final_plot)
  
}

# function to plot an outline

plot_outline <- function(coordinates) {
  
  coordinates <- rbind(
    coordinates,
    coordinates[1,]
  )
  
  plot(coordinates[,1], coordinates[,2], asp = 1, col = NULL, bty = "n",
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  lines(coordinates[,1], coordinates[,2], lwd = 3)
  
}