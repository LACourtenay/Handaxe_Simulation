
#

# Code written by Lloyd A. Courtenay #
# Lloyd A. Courtenay - ladc1995@gmail.com (Universidad de Salamanca [USAL]) #

# please ensure that the functions.R file has been sourced prior to running code in this file

#

# Define the 4 geometric samples ----------------------------------------------------------

# define shape parameters for the 4 geometric icons

shape_params <- matrix(c(
  1, 0.8, 0.1, 0.57,
  1, 0.75, 0.1, 0.45,
  1, 0.65, 0.15, 0.35,
  1, 0.58, 0.24, 0.3
), byrow = FALSE, ncol = 4)

plot_theoretical_shapes(1000, shape_params, sample_names = c("O", "C", "sT", "T"))

#

# GMM - 28 Fixed Landmark Analysis ----------------------------------------------------------

dataset <- create_28_lm_dataset(shape_params,
                                sample_names = c("O", "C", "sT", "T"),
                                sample_sizes = c(56, 56, 56, 56))

# Generalised Procrustes Analysis

gpashape <- GPA(dataset$coords)

# principal component analysis

pc_analysis <- pca_plot(vector_from_landmarks(gpashape$coordinates), dataset$labels,
                        CI_ellipse = TRUE, point_size = 4, main = "GMM",
                        label_colours = c("blue", "black", "red", "orange"))

pc_analysis$pca_plot

visualise_shape_change(gpashape$coordinates, pc_analysis$pc_scores, gmm = TRUE)

#

# GMM - semilandmark approach (with sliding using bending energy) ----------------------------

dataset <- create_experimental_dataset(28, shape_params,
                                       sample_names = c("O", "C", "sT", "T"),
                                       sample_sizes = c(56, 56, 56, 56))

# sliding semilandmarks

sliding_smlm <- rbind(geomorph::define.sliders(1:((28 / 2) + 1)),
                      geomorph::define.sliders(((28 / 2) + 1):28),
                      rbind(c(
                        28 - 1, 28, 1
                      )))

# Generalised Procrustes Analysis

gpashape <- gpagen(dataset$coords, ProcD = FALSE, curves = sliding_smlm)

# principal component analysis

pc_analysis <- pca_plot(vector_from_landmarks(gpashape$coords), dataset$labels,
                        CI_ellipse = TRUE, point_size = 4, main = "GMM",
                        label_colours = c("blue", "black", "red", "orange"))

pc_analysis$pca_plot

visualise_shape_change(gpashape$coords, pc_analysis$pc_scores, gmm = TRUE)

#

# GMM - semilandmark approach (with sliding using procrustes distances) ----------------------------

dataset <- create_experimental_dataset(28, shape_params,
                                       sample_names = c("O", "C", "sT", "T"),
                                       sample_sizes = c(56, 56, 56, 56))

# sliding semilandmarks

sliding_smlm <- rbind(geomorph::define.sliders(1:((28 / 2) + 1)),
                      geomorph::define.sliders(((28 / 2) + 1):28),
                      rbind(c(
                        28 - 1, 28, 1
                      )))

# Generalised Procrustes Analysis

gpashape <- gpagen(dataset$coords, ProcD = TRUE, curves = sliding_smlm)

# principal component analysis

pc_analysis <- pca_plot(vector_from_landmarks(gpashape$coords), dataset$labels,
                        CI_ellipse = TRUE, point_size = 4, main = "GMM",
                        label_colours = c("blue", "black", "red", "orange"))

pc_analysis$pca_plot

visualise_shape_change(gpashape$coords, pc_analysis$pc_scores, gmm = TRUE)

#

# GMM - semilandmark approach (without sliding) ----------------------------

dataset <- create_experimental_dataset(28, shape_params,
                                       sample_names = c("O", "C", "sT", "T"),
                                       sample_sizes = c(56, 56, 56, 56))

# Generalised Procrustes Analysis

gpashape <- gpagen(dataset$coords)

# principal component analysis

pc_analysis <- pca_plot(vector_from_landmarks(gpashape$coords), dataset$labels,
                        CI_ellipse = TRUE, point_size = 4, main = "GMM",
                        label_colours = c("blue", "black", "red", "orange"))

pc_analysis$pca_plot

visualise_shape_change(gpashape$coords, pc_analysis$pc_scores, gmm = TRUE)

#

# EFA ------------------------------------------------------------------------------------

dataset <- create_experimental_dataset(28, shape_params,
                                       sample_names = c("O", "C", "sT", "T"),
                                       sample_sizes = c(29, 29, 29, 29))

efa_coords <- superimpose_outlines(dataset$coords)
n_harmonics <- calculate_harmonic_power(efa_coords)
coefs_fourier <- elliptic_fourier_analysis(
  efa_coords,
  n_harmonics = n_harmonics,
  scale = TRUE
)

pc_analysis <- GraphGMM::pca_plot(coefs_fourier,
                                  dataset$labels,
                                  CI_ellipse = TRUE, point_size = 4, main = "EFA",
                                  label_colours = c("blue", "black", "red", "orange"))

pc_analysis$pca_plot

visualise_shape_change(efa_coords, pc_analysis$pc_scores, gmm = FALSE)

#