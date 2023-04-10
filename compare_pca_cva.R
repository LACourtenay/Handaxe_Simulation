
#

# Code written by Lloyd A. Courtenay #
# Lloyd A. Courtenay - ladc1995@gmail.com (Universidad de Salamanca [USAL]) #

# please ensure that the functions.R file has been sourced prior to running code in this file

#

# Functions ------------------------


compare_pca_cva <- function(shape_params,
                            sigma = 1, sample_sizes = c(30, 30, 30)) {

                                
    create_experimental_dataset <- function(
        shape_params,
        sigma = 1,
        sample_sizes = c(30, 30, 30),
        sample_names = NULL
    ) {
    
    if (length(sample_sizes) != 3) {
        stop("the sample_size vector must be of length 3")
    }
    
    if (!is.null(sample_names)) {
        if (length(sample_names) != 3) {
        stop("the sample_names vector must be of length 3")
        }
    }
    
    shape1 <- create_28_lm_shape(shape_params, param = 1)
    shape2 <- create_28_lm_shape(shape_params, param = 2)
    shape3 <- create_28_lm_shape(shape_params, param = 3)
    
    group_1 <- create_single_sample(shape1, sigma, sample_sizes[1])
    group_2 <- create_single_sample(shape2, sigma, sample_sizes[2])
    group_3 <- create_single_sample(shape3, sigma, sample_sizes[3])
    
    experimental_sample <- abind::abind(group_1, group_2, group_3, along = 3)
    
    if (is.null(sample_names)) {
        experimental_labels <- as.factor(c(
        rep("S1", dim(group_1)[3]),
        rep("S2", dim(group_2)[3]),
        rep("S3", dim(group_3)[3])
        ))
    } else {
        experimental_labels <- as.factor(c(
        rep(sample_names[1], dim(group_1)[3]),
        rep(sample_names[2], dim(group_2)[3]),
        rep(sample_names[3], dim(group_3)[3])
        ))
    }
    
    return(list(
        coords = experimental_sample,
        labels = experimental_labels
    ))
    
    }

  dataset <- create_experimental_dataset(shape_params = shape_params,
                                                sigma = sigma,
                                                sample_sizes = sample_sizes)
  
  GPAshape <- GraphGMM::GPA(dataset$coords)
  
  gmm_pca <- GraphGMM::pca_plot(GraphGMM::vector_from_landmarks(GPAshape$coordinates), dataset$labels,
                                Chull = TRUE, point_size = 4)
  
  gmm_cva <- Morpho::CVA(gmm_pca$pc_scores[,cumsum(gmm_pca$variance) <= 0.95],
                         dataset$labels)
  
  gridExtra::grid.arrange(
    create_plot(GraphGMM::pca_plot(GraphGMM::vector_from_landmarks(GPAshape$coordinates)),
                dataset$labels, "pca", "",
                cv = FALSE, centroid_plot = FALSE),
    create_plot(gmm_cva,
                dataset$labels, "cva", "",
                cv = FALSE, centroid_plot = FALSE),
    ncol = 2
  )
  
}

#

# PCA vs CVA -----------------------

shape_params <- matrix(c(
  1, 0.75, 0.1, 0.45,
  1, 0.75, 0.1, 0.45,
  1, 0.75, 0.1, 0.45
), byrow = FALSE, ncol = 3)

compare_pca_cva(shape_params, sample_sizes = c(20, 20, 20), sigma = 0.5)
compare_pca_cva(shape_params, sample_sizes = c(56, 56, 56), sigma = 0.5)
compare_pca_cva(shape_params, sample_sizes = c(112, 112, 112), sigma = 0.5)

#