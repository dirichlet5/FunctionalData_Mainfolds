library("gdata")
library("fda")
library("compiler")
library("fdapace")
library("fdaoutlier")

source("./function/JCGS/simulation-function.R")
source("./function/JCGS/simulation_parameter.R")
source("./function/JCGS/one_sample.R")
source("./function/JCGS/multi_sparse.R")
source("./function/JCGS/outl_detection.R")

JCGS_detect <- function(samp, outlier_pos) {
  bootstrapTimes <- 10
  confid_alpha <- 0.05
  argval <- seq(0, 1, length.out = dim(samp$X[, 1, ])[2])
  out_set_1 <- c()
  out_set_2 <- c()
  out_set_3 <- c()

  # dim 1
  try(
    {
      result_bs <- bootstrap_implementation(
        bootstrapTimes, argval,
        list(samp$X[, 1, ]),
        list(samp$X[, 1, ]), confid_alpha
      )

      bmfpca_list <- lapply(1:length(argval), function(k) {
        sapply(result_bs$bootstrap_fit, function(l) {
          l[, k]
        })
      })

      depth_bmfpc <- multi_depth(bmfpca_list, "MFHD")

      fit <- as.data.frame(result_bs$bootstrap_fit)
      sp <- sparse_fbplot(
        fit = fit,
        sparse = samp$X[, 1, ],
        time_index = NULL, depth = depth_bmfpc,
        two_stage = FALSE, sq_vo = sq_vo, plot = FALSE,
        xlab = NULL, ylab = NULL, title = NULL,
        yrange = NULL,
        cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3,
        medlabel = TRUE, outlabel = TRUE,
        prob = 0.5, factor = 1.5,
        color = 6, outliercolor.fb = 2, barcol = 4,
        outliercolor.dir = 3, fullout = FALSE
      )

      out_set_1 <- unique(unlist(sp$fb_outlier_index))
    },
    silent = TRUE
  )


  # dim 2
  try(
    {
      result_bs <- bootstrap_implementation(
        bootstrapTimes, argval,
        list(samp$X[, 2, ]),
        list(samp$X[, 2, ]), confid_alpha
      )

      bmfpca_list <- lapply(1:length(argval), function(k) {
        sapply(result_bs$bootstrap_fit, function(l) {
          l[, k]
        })
      })

      depth_bmfpc <- multi_depth(bmfpca_list, "MFHD")

      fit <- as.data.frame(result_bs$bootstrap_fit)
      sp <- sparse_fbplot(
        fit = fit,
        sparse = samp$X[, 2, ],
        time_index = NULL, depth = depth_bmfpc,
        two_stage = FALSE, sq_vo = sq_vo, plot = FALSE,
        xlab = NULL, ylab = NULL, title = NULL,
        yrange = NULL,
        cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3,
        medlabel = TRUE, outlabel = TRUE,
        prob = 0.5, factor = 1.5,
        color = 6, outliercolor.fb = 2, barcol = 4,
        outliercolor.dir = 3, fullout = FALSE
      )

      out_set_2 <- unique(unlist(sp$fb_outlier_index))
    },
    silent = TRUE
  )

  # dim 3
  try(
    {
      result_bs <- bootstrap_implementation(
        bootstrapTimes, argval,
        list(samp$X[, 3, ]),
        list(samp$X[, 3, ]), confid_alpha
      )

      bmfpca_list <- lapply(1:length(argval), function(k) {
        sapply(result_bs$bootstrap_fit, function(l) {
          l[, k]
        })
      })

      depth_bmfpc <- multi_depth(bmfpca_list, "MFHD")

      fit <- as.data.frame(result_bs$bootstrap_fit)
      sp <- sparse_fbplot(
        fit = fit,
        sparse = samp$X[, 3, ],
        time_index = NULL, depth = depth_bmfpc,
        two_stage = FALSE, sq_vo = sq_vo, plot = FALSE,
        xlab = NULL, ylab = NULL, title = NULL,
        yrange = NULL,
        cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3,
        medlabel = TRUE, outlabel = TRUE,
        prob = 0.5, factor = 1.5,
        color = 6, outliercolor.fb = 2, barcol = 4,
        outliercolor.dir = 3, fullout = FALSE
      )

      out_set_3 <- unique(unlist(sp$fb_outlier_index))
    },
    silent = TRUE
  )

  out_set <- unique(c(out_set_1, out_set_2, out_set_3))

  source("./function/evaluation.R")

  tpr <- evaluation(n = n, results = out_set, true = outlier_pos)$tp
  fpr <- evaluation(n = n, results = out_set, true = outlier_pos)$fp
  list(tpr = tpr, fpr = fpr)
}
