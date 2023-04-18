rm(list = ls())
load("./data/3d_plot_data.Rdata")

library(plot3D)
library(rgl)

# 绘制球面
M <- mesh(seq(-pi, pi, length = 50), seq(-pi, pi, length = 50))
u <- M$x
v <- M$y
r <- 1
x <- r * cos(u) * cos(v)
y <- r * cos(u) * sin(v)
z <- r * sin(u)
plot3d(x, y, z, col = "blue", alpha = 0.5, type = "l")

# type = 'p':点，'l':线，'s':球体，'h':线段，'n':无

# # 正常样本
# samp_num_list <- c(4, 5, 6, 7) # 取4个正常样本
# for (samp_num in samp_num_list) {
#   for (i in 1:51) {
#     points3d(normal_data[samp_num, 1, i],
#       normal_data[samp_num, 2, i],
#       normal_data[samp_num, 3, i],
#       col = "black", size = 5
#     )
#   }
# }

# mean
for (i in 1:51) {
  points3d(mu[1, i], mu[2, i], mu[3, i],
    col = "black", size = 3
  )
}

lines3d(mu[1, ], mu[2, ], mu[3, ],
  col = "black", lwd = 1
)


# type1
sample_num_list <- c(1, 2, 13)
for (sample_num in sample_num_list) {
  # ponit
  for (i in 1:51) {
    points3d(outlier_1[sample_num, 1, i],
      outlier_1[sample_num, 2, i],
      outlier_1[sample_num, 3, i],
      col = "red", size = 3
    )
    # Sys.sleep(0.25)  #  依次缓慢描点
  }
  # line
  lines3d(outlier_1[sample_num, 1, ],
    outlier_1[sample_num, 2, ],
    outlier_1[sample_num, 3, ],
    col = "red", lwd = 1
  )
}

# type2
sample_num_list <- c(2, 3, 4)
# sample_num_list <- c(2)

for (sample_num in sample_num_list) {
  # point
  for (i in 1:51) {
    points3d(outlier_2[sample_num, 1, i],
      outlier_2[sample_num, 2, i],
      outlier_2[sample_num, 3, i],
      col = "brown1", size = 3
    )
    # Sys.sleep(0.25)
  }
  # line
  lines3d(outlier_2[sample_num, 1, ],
    outlier_2[sample_num, 2, ],
    outlier_2[sample_num, 3, ],
    col = "brown1", lwd = 1
  )
}

# type3
sample_num_list <- c(5, 7, 8)
# sample_num_list <- c(2)

for (sample_num in sample_num_list) {
  # point
  for (i in 1:51) {
    points3d(outlier_3[sample_num, 1, i],
      outlier_3[sample_num, 2, i],
      outlier_3[sample_num, 3, i],
      col = "brown1", size = 3
    )
    # Sys.sleep(0.25)
  }
  # line
  lines3d(outlier_3[sample_num, 1, ],
    outlier_3[sample_num, 2, ],
    outlier_3[sample_num, 3, ],
    col = "brown1", lwd = 1
  )
}

# type4
sample_num_list <- c(8, 9, 17)
# sample_num_list <- c(2)

for (sample_num in sample_num_list) {
  # point
  for (i in 1:51) {
    points3d(outlier_4[sample_num, 1, i],
      outlier_4[sample_num, 2, i],
      outlier_4[sample_num, 3, i],
      col = "brown1", size = 3
    )
    # Sys.sleep(0.25)
  }
  # line
  lines3d(outlier_4[sample_num, 1, ],
    outlier_4[sample_num, 2, ],
    outlier_4[sample_num, 3, ],
    col = "brown1", lwd = 1
  )
}


# ---------------------------------
plot3d(x, y, z, col = "blue", alpha = 0.5, type = "l")

short_mu <- mu[, 1:20]
# mean
for (i in 1:20) {
  points3d(short_mu[1, i], short_mu[2, i], short_mu[3, i],
    col = "black", size = 3
  )
}

lines3d(short_mu[1, ], short_mu[2, ], short_mu[3, ],
  col = "black", lwd = 1
)

short_normal <- normal_data[, , 1:20]
samp_num_list <- c(4:15) # 取4个正常样本
for (samp_num in samp_num_list) {
  for (i in 1:20) {
    points3d(normal_data[samp_num, 1, i],
      normal_data[samp_num, 2, i],
      normal_data[samp_num, 3, i],
      col = "red", size = 3
    )
  }

  lines3d(short_normal[samp_num, 1, ],
    short_normal[samp_num, 2, ],
    short_normal[samp_num, 3, ],
    col = "red", lwd = 1
  )
}

#-------------------------------
# type1
plot3d(x, y, z, col = "blue", alpha = 0.5, type = "l")
short_type1 <- outlier_1[, , 1:20]
sample_num_list <- c(3, 10, 11) # (10, 14, 15)
for (sample_num in sample_num_list) {
  # ponit
  for (i in 1:20) {
    points3d(short_type1[sample_num, 1, i],
      short_type1[sample_num, 2, i],
      short_type1[sample_num, 3, i],
      col = "brown1", size = 3
    )
    # Sys.sleep(0.25)  #  依次缓慢描点
  }
  # line
  lines3d(short_type1[sample_num, 1, ],
    short_type1[sample_num, 2, ],
    short_type1[sample_num, 3, ],
    col = "brown1", lwd = 1
  )
}


#--------------------------------------
# library("RFPCA")
# library("MASS")
# source("./function/MakephiS.R")
# source("./function/MakeRotMat.R")
# source("./function/MakeSphericalProcess.R")
# source("./function/MakeSphericalOutlier.R")
# source("./function/MakeSphericalOutlier_2.R")
# source("./function/Proj.R")
#
# muList <- list(
#   function(x) x * 2,
#   function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
#   function(x) rep(0, length(x))
# )
# n <- 500 # 样本量
# outlier_n <- 20 # 异常值的数量
# m <- 51
# pts <- seq(0, 1, length.out = m)
# mfd <- structure(1, class = "Sphere")
# mu <- Makemu(mfd, muList, c(rep(0, 2), 1), pts) # 生成均值mu
#
# K <- 20
# sigma2 <- 0
# # set.seed(5)
# # 生成样本数据
# samp <- MakeSphericalProcess(n, mu, pts,
#   K = K,
#   lambda = 0.07^(seq_len(K) / 2),
#   basisType = "legendre01",
#   sigma2 = sigma2
# )
#
# library(plot3D)
# library(rgl)
#
# M <- mesh(seq(-pi, pi, length = 50), seq(-pi, pi, length = 50))
# u <- M$x
# v <- M$y
# r <- 1
# x <- r * cos(u) * cos(v)
# y <- r * cos(u) * sin(v)
# z <- r * sin(u)
# plot3d(x, y, z, col = "blue", alpha = 0.5, type = "l")
#
# for (samp_num in 4:7) {
#   for (i in 1:m) {
#     points3d(samp$X[samp_num, 1, i], samp$X[samp_num, 2, i], samp$X[samp_num, 3, i],
#       col = "black", size = 5
#     )
#   }
# }
#
# source("data_type1.R")
# outlier_1 <- data_type1(para_k = 10, outlier_n = 20)$X[1:20, , ]
#
# sample_num <- 7
# for (i in 1:m) {
#   points3d(outlier_1[sample_num, 1, i], outlier_1[sample_num, 2, i], outlier_1[sample_num, 3, i],
#     col = "red", size = 5
#   )
# }
#
# source("data_type2.R")
# outlier_2 <- data_type2(para_k = 2.5)$X[1:20, , ]
#
# sample_num <- 7
# for (i in 1:m) {
#   points3d(outlier_2[sample_num, 1, i], outlier_2[sample_num, 2, i], outlier_2[sample_num, 3, i],
#     col = "gold", size = 5
#   )
# }
#
# source("data_type3.R")
# outlier_3 <- data_type3(para_b = 1.6)$X[1:20, , ]
# sample_num <- 7
# for (i in 1:m) {
#   points3d(outlier_3[sample_num, 1, i], outlier_3[sample_num, 2, i], outlier_3[sample_num, 3, i],
#     col = "maroon1", size = 5
#   )
# }
#
# source("data_type4.R")
# outlier_4 <- data_type4(para_k = 2)$X[1:20, , ]
# sample_num <- 15
# for (i in 1:m) {
#   points3d(outlier_4[sample_num, 1, i], outlier_4[sample_num, 2, i], outlier_4[sample_num, 3, i],
#     col = "springgreen", size = 5
#   )
# }
#
# normal_data <- samp$X
#
# save(normal_data, outlier_1, outlier_2, outlier_3, outlier_4, mu,
#   file = "./3d_plot_data/3d_plot_data.Rdata"
# )
