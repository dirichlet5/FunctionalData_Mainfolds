rm(list = ls())

library(RFPCA)
source("data_type1.R")

load(file = "./data/flightaware_data/CPA250_data.RData")
load(file = "./data/flightaware_data/CES7551_data.RData")

df <- array(c(
  aperm(CPA250_data, c(3, 2, 1)),
  aperm(CES7551_data, c(3, 2, 1))
),
dim = c(3, 51, 146)
)

source("outlier_detect.R")
tpr <- c()
fpr <- c()
for (i in 1:100) {
  sample_all <- c(137, 138, 133, 125, 131, 143, 130, 142, 135, 127, 128)
  sample_all <- c(125, 126, 127, 128, 130, 131, 132, 139, 140, 141, 143, 146)
  sample_num <- sample(sample_all, 5, replace = FALSE)
  length_m <- 51
  df_two <- df[, 1:length_m, c(1:124, sample_num)]

  samp_data_X <- aperm(df_two, c(3, 1, 2)) # 交换维度

  samp_data <- list(X = samp_data_X, T = seq(0, 1, length = length_m))


  result <- outlier_detect(samp_data, outlier_pos = 125:129)
  tpr <- c(tpr, result$tpr)
  fpr <- c(fpr, result$fpr)
  # print(result)
  if (i %% 10 == 0) {
    print(paste0("finish:", i))
  }
}
print(tpr)
mean(tpr)
print(fpr)
flight_result <- cbind(tpr, fpr)
flight_result_path <- paste0("./result_flight/outlier_detect", ".csv")

# write.csv(flight_result, flight_result_path,
#           row.names = FALSE)


#------------------------------
source("JCGS_detect.R")
tpr <- c()
fpr <- c()
for (i in 1:100) {
  sample_all <- c(125, 126, 127, 128, 130, 131, 132, 139, 140, 141, 143, 146)
  sample_num <- sample(sample_all, 5, replace = FALSE)
  length_m <- 51
  df_two <- df[, 1:length_m, c(1:124, sample_num)]

  samp_data_X <- aperm(df_two, c(3, 1, 2)) # 交换维度

  samp_data <- list(X = samp_data_X, T = seq(0, 1, length = length_m))


  result <- JCGS_detect(samp_data, outlier_pos = 125:129)
  tpr <- c(tpr, result$tpr)
  fpr <- c(fpr, result$fpr)
  # print(result)
  if (i %% 10 == 0) {
    print(paste0("finish:", i))
  }
}
print(tpr)
print(fpr)
flight_result <- cbind(tpr, fpr)
flight_result_path <- paste0("./result_flight/JCGS_detect", ".csv")

# write.csv(flight_result, flight_result_path,
#           row.names = FALSE)



# box-plot
library(ggplot2)
df1 <- read.csv("./result_flight/JCGS_detect.csv")
df2 <- read.csv("./result_flight/outlier_detect.csv")
# df3 = read.csv("./result_flight2/ord_detect.csv")

par(mfrow = c(1, 2))
boxplot(df2$tpr, df1$tpr,
  names = c("outlier", "JCGS"),
  main = "TPR"
)

boxplot(df2$fpr, df1$fpr,
  names = c("outlier", "JCGS"),
  main = "FPR"
)
