
rm(list = ls())
library(RFPCA)
source("data_type1.R")
source("data_type2.R")
source("data_type3.R")
source("data_type4.R")
source("outlier_detect.R")

#----chap1 加在方差上的异常---------
## type1 (乘以k)
para_k_list <- seq(1.5, 5, 0.5) # k = 5

outlier_rate_list <- c(0.04, 0.08, 0.12, 0.16, 0.2)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type1_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type1_result)[1] <- paste("k =", 1)

  for (para_k in para_k_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type1(para_k = para_k, outlier_n = outlier_n)
      result <- outlier_detect(samp = samp_data, outlier_pos = 1:outlier_n)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type1_our/type1_our_para_", para_k,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type1_result <- rbind(type1_result, temp)
    rownames(type1_result)[nrow(type1_result)] <- paste("k =", para_k)
    print(paste("finish the k = ", para_k))
  }

  type1_result_path <- paste0("./result_mean/type1_result_our_rate_", outlier_rate, ".csv")

  # write.csv(type1_result, type1_result_path)
}


# type1_result <- read.csv("./result_mean/type1_result_our.csv")
# plot(type1_result[, 2], type = "o", ylim = c(0, 1))
# plot(type1_result[, 3], type = "o", ylim = c(0, 1))



#----chap2 加在均值mu上的异常---------
## type2 (在mu的第一个维度上乘以k，拉伸)
para_k_list <- seq(1.1, 2, 0.1)

outlier_rate_list <- c(0.04)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type2_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type2_result)[1] <- paste("k =", 1)

  for (para_k in para_k_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type2(para_k = para_k)
      result <- outlier_detect(samp = samp_data, outlier_pos = 1:20)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type2_our/type2_our_para_", para_k,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type2_result <- rbind(type2_result, temp)
    rownames(type2_result)[nrow(type2_result)] <- paste("k =", para_k)
    print(paste("finish the k = ", para_k))
  }

  type2_result_path <- paste0("./result_mean/type2_result_our_rate_", outlier_rate, ".csv")

  # write.csv(type2_result, type2_result_path)
}


# type2_result <- read.csv("./result/type2_result.csv")
# plot(type2_result[, 2], type = "o", ylim = c(0, 1))
# plot(type2_result[, 3], type = "o", ylim = c(0, 1))



## type3（在mu的第一个维度上进行平移）
para_b_list <- seq(0.2, 1.6, 0.2)

outlier_rate_list <- c(0.04)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type3_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type3_result)[1] <- paste("b =", 0)

  for (para_b in para_b_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type3(para_b = para_b)
      result <- outlier_detect(samp = samp_data, outlier_pos = 1:20)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type3_our/type3_our_para_", para_b,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type3_result <- rbind(type3_result, temp)
    rownames(type3_result)[nrow(type3_result)] <- paste("b =", para_b)
    print(paste("finish the b = ", para_b))
  }

  type3_result_path <- paste0("./result_mean/type3_result_our_rate_", outlier_rate, ".csv")

  # write.csv(type3_result, type3_result_path)
}



# type3_result <- read.csv("./result/type3_result.csv")
# plot(type3_result[, 2], type = "o", ylim = c(0, 1))
# plot(type3_result[, 3], type = "o", ylim = c(0, 1))


## type4： 形状变换（在结果中加上偏移量h，）
para_k_list <- seq(1.1, 2, 0.1)

outlier_rate_list <- c(0.04)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type4_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type4_result)[1] <- paste("k =")

  for (para_k in para_k_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type4(para_k = para_k)
      result <- outlier_detect(samp = samp_data, outlier_pos = 1:20)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type4_our/type4_our_para_", para_k,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type4_result <- rbind(type4_result, temp)
    rownames(type4_result)[nrow(type4_result)] <- paste("k =", para_k)
    print(paste("finish the k = ", para_k))
  }

  type4_result_path <- paste0("./result_mean/type4_result_our_rate_", outlier_rate, ".csv")

  # write.csv(type4_result, type4_result_path)
}


# type4_result <- read.csv("./result/type4_result.csv")
# plot(type4_result[, 2], type = "o", ylim = c(0, 1))
# plot(type4_result[, 3], type = "o", ylim = c(0, 1))



#----chap3 模拟与其他方法的比较---------
## 3.1 与JCGS中的方法比较
source("JCGS_detect.R")


source("data_type1.R")
## type1 (乘以k)
para_k_list <- seq(1.5, 5, 0.5) # k = 5

outlier_rate_list <- c(0.04, 0.08, 0.12, 0.16, 0.2)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type1_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type1_result)[1] <- paste("k =", 1)

  for (para_k in para_k_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type1(para_k = para_k, outlier_n = outlier_n)
      result <- JCGS_detect(samp = samp_data, outlier_pos = 1:outlier_n)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type1_JCGS/type1_JCGS_para_", para_k,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type1_result <- rbind(type1_result, temp)
    rownames(type1_result)[nrow(type1_result)] <- paste("k =", para_k)
    print(paste("finish the k = ", para_k))
  }

  type1_result_path <- paste0("./result_mean/type1_result_JCGS_rate_", outlier_rate, ".csv")

  # write.csv(type1_result, type1_result_path)
}


# type1_result_JCGS <- read.csv("./result/type1_result_JCGS.csv")
# plot(type1_result_JCGS[, 2], type = "o", ylim = c(0, 1))
# plot(type1_result_JCGS[, 3], type = "o", ylim = c(0, 1))




source("data_type2.R")
## type2 (在mu的第一个维度上乘以k，拉伸)
para_k_list <- seq(1.1, 2, 0.1)

outlier_rate_list <- c(0.04)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type2_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type2_result)[1] <- paste("k =", 1)

  for (para_k in para_k_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type2(para_k = para_k)
      result <- JCGS_detect(samp = samp_data, outlier_pos = 1:20)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type2_JCGS/type2_JCGS_para_", para_k,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type2_result <- rbind(type2_result, temp)
    rownames(type2_result)[nrow(type2_result)] <- paste("k =", para_k)
    print(paste("finish the k = ", para_k))
  }

  type2_result_path <- paste0("./result_mean/type2_result_JCGS_rate_", outlier_rate, ".csv")

  # write.csv(type2_result, type2_result_path)
}


# type2_result_JCGS <- read.csv("./result/type2_result_JCGS.csv")
# plot(type2_result_JCGS[, 2], type = "o", ylim = c(0, 1))
# plot(type2_result_JCGS[, 3], type = "o", ylim = c(0, 1))





source("data_type3.R")
## type3（在mu的第一个维度上进行平移）
para_b_list <- seq(0.2, 1.6, 0.2)

outlier_rate_list <- c(0.04)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type3_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type3_result)[1] <- paste("b =", 0)

  for (para_b in para_b_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type3(para_b = para_b)
      result <- JCGS_detect(samp = samp_data, outlier_pos = 1:20)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type3_JCGS/type3_JCGS_para_", para_b,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type3_result <- rbind(type3_result, temp)
    rownames(type3_result)[nrow(type3_result)] <- paste("b =", para_b)
    print(paste("finish the b = ", para_b))
  }

  type3_result_path <- paste0("./result_mean/type3_result_JCGS_rate_", outlier_rate, ".csv")

  # write.csv(type3_result, type3_result_path)
}


# type3_result_JCGS <- read.csv("./result/type3_result_JCGS.csv")
# plot(type3_result_JCGS[, 2], type = "o", ylim = c(0, 1))
# plot(type3_result_JCGS[, 3], type = "o", ylim = c(0, 1))




source("data_type4.R")
## type4： 形状变换（在结果中加上偏移量h，）
para_k_list <- seq(1.1, 2, 0.1)

outlier_rate_list <- c(0.04)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type4_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type4_result)[1] <- paste("k =")

  for (para_k in para_k_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type4(para_k = para_k)
      result <- JCGS_detect(samp = samp_data, outlier_pos = 1:20)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type4_JCGS/type4_JCGS_para_", para_k,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type4_result <- rbind(type4_result, temp)
    rownames(type4_result)[nrow(type4_result)] <- paste("k =", para_k)
    print(paste("finish the k = ", para_k))
  }

  type4_result_path <- paste0("./result_mean/type4_result_JCGS_rate_", outlier_rate, ".csv")

  # write.csv(type4_result, type4_result_path)
}



# type4_result_JCGS <- read.csv("./result/type4_result_JCGS.csv")
# plot(type4_result_JCGS[, 2], type = "o", ylim = c(0, 1))
# plot(type4_result_JCGS[, 3], type = "o", ylim = c(0, 1))









## 3.2 与mainfold data方法比较
source("ord_detect.R")

## type1 (乘以k)
para_k_list <- seq(1.5, 5, 0.5) # k = 5

outlier_rate_list <- c(0.04, 0.08, 0.12, 0.16, 0.2)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type1_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type1_result)[1] <- paste("k =", 1)

  for (para_k in para_k_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type1(para_k = para_k, outlier_n = outlier_n)
      result <- ord_detect(samp = samp_data, outlier_pos = 1:outlier_n)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type1_ord/type1_ord_para_", para_k,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type1_result <- rbind(type1_result, temp)
    rownames(type1_result)[nrow(type1_result)] <- paste("k =", para_k)
    print(paste("finish the k = ", para_k))
  }

  type1_result_path <- paste0("./result_mean/type1_result_ord_rate_", outlier_rate, ".csv")

  # write.csv(type1_result, type1_result_path)
}


# type1_result_ord <- read.csv("./result/type1_result_ord.csv")
# plot(type1_result[, 1], type = "o", ylim = c(0, 1))
# plot(type1_result[, 2], type = "o", ylim = c(0, 1))


## type2 (在mu的第一个维度上乘以k，拉伸)
para_k_list <- seq(1.1, 2, 0.1)

outlier_rate_list <- c(0.04)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type2_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type2_result)[1] <- paste("k =", 1)

  for (para_k in para_k_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type2(para_k = para_k)
      result <- ord_detect(samp = samp_data, outlier_pos = 1:20)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type2_ord/type2_ord_para_", para_k,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type2_result <- rbind(type2_result, temp)
    rownames(type2_result)[nrow(type2_result)] <- paste("k =", para_k)
    print(paste("finish the k = ", para_k))
  }

  type2_result_path <- paste0("./result_mean/type2_result_ord_rate_", outlier_rate, ".csv")

  # write.csv(type2_result, type2_result_path)
}



# type2_result_ord <- read.csv("./result/type2_result_ord.csv")
# plot(type2_result_ord[, 2], type = "o", ylim = c(0, 1))
# plot(type2_result_ord[, 3], type = "o", ylim = c(0, 1))


## type3（在mu的第一个维度上进行平移）
para_b_list <- seq(0.2, 1.6, 0.2)

outlier_rate_list <- c(0.04)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type3_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type3_result)[1] <- paste("b =", 0)

  for (para_b in para_b_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type3(para_b = para_b)
      result <- ord_detect(samp = samp_data, outlier_pos = 1:20)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type3_ord/type3_ord_para_", para_b,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type3_result <- rbind(type3_result, temp)
    rownames(type3_result)[nrow(type3_result)] <- paste("b =", para_b)
    print(paste("finish the b = ", para_b))
  }

  type3_result_path <- paste0("./result_mean/type3_result_ord_rate_", outlier_rate, ".csv")

  # write.csv(type3_result, type3_result_path)
}



# type3_result_ord <- read.csv("./result/type3_result_ord.csv")
# plot(type3_result_ord[, 2], type = "o", ylim = c(0, 1))
# plot(type3_result_ord[, 3], type = "o", ylim = c(0, 1))


## type4： 形状变换（在结果中加上偏移量h，）
para_k_list <- seq(1.1, 2, 0.1)

outlier_rate_list <- c(0.04)

for (outlier_rate in outlier_rate_list) {
  print(outlier_rate)
  outlier_n <- 500 * outlier_rate

  type4_result <- data.frame(tpr = 0, fpr = 0)
  rownames(type4_result)[1] <- paste("k =")

  for (para_k in para_k_list) {
    tpr <- c()
    fpr <- c()
    ## 重复实验取均值
    for (i in 1:200) {
      samp_data <- data_type4(para_k = para_k)
      result <- ord_detect(samp = samp_data, outlier_pos = 1:20)

      tpr <- c(tpr, result$tpr)
      fpr <- c(fpr, result$fpr)
    }

    process_result <- cbind(tpr, fpr)
    process_result_path <- paste0(
      "./result_process/type4_ord/type4_ord_para_", para_k,
      "_rate_", outlier_rate, ".csv"
    )

    # write.csv(process_result, process_result_path,
    #   row.names = FALSE
    # )

    temp <- c(mean(tpr), mean(fpr))
    type4_result <- rbind(type4_result, temp)
    rownames(type4_result)[nrow(type4_result)] <- paste("k =", para_k)
    print(paste("finish the k = ", para_k))
  }

  type4_result_path <- paste0("./result_mean/type4_result_ord_rate_", outlier_rate, ".csv")

  # write.csv(type4_result, type4_result_path)
}




# type4_result_ord <- read.csv("./result/type4_result_ord.csv")
# plot(type4_result_ord[, 2], type = "o", ylim = c(0, 1))
# plot(type4_result_ord[, 3], type = "o", ylim = c(0, 1))
