rm(list = ls())

data_pre <- function(type, outlier_rate = 0.04) {
  if (type == 1) {
    para_list <- seq(1, 5, 0.5)
  } else if (type == 2) {
    para_list <- seq(1, 2, 0.1)
  } else if (type == 3) {
    para_list <- seq(0, 1.6, 0.2)
  } else if (type == 4) {
    para_list <- seq(1.1, 2, 0.1)
  }
  result <- read.csv(paste0("./result_mean/type", type, "_result_our_rate_", outlier_rate, ".csv"))
  result_JCGS <- read.csv(paste0("./result_mean/type", type, "_result_JCGS_rate_", outlier_rate, ".csv"))
  result_ord <- read.csv(paste0("./result_mean/type", type, "_result_ord_rate_", outlier_rate, ".csv"))

  if (type == 4) {
    result <- result[2:11, ]
    result_JCGS <- result_JCGS[2:11, ]
    result_ord <- result_ord[2:11, ]
  }

  # TPR
  df <- data.frame(
    para = para_list,
    outlier_detect = result$tpr,
    JCGS_detect = result_JCGS$tpr,
    ord_detect = result_ord$tpr
  )

  # # FPR
  # df <- data.frame(
  #   para = para_list,
  #   outlier_detect = result$fpr,
  #   JCGS_detect = result_JCGS$fpr,
  #   ord_detect = result_ord$fpr
  # )

  return(df)
}

type1_df <- data_pre(type = 1)
type2_df <- data_pre(type = 2)
type3_df <- data_pre(type = 3)
type4_df <- data_pre(type = 4)


plot_fun <- function(type, outlier_rate = 0.04) {
  if (type == 1) {
    y <- type1_df[, -1]
    axis_at <- c(1, 3, 5, 7, 9)
    axis_labels <- c(1, 2, 3, 4, 5)
  } else if (type == 2) {
    y <- type2_df[, -1]
    axis_at <- c(1, 4, 8, 11)
    axis_labels <- c(1, 1.3, 1.7, 2)
  } else if (type == 3) {
    y <- type3_df[, -1]
    axis_at <- c(1, 2, 4, 6, 8, 9)
    axis_labels <- c(0, 0.2, 0.6, 1, 1.4, 1.6)
  } else if (type == 4) {
    y <- type4_df[, -1]
    axis_at <- c(1, 2, 4, 6, 8, 10)
    axis_labels <- c(1.1, 1.2, 1.4, 1.6, 1.8, 2)
  }
  par(mai = c(1.5, 1, 1, 0.5))
  matplot(y,
    type = "o",
    # col = c("black", "blue", "red"),
    col = c("blue", "red", "darkolivegreen3", "#af9a8b"),
    # col = c("black","blue","red","darkolivegreen3","#af9a8b","#594675"),
    lty = "solid",
    lwd = 2.5,
    # ylim = range(y) + c(-.15*sd(y[p,]), .15*sd(y[p,])),
    ylim = c(0, 1),
    ylab = "TPR",
    xlab = "para",
    axes = F, pch = 13:18,
    cex = 1.3, cex.lab = 2,
    col.lab = "gray20"
  )
  grid(col = "grey75", lwd = .3)
  axis(1,
    col = "white", at = axis_at, labels = axis_labels,
    col.ticks = "grey61",
    lwd.ticks = .3, tck = -0.025,
    cex.axis = 1.5, cex.lab = 8, col.axis = "gray30"
  )
  # at=c(0,0.2,0.4,0.6,0.8,1),
  # labels=c(0,0.2,0.4,0.6,0.8,1),
  axis(2,
    col = "white",
    col.ticks = "grey61",
    lwd.ticks = .2, tck = -0.025,
    cex.axis = 1.2, col.axis = "gray30"
  )

  box(col = "grey51")
  # ,"FB-MBD","ADA","OUG"
  # ,"darkolivegreen3","#af9a8b","#594675"
  show_legend <- TRUE
  if (show_legend) {
    legend("topleft",
      legend = c("outlier", "JCGS", "ord"),
      lty = c("solid", "solid", "solid"),
      lwd = c(2, 2),
      col = c("blue", "red", "darkolivegreen3", "#af9a8b"),
      text.col = "gray40", bty = "n",
      box.lwd = .1, pch = 13:18, xjust = -1, cex = 1, inset = .01
    )
  }
  mtext(paste0("TYPE_", type, "_outlier_rate_", outlier_rate), 3,
    adj = 0.5, line = 1, cex = 1.5,
    col = "gray20"
  )
}

plot_fun(type = 1)
plot_fun(type = 2)
plot_fun(type = 3)
plot_fun(type = 4)

type1_df <- data_pre(type = 1, outlier_rate = 0.08)
plot_fun(type = 1, outlier_rate = 0.08)

type1_df <- data_pre(type = 1, outlier_rate = 0.12)
plot_fun(type = 1, outlier_rate = 0.12)

type1_df <- data_pre(type = 1, outlier_rate = 0.16)
plot_fun(type = 1, outlier_rate = 0.16)

type1_df <- data_pre(type = 1, outlier_rate = 0.2)
plot_fun(type = 1, outlier_rate = 0.2)
