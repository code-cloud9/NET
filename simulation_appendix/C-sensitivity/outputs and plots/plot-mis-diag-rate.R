library(ggplot2)

c_list <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6)

a_alt_200  <- as.matrix(read.csv("a-alt-200-Net1.2.csv"))
a_alt_400  <- as.matrix(read.csv("a-alt-400-Net1.2.csv"))
a_null_200 <- as.matrix(read.csv("a-null-200-Net1.2.csv"))
a_null_400 <- as.matrix(read.csv("a-null-400-Net1.2.csv"))

c_alt_200  <- as.matrix(read.csv("c-alt-200-Net1.2.csv"))
c_alt_400  <- as.matrix(read.csv("c-alt-400-Net1.2.csv"))
c_null_200 <- as.matrix(read.csv("c-null-200-Net1.2.csv"))
c_null_400 <- as.matrix(read.csv("c-null-400-Net1.2.csv"))

get_curve_leq <- function(output_table, c_list) {
  sapply(c_list, function(cc) {
    mean(output_table[, 1] <= cc * output_table[, 2])
  })
}

get_curve_gt <- function(output_table, c_list) {
  sapply(c_list, function(cc) {
    mean(output_table[, 1] > cc * output_table[, 2])
  })
}

plot_df <- rbind(
  data.frame(c = c_list, y = get_curve_leq(a_alt_200,  c_list),  group = "a-alt",  n = "n = 200"),
  data.frame(c = c_list, y = get_curve_leq(a_null_200, c_list),  group = "a-null", n = "n = 200"),
  data.frame(c = c_list, y = get_curve_leq(c_alt_200,  c_list),  group = "c-alt",  n = "n = 200"),
  data.frame(c = c_list, y = get_curve_gt(c_null_200,  c_list),  group = "c-null", n = "n = 200"),
  
  data.frame(c = c_list, y = get_curve_leq(a_alt_400,  c_list),  group = "a-alt",  n = "n = 400"),
  data.frame(c = c_list, y = get_curve_leq(a_null_400, c_list),  group = "a-null", n = "n = 400"),
  data.frame(c = c_list, y = get_curve_leq(c_alt_400,  c_list),  group = "c-alt",  n = "n = 400"),
  data.frame(c = c_list, y = get_curve_gt(c_null_400,  c_list),  group = "c-null", n = "n = 400")
)

plot_df$group <- factor(plot_df$group, levels = c("a-alt", "a-null", "c-alt", "c-null"))
plot_df$n <- factor(plot_df$n, levels = c("n = 200", "n = 400"))

ggplot(plot_df, aes(x = c, y = y, color = group, linetype = group, shape = group)) +
  geom_line(linewidth = 0.8, alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ n, nrow = 1) +
  scale_x_continuous(breaks = c_list) +
  scale_y_continuous(limits = c(-0.05, 0.1)) +
  scale_linetype_manual(
    values = c(
      "a-alt"  = "solid",
      "a-null" = "dashed",
      "c-alt"  = "dotdash",
      "c-null" = "twodash"
    ),
    labels = c(
      "a-alt"  = "(a)-alt",
      "a-null" = "(a)-null",
      "c-alt"  = "(c)-alt",
      "c-null" = "(c)-null"
    )
  ) +
  scale_shape_manual(
    values = c(
      "a-alt"  = 16,
      "a-null" = 17,
      "c-alt"  = 4,
      "c-null" = 3
    ),
    labels = c(
      "a-alt"  = "(a)-alt",
      "a-null" = "(a)-null",
      "c-alt"  = "(c)-alt",
      "c-null" = "(c)-null"
    )
  ) +
  scale_color_discrete(
    labels = c(
      "a-alt"  = "(a)-alt",
      "a-null" = "(a)-null",
      "c-alt"  = "(c)-alt",
      "c-null" = "(c)-null"
    )
  ) +
  labs(
    x = "Value of C",
    y = "Proportion of misdiagnosis",
    color = NULL,
    linetype = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    plot.margin = margin(5.5, 5.5, 2, 5.5),
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 11),
    axis.title.x = element_text(size = 11)
  )
