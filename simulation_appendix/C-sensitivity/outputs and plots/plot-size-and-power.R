library(ggplot2)

C_list <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6)

a_null_200 <- as.matrix(read.csv("a-null-200-Net1.2.csv"))
a_null_400 <- as.matrix(read.csv("a-null-400-Net1.2.csv"))
c_null_200 <- as.matrix(read.csv("c-null-200-Net1.2.csv"))
c_null_400 <- as.matrix(read.csv("c-null-400-Net1.2.csv"))
a_alt_200  <- as.matrix(read.csv("a-alt-200-Net1.2.csv"))
a_alt_400  <- as.matrix(read.csv("a-alt-400-Net1.2.csv"))
c_alt_200  <- as.matrix(read.csv("c-alt-200-Net1.2.csv"))
c_alt_400  <- as.matrix(read.csv("c-alt-400-Net1.2.csv"))

get_rej_rate <- function(output_table, C_list) {
  sapply(C_list, function(C) {
    vc <- ifelse(output_table[, 1] > C * output_table[, 2],
                 output_table[, 3],
                 output_table[, 4])
    sum(abs(vc) > qnorm(0.975, 0, 1)) / length(vc)
  })
}

plot_df <- rbind(
  data.frame(C = C_list, y = get_rej_rate(a_null_200, C_list), panel = "(a)-null", n = "n=200"),
  data.frame(C = C_list, y = get_rej_rate(a_null_400, C_list), panel = "(a)-null", n = "n=400"),
  data.frame(C = C_list, y = get_rej_rate(c_null_200, C_list), panel = "(c)-null", n = "n=200"),
  data.frame(C = C_list, y = get_rej_rate(c_null_400, C_list), panel = "(c)-null", n = "n=400"),
  data.frame(C = C_list, y = get_rej_rate(a_alt_200,  C_list), panel = "(a)-alt",  n = "n=200"),
  data.frame(C = C_list, y = get_rej_rate(a_alt_400,  C_list), panel = "(a)-alt",  n = "n=400"),
  data.frame(C = C_list, y = get_rej_rate(c_alt_200,  C_list), panel = "(c)-alt",  n = "n=200"),
  data.frame(C = C_list, y = get_rej_rate(c_alt_400,  C_list), panel = "(c)-alt",  n = "n=400")
)

plot_df$panel <- factor(plot_df$panel, levels = c("(a)-null", "(c)-null", "(a)-alt", "(c)-alt"))
plot_df$n <- factor(plot_df$n, levels = c("n=200", "n=400"))

lab_map <- c(
  "(a)-null" = "(a)-null\nType-I error",
  "(c)-null" = "(c)-null\nType-I error",
  "(a)-alt"  = "(a)-alt\nPower",
  "(c)-alt"  = "(c)-alt\nPower"
)

# invisible points to force panel-specific y-ranges
ylim_df <- rbind(
  data.frame(panel = "(a)-null", C = c(min(C_list), max(C_list)), y = c(0.04, 0.065)),
  data.frame(panel = "(c)-null", C = c(min(C_list), max(C_list)), y = c(0.04, 0.065)),
  data.frame(panel = "(a)-alt",  C = c(min(C_list), max(C_list)), y = c(0.90, 1.05)),
  data.frame(panel = "(c)-alt",  C = c(min(C_list), max(C_list)), y = c(0.90, 1.05))
)

ylim_df$panel <- factor(ylim_df$panel, levels = levels(plot_df$panel))

ggplot(plot_df, aes(x = C, y = y, color = n, linetype = n, shape = n)) +
  geom_blank(data = ylim_df, inherit.aes = FALSE, aes(x = C, y = y)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y", labeller = as_labeller(lab_map)) +
  scale_x_continuous(breaks = C_list) +
  scale_color_manual(values = c("n=200" = "skyblue", "n=400" = "steelblue")) +
  scale_linetype_manual(values = c("n=200" = "solid", "n=400" = "dashed")) +
  scale_shape_manual(values = c("n=200" = 16, "n=400" = 17)) +
  labs(
    x = "Value of C",
    y = NULL,
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
    axis.title.x = element_text(size = 11)
  )
