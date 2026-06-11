library(readr)
library(dplyr)

mat2 <- as.matrix(read_csv("ready_to_use_data_matrices/Business General Adjacency Matrix.csv",
                           col_names = FALSE))
N <- ncol(mat2)
colnames(mat2) <- NULL
mat2 <- unname(mat2)
Bs112_mat <- mat2
diag(Bs112_mat) <- NA
mu_e_empirical <- mean(na.omit(as.vector(Bs112_mat)))

# local reciprocity effect for i=1--n
local_2 <- rep(NA,N)
for (i in 1:N) {
  local_2[i] <- mean(na.omit((Bs112_mat[i,]-mu_e_empirical) * (Bs112_mat[,i]-mu_e_empirical)))
}

# local same-sender effect for i=1--n
local_3 <- rep(NA,N)
for (i in 1:N) {
  local_3[i] <- mean(apply(combn(na.omit(Bs112_mat[i,]-mu_e_empirical), 2), 2, prod))
}

# local same-receiver effect for i=1--n
local_4 <- rep(NA,N)
for (i in 1:N) {
  local_4[i] <- mean(apply(combn(na.omit(Bs112_mat[,i]-mu_e_empirical), 2), 2, prod))
}

# local sender-receiver effect for i=1--n
local_5 <- rep(NA,N)
for (i in 1:N) {
  local_5[i] <- (sum(na.omit(Bs112_mat[,i] - mu_e_empirical)) * sum(na.omit(Bs112_mat[i,] - mu_e_empirical)) -
                   sum(na.omit((Bs112_mat[i,]-mu_e_empirical) * (Bs112_mat[,i]-mu_e_empirical)))) / (N-1) / (N-2)
}

df <- data.frame(
  "Reciprocity"=local_2,
  "Same-sender"=local_3,
  "Same-receiver"=local_4,
  "Sender-receiver"=local_5,
  check.names=FALSE
)
library(tidyr)
long_df <- gather(df, key = "Local network effects", value = "value")
long_df <- subset(long_df, value >= -0.4 & value <= 0.6)

library(ggplot2)

# BS
ggplot(long_df, aes(x=value, fill=`Local network effects`, color = `Local network effects`)) +
  geom_density(alpha=0.5) +
  facet_grid(`Local network effects` ~ .) +
  xlim(-0.4, 0.6) + 
  scale_fill_manual(values=c("orangered1", "gold1", "cyan2","orchid2")) +
  scale_color_manual(values=c("orangered1", "orange1", "cyan3","orchid2")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_light() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text = element_text(size=13),
        legend.position = c(0.2, 0.66),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.spacing.y = unit(1.8, 'cm'),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 13), breaks = c(0,5,10)) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x = "Local network effects", y = "Density") 



#'@_scatter_plot_of_eta_2&3

Dataset_2_Business_vertexlist <- read_delim("ready_to_use_data_matrices/Dataset 2. Business_vertexlist.txt", 
                                            "\t", escape_double = FALSE, trim_ws = TRUE)

schools <- Dataset_2_Business_vertexlist$institution
schools2 <- rep(NA,112)
schools2[c(1,2,3,4,6,8,9,10,12,15,16,19)] <- schools[c(1,2,3,4,6,8,9,10,12,15,16,19)]
schools2[c(1,2,3,4,6,8,9,10,12,15,16,19)] <- c("Stanford","MIT","Harvard","Berkeley","Chicago","Northwestern","Michigan","Columbia",
                                               "UPenn","UCLA","NYU","Duke")
schoolindx <- Dataset_2_Business_vertexlist$'# u'

df <- data.frame(
  "Reciprocity"=local_2,
  "Same-sender"=local_3,
  "name"=schools2[1:112],
  "nameindx" = schoolindx[1:112],
  check.names=FALSE
)

library(ggplot2)


vjust_1 <- rep(-1.2,112)
vjust_1[3] <- 2
vjust_1[4] <- -1
vjust_1[6] <- -1
vjust_1[8] <- 0.5
vjust_1[10] <- 0
vjust_1[12] <- -1.6
vjust_1[15] <- 1
vjust_1[16] <- 0.5
vjust_1[19] <- 2
hjust_1 <- rep(0.5,112)
hjust_1[4] <- 0.82
hjust_1[8] <- -0.1
hjust_1[10] <- -0.1
hjust_1[15] <- -0.1
hjust_1[16] <- 1.28

color_temp <- rep("darkgray",112)
color_temp[c(1,2,3,4,6,8,9,10,12,15,16,19)] <- "red"
size_temp <- rep(0.5,112)
size_temp[c(1,2,3,4,6,8,9,10,12,15,16,19)] <- 1

ggplot(df, aes(Reciprocity, `Same-sender`)) +
  geom_point(size = size_temp, color = color_temp) +
  xlab("Reciprocity") +
  ylab("Same-sender") +
  xlim(-0.5, 10) +
  ylim(-0.4, 9) +
  geom_text(aes(label = name), size = 3.8,
            vjust = vjust_1, hjust = hjust_1,
            na.rm = TRUE) +
  theme_light() +
  theme(axis.text.y = element_text(size = 11),
        axis.title.y = element_text(vjust = 0.8, size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(vjust = 0.8, size = 12))

