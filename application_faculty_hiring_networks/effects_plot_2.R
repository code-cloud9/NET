library(readr)
library(dplyr)

# History
mat2 <- as.matrix(read_csv("ready_to_use_data_matrices/History General Adjacency Matrix.csv",col_names = FALSE))
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
long_df <- subset(long_df, value >= -0.06 & value <= 0.1)

library(ggplot2)
ggplot(long_df, aes(x=value, fill=`Local network effects`, color = `Local network effects`)) +
  geom_density(alpha=0.5) +
  facet_grid(`Local network effects` ~ .) +
  xlim(-0.06, 0.1) +
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
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x = "Local network effects", y = "Density", title = "History") 


# Computer Science
mat2 <- as.matrix(read_csv("ready_to_use_data_matrices/CS General Adjacency Matrix.csv",col_names = FALSE))
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
long_df <- subset(long_df, value >= -0.015 & value <= 0.025)

ggplot(long_df, aes(x=value, fill=`Local network effects`, color = `Local network effects`)) +
  geom_density(alpha=0.5) +
  facet_grid(`Local network effects` ~ .) +
  xlim(-0.015, 0.025) +  
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
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(limits = c(0, 260), breaks = c(0,80,160,240)) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x = "Local network effects", y = "Density", title = "Computer science") 






