library(readr)
library(dplyr)

mat2 <- as.matrix(read_csv("ready_to_use_data_matrices/Business General Adjacency Matrix.csv",col_names = FALSE))
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

#'@_radar_plot_of_10+_schools

df <- data.frame(
  "Reciprocity"=local_2,
  "Same- \n sender"=local_3,
  "Same-receiver"=local_4,
  "Sender- \n receiver"=local_5,
  check.names=FALSE
)

Dataset_2_Business_vertexlist <- read_delim("ready_to_use_data_matrices/Dataset 2. Business_vertexlist.txt", 
                                            "\t", escape_double = FALSE, trim_ws = TRUE)
ranks_temp <- as.numeric(Dataset_2_Business_vertexlist$USN2012[c(1,2,3,4,6,8,9,10,12,15,16,19)])


df_radar <- df[c(c(1,2,3,4,6,8,9,10,12,15,16,19)),]
df_radar$rank <- ranks_temp

sc_names <- c("Stanford","MIT","Harvard","Berkeley","Chicago","Northwestern","Michigan","Columbia",
              "UPenn","UCLA","NYU","Duke")
rownames(df_radar) <- c("Stanford","MIT","Harvard","Berkeley","Chicago","Northwestern","Michigan","Columbia",
                        "UPenn","UCLA","NYU","Duke")
# add rank in name
rownames(df_radar) <- c("Stanford \nRank:1","MIT \nRank:3","Harvard \nRank:2","Berkeley \nRank:7","Chicago \nRank:5","Northwestern \nRank:5","Michigan \nRank:14","Columbia \nRank:9",
                        "UPenn \nRank:3","UCLA \nRank:14","NYU \nRank:10","Duke \nRank:12")

df_radar <- df_radar[order(df_radar$rank), ]
df_radar <- df_radar[,-5]

colnames(df_radar) <- c("Reciprocity","Same- \n sender","Same-receiver","Sender- \n receiver")

max_min <- data.frame(
  "Reciprocity"=c(9, -0.5),
  "Same- \n sender"=c(9, -0.5),
  "Same-receiver"=c(9, -0.5),
  "Sender- \n receiver"=c(9, -0.5),
  check.names=FALSE
)
rownames(max_min) <- c("Max", "Min")

df_radar <- rbind(max_min, df_radar)

library(fmsb)

create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.9,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1, 
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 1, plty = 1,
    # Customize the grid
    cglcol = "grey68", cglty = 1, cglwd = 1, 
    # Customize the axis
    axislabcol = "grey6", seg = 3, 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, calcex =1.068, palcex=1.068, cex.main=1.1
  )
}



op <- par(mar = c(0.5,0.5,3.5,0.5))
par(mfrow = c(3,4))

# Create the radar chart
for(i in 1:12){
  create_beautiful_radarchart(
    data = df_radar[c(1, 2, i+2), ], caxislabels = c(-0.5,2.5,5.5,8.5),
    color = "skyblue", title = rownames(df_radar)[i+2]
  )
}
par(op)


















