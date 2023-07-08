# 需要先把普通spatial.csv数据格式调整为index: spot_id;  column1: x;  column2:y
library(ggplot2)
source("D:\\Users\\lqlu\\work\\Codes\\STABox\\DataClipping\\select_cells.R")
# pd <- read.csv("D:\\Users\\lqlu\\work\\Codes\\Puck_200127_15_bead_locations.csv", sep = ',', row.names = 1)
pd <- read.csv("D:\\Users\\lqlu\\work\\Codes\\locations_clean.csv", sep = ',',row.names = 1)
colnames(pd) = c("x","y")
out_path <- "D:\\Users\\lqlu\\work\\Codes\\STABox\\DataClipping"
pd$ce11 = rownames(pd)
runGadget(ui, server)


#source("D:\\Users\\lqlu\\work\\Codes\\STABox\\DataClipping\\select_cells.R")
#new_pd <- read.csv("D:\\Users\\lqlu\\work\\Codes\\locations_clean.csv", sep = ',',row.names = 1)
#out_path <- "D:\\Users\\lqlu\\work\\Codes\\STABox\\DataClipping"
#runGadget(ui, server)

# 需要先把普通spatial.csv数据格式调整为index: spot_id;  column1: x;  column2:y
library(ggplot2)
source("D:\\Users\\lqlu\\work\\Codes\\STABox\\DataClipping\\select_cells.R")
# pd <- read.csv("D:\\Users\\lqlu\\work\\Codes\\Puck_200127_15_bead_locations.csv", sep = ',', row.names = 1)
pd <- read.csv("D:\\Users\\lqlu\\download\\data_coor.csv", sep = ',',row.names = 1)
colnames(pd) = c("x","y")
out_path <- "D:\\Users\\lqlu\\download"
pd$ce11 = rownames(pd)
runGadget(ui, server)