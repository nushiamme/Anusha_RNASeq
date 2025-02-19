## Code to plot phylogeny of KEGG/Biomart species with available annotations
## Code by: Anusha Shankar, github/nushiamme; 
# contact: nushiamme<at>gmail<dot>com


### Contents
## Setup, read files in, format data
## Figures

library(MCMCglmm)
library(nlme)
library(ape)
library(geiger) # for treedata() function
library(caper)
library(phytools)
library(RColorBrewer)
library(ggplot2)

#### Setup ####
setwd("/Users/ashankar/Library/CloudStorage/OneDrive-CornellUniversity/Published_papers_AShankar/JAnimEcol_Allometry")
fmr_data <- read.csv("DLW_TableS1_final.csv") # sub this with simple list of species names you want

## Read in McGuire et al. 2014 hummingbird phylogeny; contact nushiamme<at>gmail<dot>com if you cannot get access.
all_tree <-read.tree("hum294.tre") ## sub this with hackett tree
#tre_ou_edited <- read.tree("OU_hummer_tree_FMR_edit.txt")

## General plotting functions
## Generic theme
my_theme <- theme_classic(base_size = 35) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) 
# 
# ## To add linear regression equation to plot
# lm_eqn <- function(y, x){
#   m <- lm(y ~ x);
#   eq <- substitute(italic(y) == 
#                      a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                    list(a = format(coef(m)[1], digits = 2), 
#                         b = format(coef(m)[2], digits = 2), 
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));                 
# }
# 
# ## Aggregating dataset by Species and site, to get species means of mass and daily energy expenditure (DEE)
# fmr_data$Lit_data2 <- fmr_data$Lit_data
# levels(fmr_data$Lit_data2)[match("Data_2016",levels(fmr_data$Lit_data2))] <- "Data"
# levels(fmr_data$Lit_data2)
# dlw_mean <- data.frame()
# dlw_mean <- aggregate(fmr_data$kJ_day, by=list(fmr_data$Species, fmr_data$Big_site, fmr_data$Lit_data2), FUN="mean", na.omit=T)
# dlw_mass <- aggregate(fmr_data$Mass_g, by=list(fmr_data$Species, fmr_data$Big_site, fmr_data$Lit_data2), FUN="mean", na.omit=T)
# dlw_mean <- merge(dlw_mean, dlw_mass, by = c("Group.1", "Group.2", "Group.3"))
# names(dlw_mean) <- c("Species", "Region", "Lit_data2", "kJ_day", "Mass_g")

## Trimming tree to DLW dataset
## Manually replacing because it's a manageable number
tree_dlw$tip.label[1]<-"FLME"
tree_dlw$tip.label[15]<-"PHYA"
tree_dlw$tip.label[83]<-"URBE"
tree_dlw$tip.label[92]<-"HEIM"
tree_dlw$tip.label[93]<-"HERU"
tree_dlw$tip.label[95]<-"HEJA"
tree_dlw$tip.label[128]<-"AGCO"
#tree_dlw$tip.label[154]<-"PAGI" ## Doing this one separately later
tree_dlw$tip.label[156]<-"EUFU"
tree_dlw$tip.label[163]<-"LACL"
tree_dlw$tip.label[185]<-"ARAL"
tree_dlw$tip.label[188]<-"CAAN"
tree_dlw$tip.label[219]<-"CYLA"
tree_dlw$tip.label[230]<-"CHUR"
tree_dlw$tip.label[234]<-"THCO"
tree_dlw$tip.label[235]<-"THFA"
tree_dlw$tip.label[269]<-"AMTZ"


## Tree without the Giant hummingbird
tree_no_Pgigas <- tree_dlw

tree_dlw$tip.label[154]<-"PAGI"

tips<-data.frame(levels(fmr_data$Species))
colnames(tips) <- "tips"
rownames(tips)<-tips$tips

#match tree to data, prune tree, species names should be in rownnames of "data" 
tre1<-treedata(tree_dlw, tips)$phy
#To check that the relationships between species in the trimmed tree look right
plot(tre1, cex=1.5, edge.width = 3) 
## Matching tree without P. gigas and trimming
tips2<-data.frame(levels(droplevels(fmr_data$Species[fmr_data$Species != "PAGI"])))
colnames(tips2) <- "tips"
rownames(tips2)<-tips2$tips
