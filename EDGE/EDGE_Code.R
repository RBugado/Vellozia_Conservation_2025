#### Calculating GE2 and EDGE2 for Vellozia species ####
# July 2024
library(ape)

wd <- "~/VelloziaProject2024/R Things Vellozia 2024/EDGE_R"
setwd(wd)

### GE2 ###
# generate 1,000 GE2 values distributed across all Red List categories 
# and assign a distribution of GE2 values to each RL category to capture uncertainty

GE.2.calc <- function(pext){
  require(geiger)
  treesim <- sim.bdtree(n=10000)
  iucn <- sample(1:5, size=length(treesim$tip.label), replace=TRUE)
  data <- data.frame(species=treesim$tip.label, pext=pext[iucn])
  data <- data[order(data$pext),]
  data$rank <- seq_len(nrow(data))
  rank <- c(0, with(data, tapply(rank, pext, median)))
  pext <- c(0, pext)
  rank.sq <- rank^2; rank.cub <- rank^3; rank.qu <- rank^4; rank.quu <- rank^5
  model <- lm(pext ~ rank + rank.sq + rank.cub + rank.qu)
  data$rank.sq <- data$rank^2; data$rank.cub <- data$rank^3; data$rank.qu <- data$rank^4; data$rank.quu <- data$rank^5
  data$rank.pext <- predict(model, data)
  data$rank.pext[data$rank.pext <= 0] <- 0.0001
  data$rank.pext[data$rank.pext >= 1] <- 0.9999
  pext.LC <- data.frame(RL.cat = "LC", pext =data$rank.pext[data$pext == pext[2]])
  pext.NT <- data.frame(RL.cat = "NT", pext =data$rank.pext[data$pext == pext[3]])
  pext.VU <- data.frame(RL.cat = "VU", pext =data$rank.pext[data$pext == pext[4]])
  pext.EN <- data.frame(RL.cat = "EN", pext =data$rank.pext[data$pext == pext[5]])
  pext.CR <- data.frame(RL.cat = "CR", pext =data$rank.pext[data$pext == pext[6]])
  return(rbind(pext.CR,pext.EN, pext.VU, pext.NT, pext.LC))
}

EDGE.2.calc <- function(tree, pext){
  require(phylobase)
  require(data.table)
  if(!class(tree) == "phylo"){
    tree <- as(tree, "phylo")
  }
  tree_dat <- data.frame(Species = as.character(unique(tree$tip.label)),
                         TBL = tree$edge.length[sapply(c(1:length(tree$tip.label)),
                                                       function(x,y) which(y==x),y=tree$edge[,2])], 
                         pext = NA, ED = NA, EDGE = NA)
  ePD.dat <- data.frame(PD = sum(tree$edge.length),ePDloss = NA)
  tree <- as(tree, "phylo4")
  names(pext) <- c("species","pext")
  for(i in 1:length(tree_dat$Species)){
    tree_dat$pext[i] <- pext$pext[pext$species == tree_dat$Species[i]]
  }
  nodes <- descendants(tree, rootNode(tree), "all")
  for(i in 1:length(nodes)){
    tips <- descendants(tree, nodes[i], "tips")
    tips <- names(tips)
    tipscores <- which(pext$species %in% tips)
    tree@edge.length[which(tree@edge[,2] == nodes[i])] <- edgeLength(tree, nodes[i])*prod(as.numeric(pext$pext[tipscores]))
  }
  for(i in 1:length(tree_dat$Species)){
    tree_dat$EDGE[i] <- sum(tree@edge.length[which(tree@edge[,2] %in% ancestors(tree, 
                                                                                which(tipLabels(tree) == tree_dat$Species[i]), "ALL"))], na.rm=T)
    tree_dat$ED[i] <- tree_dat$EDGE[i] / as.numeric(tree_dat$pext[i])
  }
  tree <- as(tree, "phylo")
  ePD.dat$ePDloss <- sum(tree$edge.length)
  edge.res <- list(tree_dat,tree,ePD.dat)
  return(edge.res)
}


# calculating ge2
pext <- rev(c(0.97, 0.97/2, 0.97/4,0.97/8,0.97/16))
ge2 <- GE.2.calc(pext)


### EDGE2 ###
# provide dated phylogenetic tree and dataframe with two columns: 
# the first comprising species names (column name: species), the with the letter IUCN category (column name: iucn)
# function returns three objects: 
# 1. dataframe with terminal branch length, GE2, ED2 and EDGE2 scores for each species
# 2. expected PD loss tree
# 3. PD and expected PD loss in MY for the clade

tree_c <- read.tree("Vellozia_EDGE_Phylogeny_BEAST_74.tre")
data <- read.csv("Vellozia_Species_vs_Extinction_Risk_letters_74.csv")
data_to_edge <- data
results <- list()
for(i in 1:1000) {
  for(species_index in 1:nrow(data)){
    one_ge_value <- data[,2][species_index]
    if(one_ge_value=="CR") {
      one_sample <- round(sample(ge2$pext[which(ge2$RL.cat=="CR")], 1),3)
    }
    if(one_ge_value=="EN") {
      one_sample <- round(sample(ge2$pext[which(ge2$RL.cat=="EN")], 1),3)
    }
    if(one_ge_value=="VU") {
      one_sample <- round(sample(ge2$pext[which(ge2$RL.cat=="VU")], 1),3)
    }
    if(one_ge_value=="LC") {
      one_sample <- round(sample(ge2$pext[which(ge2$RL.cat=="LC")], 1),3)
    }
    data_to_edge$iucn[species_index]<-one_sample
  }
  # colnames(data_to_edge)[2] <- "GE2"
  one_calc <- EDGE.2.calc(tree=tree_c, pext=data_to_edge)
  results[[i]] <- one_calc[[1]]
  cat(i, "\r")
}


all_results <- do.call(rbind, results)

#### Plotting the results ####

# Creating range bars
# Miminum and maximum EDGE value in a dataframe

sp <- rev(tree_c$tip.label)
plot_data <- list()
for(i in 1:length(unique(all_results$Species))){
  sp_now <- sp[i]
  all_values <- all_results[all_results$Species == sp[i],]
  Max <- max(all_values$EDGE)
  Min <- min(all_values$EDGE)
  plot_data[[i]] <- data.frame(sp_now, Min, Max)
}
final_plot_data <- do.call(rbind, plot_data)

x <- 1:74  # number of species in the tree
High <-  final_plot_data$Max
Low  <-  final_plot_data$Min

## Blank plot
plot(1, type="n", xlab="species", ylab="EDGE2 range", xlim=range(x), 
     ylim=c(min(Low), max(High)))

## Add rectangles
# This plot is the range of EDGE values for each species
rect(x - 0.2, Low, x + 0.2, High, col="red")  


# This plot is a barplot of the mean EDGE value for each species
sp <- rev(tree_c$tip.label)
plot_data <- list()
for(i in 1:length(unique(all_results$Species))){
  sp_now <- sp[i]
  all_values <- all_results[all_results$Species == sp[i],]
  Mean <- mean(all_values$EDGE)
  plot_data[[i]] <- data.frame(sp_now, Mean)
}
final_plot_data2 <- do.call(rbind, plot_data)
barplot(final_plot_data2$Mean)

# This is a horizontal bar chart with the mean EDGE value for each species
# It is in vaguely in order of species for the descending phylogenetic tree, but should be checked to match the tree order
barplot(rev(final_plot_data2$Mean), horiz=TRUE, width = 1, space=0, xlab = "Mean EDGE2 Value",
        ylim=c(1,length(tree_c$tip.label))-0.75, xlim = c(0, 14), names="", col="azure3")

write.csv(final_plot_data2, "final_plot_data2.csv") ## this was to reorder the species to match the tree
final_plot_data2 <- read.csv("final_plot_data2.csv")

## Matches colors with the IUCN category for each species
coltips <- read.csv("Vellozia_Species_vs_Extinction_Risk_colors.csv",
                      header = TRUE, sep = ",", row.names = 1)
coltips <- coltips[tree_c$tip.label,] # For tip labels

colplot <- read.csv("Vellozia_Species_vs_Extinction_Risk_colors.csv",
                    header = TRUE, sep = ",", row.names = 1)
colplot <- colplot[final_plot_data2$sp_now,] # For graph bars

## pdf version of the tree and the mean EDGE bar graph next to each other

pdf("Vellozia_EDGE_Phylogeny_with_EDGE_plot.pdf")
par(mai=c(0.5,0.1,0.1,0.2)) #make margins smaller
par(fig=c(0,0.8,0,1)) #define area
  plot(ladderize(tree_c,TRUE), cex = 0.6, label.offset = 1)
  axisPhylo()
  tiplabels(pch = 21, cex = 0.8, 
            bg=as.character(coltips), adj = 1)
par(fig=c(0.65,1,0,1), new=TRUE)
  barplot(rev(final_plot_data2$Mean), horiz=TRUE, width = 1, space=0, 
          xlab = "Mean EDGE2 Value", col = rev(colplot),
          ylim=c(1,length(tree_c$tip.label))-0.75, xlim = c(0, 14), names="")
  dev.off()

## Barplot with median instead of mean - median is the correct way to display this data
sp <- rev(tree_c$tip.label)
plot_data <- list()
for(i in 1:length(unique(all_results$Species))){
    sp_now <- sp[i]
    all_values <- all_results[all_results$Species == sp[i],]
    Median <- median(all_values$EDGE)
    plot_data[[i]] <- data.frame(sp_now, Median)
  }
final_plot_data3 <- do.call(rbind, plot_data)
barplot(final_plot_data3$Median)

barplot(rev(final_plot_data3$Median), horiz=TRUE, width = 1, space=0, 
        xlab = "Median EDGE2 Value", col = rev(colplot),
        ylim=c(1,length(tree_c$tip.label))-0.75, xlim = c(0, 14), names="")

write.csv(final_plot_data3, "final_plot_data3.csv") # Again, this is reordering the data to match the tree
final_plot_data3 <- read.csv("final_plot_data3.csv")

## PDF of the tree and median bar plot together
pdf("Vellozia_EDGE_Phylogeny_with_EDGE_median.pdf")
par(mai=c(0.6,0.1,0.1,0.2)) #make margins smaller
par(fig=c(0,0.78,0,1)) #define area
  plot(ladderize(tree_c,TRUE), cex = 0.6, label.offset = 1)
  axisPhylo(cex.axis = 0.7)
  tiplabels(pch = 21, cex = 0.8, 
          bg=as.character(coltips), adj = 1)
par(fig=c(0.63,1,0,1), new=TRUE)
  barplot(rev(final_plot_data2$Mean), horiz=TRUE, width = 1, 
          space=0, col = rev(colplot), cex.axis = 0.7, 
          ylim=c(1,length(tree_c$tip.label))-0.75, xlim = c(0, 14), names="")
  dev.off()

  
##Boxplots  

myplot <- boxplot(EDGE ~ Species, 
                  data = all_results)

summary(final_plot_data3)
myplot2 <- boxplot(final_plot_data3$Median, horizontal = TRUE)

                                                
## Create summary of the EDGE table, 
library(dplyr)

# mean by species for each column
numeric_all_results <- transform(all_results, pext = as.numeric(pext))

edge2_table_mean <- numeric_all_results %>% 
  group_by(Species) %>% 
  summarize(TBL = mean(TBL),
            pext = mean(pext),
            ED = mean(ED),
            EDGE = mean(EDGE))

write.csv(edge2_table_mean, "edge2_table_mean_Vellozia.csv")


# median by species for each column
numeric_all_results <- transform(all_results, pext = as.numeric(pext))

edge2_table_median <- numeric_all_results %>% 
  group_by(Species) %>% 
  summarize(TBL = median(TBL),
            pext = median(pext),
            ED = median(ED),
            EDGE = median(EDGE))

write.csv(edge2_table_median, "edge2_table_median_Vellozia.csv")

