#LBF Community Shift last 20 years data analysis (not PCA, not NMDS)

#TO RUN BEFORE ANYTHING ELSE - all packages needed for all analysis
library(dplyr)
library(reshape2)
library(vegan)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(broom)
library(forcats)
library(alluvial)
library(pheatmap)
library(indicspecies)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(ggpmisc)
library(naniar)
library(visNetwork)
library(cooccur)
library(igraph)
library(ggpubr)

"%not%" <- Negate("%in%") #to create a "not in" sign, the opposite of %in%

#function to plot co-occurrence matrix
plotforcooccur <- function (x, ...) 
{
  allargs <- match.call(expand.dots = TRUE)
  plotrand <- allargs$plotrand
  plotrand <- ifelse(test = is.null(plotrand), yes = FALSE, 
                     no = plotrand)
  randsummary <- allargs$randsummary
  randsummary <- ifelse(test = is.null(randsummary), yes = FALSE, 
                        no = randsummary)
  dim <- x$species
  comat_pos <- comat_neg <- matrix(nrow = dim, ncol = dim)
  co_tab <- x$result
  for (i in 1:nrow(co_tab)) {
    comat_pos[co_tab[i, "sp1"], co_tab[i, "sp2"]] <- co_tab[i, 
                                                            "p_gt"]
    comat_pos[co_tab[i, "sp2"], co_tab[i, "sp1"]] <- co_tab[i, 
                                                            "p_gt"]
    row.names(comat_pos[co_tab[i, "sp2"], co_tab[i, 
                                                 "sp1"]])
  }
  for (i in 1:nrow(co_tab)) {
    comat_neg[co_tab[i, "sp1"], co_tab[i, "sp2"]] <- co_tab[i, 
                                                            "p_lt"]
    comat_neg[co_tab[i, "sp2"], co_tab[i, "sp1"]] <- co_tab[i, 
                                                            "p_lt"]
  }
  comat <- ifelse(comat_pos >= 0.05, 0, 1) + ifelse(comat_neg >= 
                                                      0.05, 0, -1)
  colnames(comat) <- 1:dim
  row.names(comat) <- 1:dim
  if ("spp_key" %in% names(x)) {
    sp1_name <- merge(x = data.frame(order = 1:length(colnames(comat)), 
                                     sp1 = colnames(comat)), y = x$spp_key, by.x = "sp1", 
                      by.y = "num", all.x = T)
    sp2_name <- merge(x = data.frame(order = 1:length(row.names(comat)), 
                                     sp2 = row.names(comat)), y = x$spp_key, by.x = "sp2", 
                      by.y = "num", all.x = T)
    colnames(comat) <- sp1_name[with(sp1_name, order(order)), 
                                "spp"]
    row.names(comat) <- sp2_name[with(sp2_name, order(order)), 
                                 "spp"]
  }
  comat[is.na(comat)] <- 0
  origN <- nrow(comat)
  
  #if (plotrand == FALSE) {
  #  ind <- apply(comat, 1, function(x) all(x == 0))
  #  comat <- comat[!ind, ]
  #  ind <- apply(comat, 2, function(x) all(x == 0))
  #  comat <- comat[, !ind]
  #}
  postN <- nrow(comat)
  #ind <- apply(comat, 1, function(x) all(x == 0))
  #comat <- comat[names(sort(ind)), ]
  #ind <- apply(comat, 2, function(x) all(x == 0))
  #comat <- comat[, names(sort(ind))]
  data.m = melt(comat)
  colnames(data.m) <- c("X1", "X2", "value")
  data.m$X1 <- as.character(data.m$X1)
  data.m$X2 <- as.character(data.m$X2)
  meas <- as.character(unique(data.m$X2))
  dfids <- subset(data.m, X1 == X2)
  X1 <- data.m$X1
  X2 <- data.m$X2
  df.lower = subset(data.m[lower.tri(comat), ], X1 != X2)
  if (randsummary == FALSE) {
  }
  else {
    dim <- nrow(comat)
    ext.dim <- round(dim * 0.2, digits = 0)
    if (ext.dim < 0) {
      ext.dim <- 1
    }
    placehold <- paste("ext_", rep(c(1:ext.dim), each = dim), 
                       sep = "")
    randcol.df <- data.frame(X1 = placehold, X2 = rep(meas, 
                                                      times = ext.dim), value = rep(x = c(-2), times = dim * 
                                                                                      ext.dim))
    df.lower <- rbind(df.lower, randcol.df)
    meas <- c(meas, unique(placehold))
  }
  X1 <- df.lower$X1
  X2 <- df.lower$X2
  value <- df.lower$value
  
  p <- ggplot(df.lower, aes(X1, X2)) + geom_tile(aes(fill = factor(value, 
                                                                   levels = c(-1, 0, 1))), colour = "white")
  p <- p + scale_fill_manual(values = c("navy", 
                                        "gray90", "indianred1"), name = "", 
                             labels = c("negative", "random", "positive"), 
                             drop = FALSE) + theme(axis.text.x = element_blank(), 
                                                   axis.text.y = element_blank(), axis.ticks = element_blank(), 
                                                   plot.title = element_text(vjust = -4, size = 20, 
                                                                             face = "bold"), panel.background = element_rect(fill = "white", 
                                                                                                                             colour = "white"), panel.grid.major = element_blank(), 
                                                   panel.grid.minor = element_blank(), legend.position = c(0.9, 
                                                                                                           0.5), legend.text = element_text(size = 18)) + 
    ggtitle("Species Co-occurrence Matrix") + xlab("") + 
    ylab("") + scale_x_discrete(limits = meas, 
                                expand = c(0.3, 0), drop = FALSE) + scale_y_discrete(limits = meas, 
                                                                                     expand = c(0.3, 0), drop = FALSE)
  p <- p + geom_text(data = dfids, aes(label = X1), hjust = 1, 
                     vjust = 0, angle = -22.5)
  
  p
}



#---------To write new datasets: create depth groups, merge Elphidium, isolate 12 islands, remove rare species and add relative abundance---------------------------
df <- read.csv("~/Downloads/Dataset/20210217Living_LBF_1997_2010_2012_2013_2018-corrected.csv", sep = ";")

df$Sample_Year <- as.character(as.numeric(df$Sample_Year))
df$Zone <- as.character(as.numeric(df$Zone))
df$Species_Code <- as.character(df$Species_Code)
df$Species <- as.character(df$Species)

#create groups for reef slope depth: shallow slope, mid-slope, deep-slope based on relative depth normalized for whole Spermonde
df$Depth_slope <- c(1:nrow(df))
df <- df %>% mutate(Depth_slope=cut(Relative_Depth, breaks=c(0.0, 0.25, 0.75, Inf), labels=c("shallow","mid","deep")))

#Rename Calcarina sp2 for Calcarina hispida
df$Species_Code <- gsub("Calcarina_hispida", "Calcarina_hispida", df$Species_Code)
df$Species<- ifelse(df$Species_Code == "Calcarina_hispida", gsub("sp2", "hispida", df$Species), df$Species)

df$Species_Code <- as.factor(as.character(df$Species_Code))


#Selecting specific data from main dataframe using subset function
df_west <- subset(df, Island_Side == "West")
df <- df_west

#12 islands of interest for the study
df <- subset(df, subset = Sampling_Site %in% c("LL1","BN1", "LD2", "KG2", "PO2", "SA3", "PJ3", "LO3", "KU4", "LU4", "BA4", "LA5"))

#merge all Elphidium spp together
df_noelphi <- subset(df, subset = Genus %not% "Elphidium")
df_elphi <- subset(df, subset = Genus %in% "Elphidium")
df_elphi$Species <- "spp"
df_elphi$Species_Code <- "Elphidium_spp"
df_elphi <- df_elphi %>% group_by(Sample_Code,Station,Island,Sampling_Site,Island_Side,Genus,Species,Species_Code,Sample_Year, 
                                  Zone,Habitat,Inhabitated,Distance_Makassar_km,Depth_m,Relative_Depth,Substrate, Depth_slope) %>% summarise(N=sum(N))
df_elphi <- df_elphi[, c(1:13, 18, 14:17)]
df <- rbind.data.frame(df_noelphi[,1:18], df_elphi[,1:18])

#merge all duplicate species by summing the counts for each sample
df <- df %>% group_by(Sample_Code,Station,Island,Sampling_Site,Island_Side,Genus,Species,Species_Code,Sample_Year, 
                      Zone,Habitat,Inhabitated,Distance_Makassar_km,Depth_m,Relative_Depth,Substrate, Depth_slope) %>% summarise(N=sum(N))
df <- df[, c(1:13, 18, 14:17)]


#Create intermediate dataset - all species, with 12 islands, depth group and merged Elphidium spp.
write.csv(df, "~/Downloads/Dataset/20210217_LBFdataset_12islands_withElphidiummerged.csv")

#rare spieces removed from the dataset with < 0.1% mean relative abundance
#write new dataset with only all 17 most abundant species (without rare species) and with both true (N) and relative (RE) abundance
df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_withElphidiummerged.csv")

#To give the same format to this dataframe than the normalized one, but only dcast and melt - no change in data
df_N <- df

#Find in which percentage of the samples these species are found at least 1 specimen
alveo <- subset(df_N, Species_Code == "Alveolinella_quoyii")
alveo$N[alveo$N > 0] <- 1
ratio_alveo <- sum(alveo$N)/445*100
dent <- subset(df_N, Species_Code == "Dendritina_ambigua")
dent$N[dent$N > 0] <- 1
ratio_dent <- sum(dent$N)/445*100
laev <- subset(df_N, Species_Code == "Laevipeneroplis_malayensis")
laev$N[laev$N > 0] <- 1
ratio_laev <- sum(laev$N)/445*100
numm <- subset(df_N, Species_Code == "Nummulites_venosus")
numm$N[numm$N > 0] <- 1
ratio_numm <- sum(numm$N)/445*100
par1 <- subset(df_N, Species_Code == "Parasorites_sp1")
par1$N[par1$N > 0] <- 1
ratio_par1 <- sum(par1$N)/445*100
par2 <- subset(df_N, Species_Code == "Parasorites_sp2")
par2$N[par2$N > 0] <- 1
ratio_par2 <- sum(par2$N)/445*100
pene <- subset(df_N, Species_Code == "Peneroplis_pertusus")
pene$N[pene$N > 0] <- 1
ratio_pene <- sum(pene$N)/445*100
penesp <- subset(df_N, Species_Code == "Peneroplis_sp")
penesp$N[penesp$N > 0] <- 1
ratio_penesp <- sum(penesp$N)/445*100
penesp2 <- subset(df_N, Species_Code == "Peneroplis_sp2")
penesp2$N[penesp2$N > 0] <- 1
ratio_penesp2 <- sum(penesp2$N)/445*100

sphae <- subset(df_N, Species_Code == "Sphaerogypsina_globulus")
sphae$N[sphae$N > 0] <- 1
ratio_sphae <- sum(sphae$N)/445*100
spl <- subset(df_N, Species_Code == "Amphisorus_SpL")
spl$N[spl$N > 0] <- 1
ratio_spl <- sum(spl$N)/445*100
sps <- subset(df_N, Species_Code == "Amphisorus_SpS")
sps$N[sps$N > 0] <- 1
ratio_sps <- sum(sps$N)/445*100
radia <- subset(df_N, Species_Code == "Amphistegina_radiata")
radia$N[radia$N > 0] <- 1
ratio_radia <- sum(radia$N)/445*100
lobi <- subset(df_N, Species_Code == "Amphistegina_lobifera")
lobi$N[lobi$N > 0] <- 1
ratio_lobi <- sum(lobi$N)/445*100
speng <- subset(df_N, Species_Code == "Calcarina_spengleri")
speng$N[speng$N > 0] <- 1
ratio_speng <- sum(speng$N)/445*100
csp2 <- subset(df_N, Species_Code == "Calcarina_hispida")
csp2$N[csp2$N > 0] <- 1
ratio_csp2 <- sum(csp2$N)/445*100
mayo <- subset(df_N, Species_Code == "Calcarina_mayori")
mayo$N[mayo$N > 0] <- 1
ratio_mayo <- sum(mayo$N)/445*100
less <- subset(df_N, Species_Code == "Amphistegina_lessonii")
less$N[less$N > 0] <- 1
ratio_less <- sum(less$N)/445*100
elph <- subset(df_N, Species_Code == "Elphidium_spp")
elph$N[elph$N > 0] <- 1
ratio_elph <- sum(elph$N)/445*100
bacu <- subset(df_N, Species_Code == "Baculogypsinoides_spinosus")
bacu$N[bacu$N > 0] <- 1
ratio_bacu <- sum(bacu$N)/445*100
gaim <- subset(df_N, Species_Code == "Neorotalia_gaimardi")
gaim$N[gaim$N > 0] <- 1
ratio_gaim <- sum(gaim$N)/445*100
hete <- subset(df_N, Species_Code == "Heterostegina_depressa")
hete$N[hete$N > 0] <- 1
ratio_hete <- sum(hete$N)/445*100
plana <- subset(df_N, Species_Code == "Peneroplis_planatus")
plana$N[plana$N > 0] <- 1
ratio_plana <- sum(plana$N)/445*100
opp <- subset(df_N, Species_Code == "Operculina_ammonoides")
opp$N[opp$N > 0] <- 1
ratio_opp <- sum(opp$N)/445*100
calc <- subset(df_N, Species_Code == "Neorotalia_calcar")
calc$N[calc$N > 0] <- 1
ratio_calc <- sum(calc$N)/445*100
sori <- subset(df_N, Species_Code == "Sorites_orbiculus")
sori$N[sori$N > 0] <- 1
ratio_sori <- sum(sori$N)/445*100


#write new dataset with only all 17 most abundant species (without rare species) and with both true (N) and relative (RE) abundance
df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_withElphidiummerged.csv")

#To give the same format to this dataframe than the normalized one, but only dcast and melt - no change in data
df_N <- dcast(df, Sample_Code + Island + Sampling_Site + Sample_Year + 
                Zone + Habitat + Inhabitated + Distance_Makassar_km + Depth_m + Relative_Depth + Substrate + Depth_slope ~ 
                Genus + Species_Code, fun.aggregate = sum, value.var = "N")
df_N <- melt(df_N, id.vars = c("Sample_Year","Habitat","Zone","Island","Sampling_Site","Inhabitated","Distance_Makassar_km","Sample_Code","Substrate",
                               "Depth_m","Relative_Depth","Depth_slope"), 
             value.name = "N", variable.name = "Species_Code")

#calculation relative abundance
df_RE <- dcast(df, Sample_Code + Island + Sampling_Site + Sample_Year + 
                      Zone + Habitat + Inhabitated + Distance_Makassar_km + Depth_m + Relative_Depth + Substrate + Depth_slope ~ 
                        Genus + Species_Code, fun.aggregate = sum, value.var = "N")
df_west_1 <- apply(df_RE[13:ncol(df_RE)], 1, function(x) x/sum(x))
df_west_1 <- t(df_west_1)
df_RE <- cbind(df_RE[1:12], df_west_1)
df_RE <- melt(df_RE, id.vars = c("Sample_Code","Island","Sampling_Site","Sample_Year", 
                                               "Zone","Habitat","Inhabitated","Distance_Makassar_km","Depth_m","Relative_Depth","Substrate","Depth_slope"), 
                     value.name = "RE", variable.name = "Species_Code")

df <- cbind(df_N, df_RE$RE)
colnames(df)[15] <- "RE" #rename column of relative abundance to RE

ordered <- str_split_fixed(df$Species_Code, "_", 2)
colnames(ordered) <- c("Genus", "Species_Code1")

df <- cbind(df, ordered)
colnames(df)[13] <- "Genus_Species"
colnames(df)[17] <- "Species_Code"

#To remove rare species
df <- subset(df, subset = Species_Code %in% c("Amphisorus_SpL","Amphisorus_SpS", "Amphistegina_lessonii", "Amphistegina_lobifera", "Amphistegina_radiata", 
                                              "Calcarina_spengleri", "Calcarina_hispida", "Calcarina_mayori", "Elphidium_spp", "Heterostegina_depressa", 
                                              "Neorotalia_calcar", "Neorotalia_gaimardi", "Operculina_ammonoides", "Peneroplis_planatus", "Sorites_orbiculus", 
                                              "Baculogypsinoides_spinosus", "Sphaerogypsina_globulus"))
#to remove flat 2010 from the analysis
df <- df[!(df$Habitat == "flat" & df$Sample_Year == "2010"),]

df$Sampling_Time <- 0
df[df$Sample_Year <= 2013,18] <- "august"
df[df$Sample_Year >= 2018,18] <- "april"

write.csv(df, "~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv", row.names=FALSE)


#---------create working dataset with rare species in---------------------------

df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_withElphidiummerged.csv") #with rare species


df_N <- dcast(df, Sample_Code + Island + Sampling_Site + Sample_Year + 
                Zone + Habitat + Inhabitated + Distance_Makassar_km + Depth_m + Relative_Depth + Substrate + Depth_slope ~ 
                Genus + Species_Code, fun.aggregate = sum, value.var = "N")
df_N <- melt(df_N, id.vars = c("Sample_Year","Habitat","Zone","Island","Sampling_Site","Inhabitated","Distance_Makassar_km","Sample_Code","Substrate",
                               "Depth_m","Relative_Depth","Depth_slope"), 
             value.name = "N", variable.name = "Species_Code")

#calculation relative abundance
df_RE <- dcast(df, Sample_Code + Island + Sampling_Site + Sample_Year + 
                 Zone + Habitat + Inhabitated + Distance_Makassar_km + Depth_m + Relative_Depth + Substrate + Depth_slope ~ 
                 Genus + Species_Code, fun.aggregate = sum, value.var = "N")
df_west_1 <- apply(df_RE[13:ncol(df_RE)], 1, function(x) x/sum(x))
df_west_1 <- t(df_west_1)
df_RE <- cbind(df_RE[1:12], df_west_1)
df_RE <- melt(df_RE, id.vars = c("Sample_Code","Island","Sampling_Site","Sample_Year", 
                                 "Zone","Habitat","Inhabitated","Distance_Makassar_km","Depth_m","Relative_Depth","Substrate","Depth_slope"), 
              value.name = "RE", variable.name = "Species_Code")

df <- cbind(df_N, df_RE$RE)
colnames(df)[15] <- "RE" #rename column of relative abundance to RE

ordered <- str_split_fixed(df$Species_Code, "_", 2)
colnames(ordered) <- c("Genus", "Species_Code1")

df <- cbind(df, ordered)
colnames(df)[13] <- "Genus_Species"
colnames(df)[17] <- "Species_Code"

#to remove flat 2010 from the analysis
df <- df[!(df$Habitat == "flat" & df$Sample_Year == "2010"),]

df$Sampling_Time <- 0
df[df$Sample_Year <= 2013,18] <- "august"
df[df$Sample_Year >= 2018,18] <- "april"

write.csv(df, "~/Downloads/Dataset/20210905_LBFdataset_12islands_withrarespecies.csv", row.names=FALSE)



#---------To assess sample size and number of samples in all years and habitat and species distribution------------------------------------------------

df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv")

#sample size range - number of specimens counted per sample
df_samplesize <- df %>% group_by(Sample_Code) %>% summarise(Sum=sum(N))
summary(df_samplesize$Sum)

#number of samples per year
nb_all <- levels(droplevels(df$Sample_Code))
df_1997f <- subset(df, Sample_Year == "1997" & Habitat == "flat")
df_1997s <- subset(df, Sample_Year == "1997" & Habitat == "slope")
nb_1997f <- levels(droplevels(df_1997f$Sample_Code))
nb_1997s <- levels(droplevels(df_1997s$Sample_Code))
df_2010f <- subset(df, Sample_Year == "2010" & Habitat == "flat")
df_2010s <- subset(df, Sample_Year == "2010" & Habitat == "slope")
nb_2010f <- levels(droplevels(df_2010f$Sample_Code))
nb_2010s <- levels(droplevels(df_2010s$Sample_Code))
df_2012f <- subset(df, Sample_Year == "2012" & Habitat == "flat")
df_2012s <- subset(df, Sample_Year == "2012" & Habitat == "slope")
nb_2012f <- levels(droplevels(df_2012f$Sample_Code))
nb_2012s <- levels(droplevels(df_2012s$Sample_Code))
df_2013f <- subset(df, Sample_Year == "2013" & Habitat == "flat")
df_2013s <- subset(df, Sample_Year == "2013" & Habitat == "slope")
nb_2013f <- levels(droplevels(df_2013f$Sample_Code))
nb_2013s <- levels(droplevels(df_2013s$Sample_Code))
df_2018f <- subset(df, Sample_Year == "2018" & Habitat == "flat")
df_2018s <- subset(df, Sample_Year == "2018" & Habitat == "slope")
nb_2018f <- levels(droplevels(df_2018f$Sample_Code))
nb_2018s <- levels(droplevels(df_2018s$Sample_Code))

#---------species information and distribution------------------------------------
df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_withElphidiummerged.csv")
df$Zone <- as.character(df$Zone)


#To give the same format to this dataframe than the normalized one, but only dcast and melt - no change in data
df_N <- dcast(df, Sample_Code + Island + Sampling_Site + Sample_Year + 
                Zone + Habitat + Inhabitated + Distance_Makassar_km + Depth_m + Relative_Depth + Substrate + Depth_slope ~ 
                Genus + Species_Code, fun.aggregate = sum, value.var = "N")
df_N <- melt(df_N, id.vars = c("Sample_Year","Habitat","Zone","Island","Sampling_Site","Inhabitated","Distance_Makassar_km","Sample_Code","Substrate",
                               "Depth_m","Relative_Depth","Depth_slope"), 
             value.name = "N", variable.name = "Species_Code")

#calculation relative abundance
df_RE <- dcast(df, Sample_Code + Island + Sampling_Site + Sample_Year + 
                 Zone + Habitat + Inhabitated + Distance_Makassar_km + Depth_m + Relative_Depth + Substrate + Depth_slope ~ 
                 Genus + Species_Code, fun.aggregate = sum, value.var = "N")
df_west_1 <- apply(df_RE[13:ncol(df_RE)], 1, function(x) x/sum(x))
df_west_1 <- t(df_west_1)
df_RE <- cbind(df_RE[1:12], df_west_1)
df_RE <- melt(df_RE, id.vars = c("Sample_Code","Island","Sampling_Site","Sample_Year", 
                                 "Zone","Habitat","Inhabitated","Distance_Makassar_km","Depth_m","Relative_Depth","Substrate","Depth_slope"), 
              value.name = "RE", variable.name = "Species_Code")

df <- cbind(df_N, df_RE$RE)
colnames(df)[15] <- "RE" #rename column of relative abundance to RE

ordered <- str_split_fixed(df$Species_Code, "_", 2)
colnames(ordered) <- c("Genus", "Species_Code1")

df <- cbind(df, ordered)
colnames(df)[13] <- "Genus_Species"
colnames(df)[17] <- "Species_Code"

df <- df[!(df$Habitat == "flat" & df$Sample_Year == "2010"),]

df_N <- df
df1 <- df_N[1:nrow(df_N), c(3,10,14,17)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Species_Code, Zone, depthbins) %>% summarise(Sum = sum(N))

df_west_norm2 <- dcast(df1, Zone + depthbins ~ Species_Code, fun.aggregate = sum, value.var = "Sum")
df11 <- subset(df_west_norm2, Zone == "1")
df12 <- subset(df_west_norm2, Zone == "2")
df13 <- subset(df_west_norm2, Zone == "3")
df14 <- subset(df_west_norm2, Zone == "4")
df15 <- subset(df_west_norm2, Zone == "5")
df11 <- apply(df11[3:ncol(df11)], 2, function(x) x/sum(x))
df12 <- apply(df12[3:ncol(df12)], 2, function(x) x/sum(x))
df13 <- apply(df13[3:ncol(df13)], 2, function(x) x/sum(x))
df14 <- apply(df14[3:ncol(df14)], 2, function(x) x/sum(x))
df15 <- apply(df15[3:ncol(df15)], 2, function(x) x/sum(x))

df_All <- rbind(df11, df12, df13, df14, df15)
df_west_norm2 <- cbind(df_west_norm2[1:2], df_All)

df1 <- melt(df_west_norm2, id.vars = c("Zone", "depthbins"), 
             value.name = "Sum", variable.name = "Species_Code")
df1 <- df1[complete.cases(df1), ]


ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + 
  facet_grid(Zone~Species_Code)


alveo <- subset(df_N, Species_Code == "Alveolinella_quoyii")
mean(alveo$RE)
df1 <- alveo[1:nrow(alveo), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

dent <- subset(df_N, Species_Code == "Dendritina_ambigua")
mean(dent$RE)
df1 <- dent[1:nrow(dent), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

laev <- subset(df_N, Species_Code == "Laevipeneroplis_malayensis")
mean(laev$RE)
df1 <- laev[1:nrow(laev), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

numm <- subset(df_N, Species_Code == "Nummulites_venosus")
mean(numm$RE)
df1 <- numm[1:nrow(numm), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

par1 <- subset(df_N, Species_Code == "Parasorites_sp1")
mean(par1$RE)
df1 <- par1[1:nrow(par1), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

par2 <- subset(df_N, Species_Code == "Parasorites_sp2")
mean(par2$RE)
df1 <- par2[1:nrow(par2), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

pene <- subset(df_N, Species_Code == "Peneroplis_pertusus")
mean(pene$RE)
df1 <- pene[1:nrow(pene), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

penesp <- subset(df_N, Species_Code == "Peneroplis_sp")
mean(penesp$RE)
df1 <- penesp[1:nrow(penesp), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

penesp2 <- subset(df_N, Species_Code == "Peneroplis_sp2")
mean(penesp2$RE)
df1 <- penesp2[1:nrow(penesp2), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)


spl <- subset(df_N, Species_Code == "Amphisorus_SpL")
mean(spl$RE)
df1 <- spl[1:nrow(spl), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

sps <- subset(df_N, Species_Code == "Amphisorus_SpS")
mean(sps$RE)
df1 <- sps[1:nrow(sps), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

radia <- subset(df_N, Species_Code == "Amphistegina_radiata")
mean(radia$RE)
df1 <- radia[1:nrow(radia), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

lobi <- subset(df_N, Species_Code == "Amphistegina_lobifera")
mean(lobi$RE)
df1 <- lobi[1:nrow(lobi), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

speng <- subset(df_N, Species_Code == "Calcarina_spengleri")
mean(speng$RE)
df1 <- speng[1:nrow(speng), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

csp2 <- subset(df_N, Species_Code == "Calcarina_hispida")
mean(csp2$RE)
df1 <- csp2[1:nrow(csp2), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

mayo <- subset(df_N, Species_Code == "Calcarina_mayori")
mean(mayo$RE)
df1 <- mayo[1:nrow(mayo), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

less <- subset(df_N, Species_Code == "Amphistegina_lessonii")
mean(less$RE)
df1 <- less[1:nrow(less), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

elph <- subset(df_N, Species_Code == "Elphidium_spp")
mean(elph$RE)
df1 <- elph[1:nrow(elph), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

bacu <- subset(df_N, Species_Code == "Baculogypsinoides_spinosus")
mean(bacu$RE)
df1 <- bacu[1:nrow(bacu), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

gaim <- subset(df_N, Species_Code == "Neorotalia_gaimardi")
mean(gaim$RE)
df1 <- gaim[1:nrow(gaim), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

hete <- subset(df_N, Species_Code == "Heterostegina_depressa")
mean(hete$RE)
df1 <- hete[1:nrow(hete), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

plana <- subset(df_N, Species_Code == "Peneroplis_planatus")
mean(plana$RE)
df1 <- plana[1:nrow(plana), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

opp <- subset(df_N, Species_Code == "Operculina_ammonoides")
mean(opp$RE)
df1 <- opp[1:nrow(opp), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

calc <- subset(df_N, Species_Code == "Neorotalia_calcar")
mean(calc$RE)
df1 <- calc[1:nrow(calc), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

sori <- subset(df_N, Species_Code == "Sorites_orbiculus")
mean(sori$RE)
df1 <- sori[1:nrow(sori), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)

sphae <- subset(df_N, Species_Code == "Sphaerogypsina_globulus")
mean(sphae$RE)
df1 <- sphae[1:nrow(sphae), c(3,10,14)] 
df1$depthbins <- cut(df1$Depth_m, c(0,1,4,7,10,13,16,19,22,25,28,31), labels=c("1","2-4","5-7","8-10", "11-13","14-16","17-19","20-22","23-25","26-28","29-31"))
df1 <- df1 %>% group_by(Zone, depthbins) %>% summarise(Sum = sum(N))
ggplot(df1, aes(depthbins, Sum)) + 
  geom_bar(stat = "identity") + facet_grid(.~Zone)


#---------Fig. 3 and appendix Fig. A.4 - Time series-----------------------------------------

#____________________________Time series_______________________________________________________________
df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv") #without rare species
df <- read.csv("~/Downloads/Dataset/20210905_LBFdataset_12islands_withrarespecies.csv") #with rare species

df$Zone <- as.character(as.numeric(df$Zone))
df$Species_Code <- as.factor(as.character(df$Species_Code))

#Average relative abundance of a species per year per habitat
df1 <- df %>% group_by(Habitat, Sample_Year, Species_Code) %>% summarise(Mean = mean(RE))
df_peryear <- dcast(df1, Habitat + Species_Code ~ Sample_Year, value.var = "Mean")
#remove absent species over all years (zero abundance every year)
df_peryear <- df_peryear[apply(df_peryear[, 3:7], 1, function(x) any(x != 0)),]
ind <- apply(df_peryear, 1, function(x) sum(is.na(x))>2)
df_peryear <- df_peryear[!ind,]
df_peryear$slope <- 0
#to calculate slope of each regression species ~ island - Reef flat and slope
for (x in 1:nrow(df_peryear)){
  years <- c(1, 13, 15, 16, 21) #vector associated to years
  values <- as.numeric(paste(df_peryear[x, 3:7])) #vector associated to mean relative abundance each year
  reg <- lm(values ~ years)
  df_peryear[x,8] <- as.numeric(paste(reg$coef[[2]]))
}

df_peryear$slopesign <- NA
df_peryear$slopesign[df_peryear$slope > 0] <- "positive"
df_peryear$slopesign[df_peryear$slope < 0] <- "negative"

df1$slopesign <- NA

for (i in 1:nrow(df_peryear)) {
  
  df1$slopesign[df1$Species_Code == df_peryear[i,2] & df1$Habitat == df_peryear[i,1]] <- df_peryear[i,9]
  
}

#add standard error range around the mean of each species
strd_error <- function(x) sd(x)/sqrt(length(x))

#calculate standard error of the species mean relative abundance
df_SE <- df %>% group_by(Habitat, Sample_Year, Species_Code) %>% summarise(SE = strd_error(RE))

df1$SE <- df_SE$SE

df2 <- na.omit(df1)

colors <- c('#E6194B','#F58231','#FFE119','#FABED4','#FFD8B1','#BFEF45','#3CB44B','#42D4F4',
            '#000075','#4363D8','#AAFFC3','#911EB4','#F032E6','#A9A9A9','#800000','#DCBEFF','#FFFAC8')

levels(droplevels(df2$Species_Code))


ggplot(df2, aes(Sample_Year, Mean, group = Species_Code)) + 
  geom_ribbon(aes(ymin=Mean-SE, ymax=Mean+SE, fill = Species_Code), alpha=0.25)+
  geom_line(aes(col=Species_Code), size = 1) + geom_point(aes(col=Species_Code), size = 2) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + 
  facet_grid(slopesign~Habitat) + theme_bw()

#____________________________Pie charts_______________________________________________________________

df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv")

#to give the right color to all species as in the bar chart
colors <- c('#E6194B','#F58231','#FFE119','#FABED4','#FFD8B1','#BFEF45','#3CB44B','#42D4F4',
            '#FFFAC8','#4363D8','#AAFFC3','#911EB4','#F032E6','#A9A9A9','#800000','#DCBEFF','#000075')

df_flat <- subset(df, subset = Habitat %in% "flat")
df_slope <- subset(df, subset = Habitat %in% "slope")

#transform dataframe to get: Species code, and relative abundance for all years as column
flat <- df_flat %>% group_by(Sample_Year, Species_Code) %>% summarise(Mean = mean(RE))
flat <- dcast(flat, Species_Code ~ Sample_Year, fun.aggregate = sum, value.var = "Mean")

slope <- df_slope %>% group_by(Sample_Year, Species_Code) %>% summarise(Mean = mean(RE))
slope <- dcast(slope, Species_Code ~ Sample_Year, fun.aggregate = sum, value.var = "Mean")

#print 2 row of 5 piecharts
par(mfrow=c(2,5))
#to make a piechart slope 
lbls <- slope$Species_Code

slices <- slope$'1997'
pct <- round(slices/sum(slices)*100)
lbls1 <- paste(lbls, pct) # add percents to labels
lbls1 <- paste(lbls1,"%",sep="") # ad % to labels
pie(slices,labels = lbls1, col=colors,
    main="slope 1997")
slices <- slope$'2010'
pct <- round(slices/sum(slices)*100)
lbls2 <- paste(lbls, pct) # add percents to labels
lbls2 <- paste(lbls2,"%",sep="") # ad % to labels
pie(slices,labels = lbls2, col=colors,
    main="slope 2010")
slices <- slope$'2012'
pct <- round(slices/sum(slices)*100)
lbls3 <- paste(lbls, pct) # add percents to labels
lbls3 <- paste(lbls3,"%",sep="") # ad % to labels
pie(slices,labels = lbls3, col=colors,
    main="slope 2012")
slices <- slope$'2013'
pct <- round(slices/sum(slices)*100)
lbls4 <- paste(lbls, pct) # add percents to labels
lbls4 <- paste(lbls4,"%",sep="") # ad % to labels
pie(slices,labels = lbls4, col=colors,
    main="slope 2013")
slices <- slope$'2018'
pct <- round(slices/sum(slices)*100)
lbls5 <- paste(lbls, pct) # add percents to labels
lbls5 <- paste(lbls5,"%",sep="") # ad % to labels
pie(slices,labels = lbls5, col=colors,
    main="slope 2018")

#to make a piechart flat
lbls <- flat$Species_Code

slices <- flat$'1997'
pct <- round(slices/sum(slices)*100)
lbls6 <- paste(lbls, pct) # add percents to labels
lbls6 <- paste(lbls6,"%",sep="") # ad % to labels
pie(slices,labels = lbls6, col=colors,
    main="flat 1997")
slices <- flat$'2012'
pct <- round(slices/sum(slices)*100)
lbls8 <- paste(lbls, pct) # add percents to labels
lbls8 <- paste(lbls8,"%",sep="") # ad % to labels
pie(slices,labels = lbls8, col=colors,
    main="flat 2012")
slices <- flat$'2013'
pct <- round(slices/sum(slices)*100)
lbls9 <- paste(lbls, pct) # add percents to labels
lbls9 <- paste(lbls9,"%",sep="") # ad % to labels
pie(slices,labels = lbls9, col=colors,
    main="flat 2013")
slices <- flat$'2018'
pct <- round(slices/sum(slices)*100)
lbls11 <- paste(lbls, pct) # add percents to labels
lbls11 <- paste(lbls11,"%",sep="") # ad % to labels
pie(slices,labels = lbls11, col=colors,
    main="flat 2018")

#---------Fig. 4 - co-occurrence matrix--------------------------------
df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv")

df$Sample_Year <- as.factor(as.numeric(df$Sample_Year))
df$Zone <- as.factor(as.numeric(df$Zone))
df$Species_Code <- as.factor(as.character(df$Species_Code))
df$Habitat <- as.factor(df$Habitat)
df$Sampling_Site <- as.factor(df$Sampling_Site)
df$Sample_Code <- as.factor(df$Sample_Code)
df$Substrate <- as.factor(df$Substrate)
df$N <- as.numeric(df$N)
df$Inhabitated <- as.factor(df$Inhabitated)
df$Sampling_Time <- as.factor(df$Sampling_Time)

lev <- c("Amphisorus_SpL","Amphisorus_SpS","Elphidium_spp",
         "Amphistegina_lobifera",  "Calcarina_hispida","Neorotalia_calcar", "Neorotalia_gaimardi", "Peneroplis_planatus", "Sorites_orbiculus", 
         "Amphistegina_lessonii", "Amphistegina_radiata",  "Baculogypsinoides_spinosus", "Calcarina_mayori","Calcarina_spengleri", "Heterostegina_depressa", 
         "Operculina_ammonoides", "Sphaerogypsina_globulus")

df <- df %>%
  mutate(Species_Code =  factor(Species_Code, levels = lev)) %>%
  arrange(Species_Code)  

dff <- subset(df, Habitat == "flat")
dfs <- subset(df, Habitat == "slope")

df1997 <- subset(dfs, Sample_Year == "1997")
df2010 <- subset(dfs, Sample_Year == "2010")
df2012 <- subset(dfs, Sample_Year == "2012")
df2013 <- subset(dfs, Sample_Year == "2013")
df2018 <- subset(dfs, Sample_Year == "2018")

#make matrix
df_NMDS <- dcast(df, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat + Inhabitated + Sampling_Time ~ Species_Code, fun.aggregate = sum, value.var = "N")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(df, Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$Sample_Code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
countmatrix <- t(df_NMDS_vegan)
presabsmatrix <- countmatrix
presabsmatrix[presabsmatrix > 0] <- 1

x <- cooccur(presabsmatrix, spp_names = TRUE)

plotforcooccur(x)


df_NMDS <- dcast(df1997, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat + Inhabitated + Sampling_Time ~ Species_Code, fun.aggregate = sum, value.var = "N")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(df1997, Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$Sample_Code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
countmatrix <- t(df_NMDS_vegan)
presabsmatrix <- countmatrix
presabsmatrix[presabsmatrix > 0] <- 1

x <- cooccur(presabsmatrix, spp_names = TRUE)

plotforcooccur(x)

df_NMDS <- dcast(df2010, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat + Inhabitated + Sampling_Time ~ Species_Code, fun.aggregate = sum, value.var = "N")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(df2010, Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$Sample_Code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
countmatrix <- t(df_NMDS_vegan)
presabsmatrix <- countmatrix
presabsmatrix[presabsmatrix > 0] <- 1

x <- cooccur(presabsmatrix, spp_names = TRUE)

plotforcooccur(x)

df_NMDS <- dcast(df2012, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat + Inhabitated + Sampling_Time ~ Species_Code, fun.aggregate = sum, value.var = "N")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(df2012, Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$Sample_Code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
countmatrix <- t(df_NMDS_vegan)
presabsmatrix <- countmatrix
presabsmatrix[presabsmatrix > 0] <- 1

x <- cooccur(presabsmatrix, spp_names = TRUE)

plotforcooccur(x)

df_NMDS <- dcast(df2013, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat + Inhabitated + Sampling_Time ~ Species_Code, fun.aggregate = sum, value.var = "N")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(df2013, Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$Sample_Code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
countmatrix <- t(df_NMDS_vegan)
presabsmatrix <- countmatrix
presabsmatrix[presabsmatrix > 0] <- 1

x <- cooccur(presabsmatrix, spp_names = TRUE)

plotforcooccur(x)

df_NMDS <- dcast(df2018, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat + Inhabitated + Sampling_Time ~ Species_Code, fun.aggregate = sum, value.var = "N")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(df2018, Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$Sample_Code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
countmatrix <- t(df_NMDS_vegan)
presabsmatrix <- countmatrix
presabsmatrix[presabsmatrix > 0] <- 1

x <- cooccur(presabsmatrix, spp_names = TRUE)

plotforcooccur(x)














#---------Fig. 5 and appendix Fig. A.5 - Bubble charts--------------------------

df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv")
df <- read.csv("~/Downloads/Dataset/20210905_LBFdataset_12islands_withrarespecies.csv") #with rare species


df_west_norm1 <- df[, c(1, 2, 9, 15, 17)]

#to compare only 1997, 2012, 2018 (2010 has no reef flat)
df_west_norm1 <- subset(df_west_norm1, subset = Substrate %in% c("RG", "RH", "R", "S", "RS" ,"RA", "GA", "RSA", "A"))
df_west_norm1 <- subset(df_west_norm1, subset = Sample_Year %in% c("1997","2010", "2012", "2013", "2018"))

#summarizing to the mean is the best option to represent the ratio which species to which substrate
df_reordered1 <- df_west_norm1 %>% group_by(Substrate, Species_Code, Sample_Year, Habitat) %>% summarise(RE = mean(RE)) 
#order levels for the alluvial chart

df_reordered1 <- df_reordered1 %>% mutate(Substrate = fct_relevel(as.factor(Substrate), c("RG", "RH", "R", "S", "RS" ,"RA", "RSA", "GA","A")))
df_reordered1 <- df_reordered1 %>% mutate(Sample_Year = fct_relevel(as.factor(Sample_Year), c("1997","2010", "2012", "2013", "2018")))

#plotting with alluvial
df_2 <- df_reordered1[,c(3,1,2,4,5)]

#plotting with alluvial

p1 <- ggballoonplot(df_2, x = "Substrate", y = "Species_Code", size = "RE",
              fill = "Substrate",color = "Substrate", facet.by = c("Habitat", "Sample_Year"),
              ggtheme = theme_bw()) +  
  scale_fill_manual(values = c("aquamarine3", "aquamarine4","coral1","antiquewhite","bisque3", "darkseagreen1", "palegreen2","seagreen3","springgreen4")) +
  scale_color_manual(values = c("aquamarine3", "aquamarine4","coral1","antiquewhite","bisque3", "darkseagreen1", "palegreen2","seagreen3","springgreen4")) +
  theme(axis.text.x = element_text(angle = 90), panel.grid.major.x = element_blank())


df$Sample_Year <- as.character(as.numeric(df$Sample_Year))

df_substrate <- df 
df_substrate$counts <- 0
df_substrate <- df_substrate %>% group_by(Sample_Year, Habitat, Substrate, Sample_Code) %>% summarize(counts = sum(counts))
df_substrate$counts <- 1
df_ss <- df_substrate %>% group_by(Sample_Year, Habitat, Substrate) %>% summarize(counts = sum(counts))

#to remove the empty row with no identified substrate
df_ss <- df_ss[-13,]

df_ss <- dcast(df_ss, Sample_Year + Habitat ~ Substrate,  value.var = "counts")
df_ss[is.na(df_ss)] <- 0
df_ss1 <- t(apply(df_ss[3:ncol(df_ss)], 1, function(x) x/sum(x)))
colSums(df_ss1)
rowSums(df_ss1)
df_ss[,3:ncol(df_ss)] <- df_ss1

df_ss <- melt(df_ss, id.vars = c("Sample_Year", "Habitat"), variable.name = "Substrate" , value.name = "counts")

#to order the substrates
df_ss <- df_ss %>% mutate(Substrate = fct_relevel(as.factor(Substrate), c("RG", "RH", "R","S", "RS" ,"RA", "RSA", "GA", "A")))

p2 <- ggplot(df_ss, aes(x=Substrate, y=counts)) + 
  geom_bar(aes(fill = Substrate), position="dodge", stat="identity") + facet_grid(Habitat ~ Sample_Year) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), panel.grid.major.x = element_blank())+
  scale_fill_manual(values = c("aquamarine3", "aquamarine4","coral1","antiquewhite","bisque3", "darkseagreen1", "palegreen2","seagreen3","springgreen4")) +
  ggtitle("Substrate through time") +
  xlab("Sampling year") + ylab("Number of time a substrate was sampled")


lay <- rbind(c(1,1,1),
             c(1,1,1),
             c(2,2,2))

gridExtra::grid.arrange(p1,p2, layout_matrix = lay)



#---------Appendix Fig. A.1 and Tab A.2 - NMDS and ANOSIM--------------------------------------------------------

df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv")

df$Sample_Year <- as.factor(as.numeric(df$Sample_Year))
df$Zone <- as.factor(as.numeric(df$Zone))
df$Species_Code <- as.factor(as.character(df$Species_Code))
df$Habitat <- as.factor(df$Habitat)
df$Sampling_Site <- as.factor(df$Sampling_Site)
df$Sample_Code <- as.factor(df$Sample_Code)
df$Substrate <- as.factor(df$Substrate)
df$N <- as.numeric(df$N)
df$Inhabitated <- as.factor(df$Inhabitated)
df$Sampling_Time <- as.factor(df$Sampling_Time)

#Select slope or flat
df <- subset(df, Habitat == "flat")

#_______________________habitat differences in LBF____________________________________

#To build the NMDS, I have to "unmelt" the dataframe with reshape
df_NMDS <- dcast(df, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat + Inhabitated + Sampling_Time ~ Species_Code, fun.aggregate = sum, value.var = "N")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(df, Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$Sample_Code
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df_NMDS[,8:ncol(df_NMDS)]

#groups for plots
zones <- df_NMDS$Zone
years <- df_NMDS$Sample_Year
islands <- df_NMDS$Sampling_Site
habitat <- df_NMDS$Habitat
inhabitated <- df_NMDS$Inhabitated
time <- df_NMDS$Sampling_Time
colszone <- c("red", "orange", "yellow", "green", "blue") 
colsyear <- c("red", "orange", "yellow", "green", "blue")
colsisland <- c("red", "magenta", "red4",
                "darkorange", "orangered", "orangered3",
                "green", "seagreen", "green3",
                "cyan", "blue","purple")
colhabitat <- c("red", "blue")
colinhabitated <- c("red", "blue")
coltime <- c("red", "blue")

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, habitat) 
#for plot
plot(ord, disp="sites", type = "n", main = paste ("Slope vs Flat without filters - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordiellipse(ord, habitat, col= colhabitat, conf = 0.95)
points(ord, disp="sites", pch= c(19), bg=as.numeric(habitat), cex=0.7, col = colhabitat[habitat])
with(ord, legend(x = "topright", legend = levels(habitat), col = colhabitat, pch = c(19)))
plot(en, col = "black")

ano <- anosim(df_NMDS_vegan, years) 
plot(ord, disp="sites", type = "n", main = paste ("Both habitats - years - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordiellipse(ord, years, col= colsyear, conf = 0.95)
points(ord, disp="sites", pch= c(19), bg=as.numeric(years), cex=0.7, col = colsyear[years])
with(ord, legend(x = "topright", legend = levels(years), col = colsyear, pch = c(19)))
plot(en, col = "black")

ano <- anosim(df_NMDS_vegan, inhabitated) 
plot(ord, disp="sites", type = "n", main = paste ("Flat - Inhabitancy - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordiellipse(ord, inhabitated, col= colinhabitated, conf = 0.95)
points(ord, disp="sites", pch= c(19), bg=as.numeric(inhabitated), cex=0.7, col = colinhabitated[inhabitated])
with(ord, legend(x = "topright", legend = levels(inhabitated), col = colinhabitated, pch = c(19)))
plot(en, col = "black")

ano <- anosim(df_NMDS_vegan, zones) 
plot(ord, disp="sites", type = "n", main = paste ("Slope zones - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordiellipse(ord, zones, col= colszone, conf = 0.95)
points(ord, disp="sites", pch= c(19), bg=as.numeric(zones), cex=0.7, col = colszone[zones])
with(ord, legend(x = "topright", legend = levels(zones), col = colszone, pch = c(19)))
plot(en, col = "black")

ano <- anosim(df_NMDS_vegan, time) 
plot(ord, disp="sites", type = "n", main = paste ("Flat - Sampling time - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordiellipse(ord, time, col= coltime, conf = 0.95)
points(ord, disp="sites", pch= c(19), bg=as.numeric(time), cex=0.7, col = coltime[time])
with(ord, legend(x = "topright", legend = levels(time), col = coltime, pch = c(19)))
plot(en, col = "black")


#_________________Flat, slope separately - in order to attest per zone the significance per year_________________________________________

#change the zone number (1 to 5) for different analysis
df_flat <- subset(df, Habitat == "flat" & Zone == "5")
df_slope <- subset(df, Habitat == "slope" & Zone == "5")

#flat
df_NMDS_flat <- dcast(df_flat, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat ~ Species_Code, fun.aggregate = sum, value.var = "N")
df_NMDS_vegan <- dcast(df_flat, Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$Sample_Code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total")
years_flat <- df_NMDS_flat$Sample_Year
anosim(df_NMDS_vegan, years_flat) 

#slope
df_NMDS_slope <- dcast(df_slope, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat ~ Species_Code, fun.aggregate = sum, value.var = "N")
df_NMDS_vegan <- dcast(df_slope, Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$Sample_Code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total")
years_slope <- df_NMDS_slope$Sample_Year
anosim(df_NMDS_vegan, years_slope) 



#---------Appendix Fig. A.2 - Community diversity indexes-------------------------------

#identify top 3 most abundant species per sample, make a list and isolate the species to make a new dataset.
#replace dataset in df <- XXX for different analysis
df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv")

df1 <- dcast(data = df, Sample_Code + Sample_Year + Island + Zone + Habitat + Substrate ~ Species_Code, fun.aggregate = sum, value.var = "RE")

#Community informations on the samples
df2 <- df1
df2$Richness_s <- apply(df1[,7:ncol(df1)]>0,1,sum)
df2$Abundance <- apply(df1[,7:ncol(df1)],1,sum)
#Diversity indexes
df2$Shannon_H <- vegan::diversity(df1[,7:ncol(df1)], index="shannon")
df2$Simpson_D <- vegan::diversity(df1[,7:ncol(df1)], index="simpson")
df2$Evenness_J <- vegan::diversity(df1[,7:ncol(df1)], index="simpson")/log(df2$Richness)
#True diversiy indexes
df2$True_shannon <- exp(df2$Shannon_H)
df2$True_simpson <- 1/df2$Simpson_D

#melt dataset to plot diversity index per habitat, per year, per index
df3 <- melt(df2, id.vars=c("Sample_Year", "Habitat", "Zone", "Island", "Sample_Code", "Substrate", "Richness_s", "Abundance", "Shannon_H", 
                           "Simpson_D", "Evenness_J", "True_shannon", "True_simpson"),
            variable.name = "Species_Code", value.name = "N")
df4 <- melt(df3, id.vars=c("Sample_Year", "Habitat", "Zone", "Island", "Sample_Code", "Substrate", "Species_Code", "N"),
            variable.name = "Diversity_Index", value.name = "Value")

df5 <- subset(df4, subset = Diversity_Index %in% c("Richness_s", "Shannon_H", 
                                                   "Simpson_D", "Evenness_J", "True_shannon", "True_simpson"))

df6 <- df5 %>% group_by(Sample_Year, Zone, Habitat, Diversity_Index) %>% summarise(Mean = mean(Value))
df7 <- df5 %>% group_by(Sample_Year, Zone, Habitat, Diversity_Index) %>% summarise(SD = sd(Value))
df7$Mean <- df6$Mean

df7$Sample_Year <- as.character(as.numeric(df7$Sample_Year))
df5$Sample_Year <- as.character(as.numeric(df5$Sample_Year))

#plotting
ggplot(df5, aes(x=Sample_Year, y=Value)) + 
  geom_boxplot() +
  #scale_y_continuous(trans = "log10")+
  facet_grid(Diversity_Index~Habitat, scales="free_y") + theme_bw()

ggplot(df7, aes(x=Sample_Year, y=Mean, group = Diversity_Index)) + 
  geom_line(aes(col=Diversity_Index), size=1.5,position=position_dodge(0.5)) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD, col=Diversity_Index), width=0.5, size = 0.5,
                position=position_dodge(0.5),linetype = "dashed") +
  facet_grid(Habitat~Zone) + 
  theme_bw() +
  ggtitle("Diversity indexes without rare species") +
  xlab("Sampling year") + ylab("Index value")


#---------Appendix Fig. A.3 - Heatmap-----------------------------------
df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv")

df_flat <- subset(df, subset = Habitat %in% "flat")
df_slope <- subset(df, subset = Habitat %in% "slope")

df_flat <- subset(df_flat, subset = Species_Code %in% 
                    c("Amphistegina_lobifera","Sorites_orbiculus","Peneroplis_planatus",   
                      "Neorotalia_calcar","Calcarina_hispida","Neorotalia_gaimardi"))
df_slope <- subset(df_slope, subset = Species_Code %in% 
                     c("Amphistegina_radiata","Calcarina_spengleri","Heterostegina_depressa","Sphaerogypsina_globulus",  
                       "Operculina_ammonoides","Baculogypsinoides_spinosus", "Amphistegina_lessonii","Calcarina_mayori"))

#create mean relative abundance of each species according to island and year
groups_flat <- df_flat %>% group_by(Zone, Sampling_Site, Sample_Year, Species_Code) %>% summarise(Mean = mean(RE))
groups_slope <- df_slope %>% group_by(Zone, Sampling_Site, Sample_Year, Species_Code) %>% summarise(Mean = mean(RE))

df_peryear_flat <- dcast(groups_flat, Zone + Sampling_Site + Species_Code ~ 
                           Sample_Year, value.var = "Mean")
df_peryear_slope <- dcast(groups_slope, Zone + Sampling_Site + Species_Code ~ 
                            Sample_Year, value.var = "Mean")

#remove absent species over all years (zero abundance every year)
df_peryear_flat <- df_peryear_flat[apply(df_peryear_flat[, 4:7], 1, function(x) any(x != 0)),]
ind <- apply(df_peryear_flat, 1, function(x) sum(is.na(x))>2)
df_peryear_flat <- df_peryear_flat[!ind,]

df_peryear_slope <- df_peryear_slope[apply(df_peryear_slope[, 4:8], 1, function(x) any(x != 0)),]
ind <- apply(df_peryear_slope, 1, function(x) sum(is.na(x))>3)
df_peryear_slope <- df_peryear_slope[!ind,]

df_peryear_flat$slope <- 0
df_peryear_slope$slope <- 0

#to calculate slope of each regression species ~ island - Reef flat and slope

for (x in 1:nrow(df_peryear_flat)){
  years <- c(1, 15, 16, 21) #vector associated to years
  values <- as.numeric(paste(df_peryear_flat[x, 4:7])) #vector associated to mean relative abundance each year
  
  reg <- lm(values ~ years)
  
  df_peryear_flat[x,8] <- as.numeric(paste(reg$coef[[2]]))
  
  
}

for (x in 1:nrow(df_peryear_slope)){
  years <- c(1, 13, 15, 16, 21) #vector associated to years
  values <- as.numeric(paste(df_peryear_slope[x, 4:8])) #vector associated to mean relative abundance each year
  
  reg <- lm(values ~ years)
  
  df_peryear_slope[x,9] <- as.numeric(paste(reg$coef[[2]]))
  
  
}

#Heatmap 
order_island <- c("LL1","BN1", "LD2", "KG2","PO2", "SA3","LO3", "PJ3", "KU4", "LU4", "BA4", "LA5")

#species order for SLOPE
order_species_slope <- c("Amphistegina_radiata","Calcarina_spengleri","Heterostegina_depressa","Sphaerogypsina_globulus",  
                         "Operculina_ammonoides","Baculogypsinoides_spinosus", "Amphistegina_lessonii","Calcarina_mayori")

#species order for FLAT
order_species_flat <- c("Amphistegina_lobifera","Sorites_orbiculus","Peneroplis_planatus",   
                        "Neorotalia_calcar","Calcarina_hispida","Neorotalia_gaimardi")

#plot flat regression slope
df_heatmap <- dcast(df_peryear_flat, Sampling_Site ~ Species_Code, value.var = "slope")
rownames(df_heatmap) <- df_heatmap$Sampling_Site
df_heatmap <- df_heatmap[,-1]
#df_heatmap <- df_heatmap[order_island, order_species_flat]
df_heatmap <- df_heatmap[order_island,]

df_heatmap <- as.matrix(df_heatmap)
df_heatmap <- t(df_heatmap)

#set up flat colors
myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(30)

#to see where the NA are
pheatmap(df_heatmap, border_color=NA,color = myColor,
         cluster_cols = FALSE, cluster_rows = FALSE, na_col = "black")

#after finding where the NAs are, we can transform the data with NAs into 0 and then edit in illustrator
df_heatmap[is.na(df_heatmap)] <- 0 #NAs into zeros
myBreaks <- c(seq(min(df_heatmap), 0, length.out=ceiling(15) + 1), 
              seq(max(df_heatmap)/30, max(df_heatmap), length.out=floor(15)))

pheatmap(df_heatmap, border_color=NA,color = myColor, breaks = myBreaks,
         cluster_cols = FALSE, cluster_rows = FALSE, na_col = "black")

#plot slope regression slope
df_heatmap <- dcast(df_peryear_slope, Sampling_Site ~ Species_Code, value.var = "slope")
rownames(df_heatmap) <- df_heatmap$Sampling_Site
df_heatmap <- df_heatmap[,-1]
#df_heatmap <- df_heatmap[order_island, order_species_slope]
df_heatmap <- df_heatmap[order_island,]
df_heatmap <- as.matrix(df_heatmap)
df_heatmap <- t(df_heatmap)

myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(30)
pheatmap(df_heatmap, border_color=NA, color = myColor,
         cluster_cols = FALSE, cluster_rows = FALSE, na_col = "black")

df_heatmap[is.na(df_heatmap)] <- 0 #NAs into zeros
myBreaks <- c(seq(min(df_heatmap), 0, length.out=ceiling(15) + 1), 
              seq(max(df_heatmap)/30, max(df_heatmap), length.out=floor(15)))

pheatmap(df_heatmap, border_color=NA, color = myColor, breaks = myBreaks,
         cluster_cols = FALSE, cluster_rows = FALSE, na_col = "black")


#to calculate the constant a in y = a + bx of each regression species ~ island
df_peryear_flat$constant <- 0
df_peryear_slope$constant <- 0

for (x in 1:nrow(df_peryear_flat)){
  years <- c(1, 15, 16, 21) #vector associated to years
  values <- as.numeric(paste(df_peryear_flat[x, 4:7])) #vector associated to mean relative abundance each year
  
  reg <- lm(values ~ years)
  
  df_peryear_flat[x,9] <- as.numeric(paste(reg$coef[[1]]))
  
  
}

for (x in 1:nrow(df_peryear_slope)){
  years <- c(1, 13, 15, 16, 21) #vector associated to years
  values <- as.numeric(paste(df_peryear_slope[x, 4:8])) #vector associated to mean relative abundance each year
  
  reg <- lm(values ~ years)
  
  df_peryear_slope[x,10] <- as.numeric(paste(reg$coef[[1]]))
  
  
}


#to calculate correlation coefficient from Pearson's formula
df_peryear_flat$pearson <- 0
df_peryear_slope$pearson <- 0

for (x in 1:nrow(df_peryear_flat)){
  years <- c(1, 15, 16, 21)
  abundance <- c()
  
  abundance <- as.numeric(df_peryear_flat[x, 4:7])
  corcoef <- cor(years, abundance, use = "complete.obs", method = c("pearson"))
  
  df_peryear_flat$pearson[x] <- corcoef
  
}

for (x in 1:nrow(df_peryear_slope)){
  years <- c(1, 13, 15, 16, 21)
  abundance <- c()
  
  abundance <- as.numeric(df_peryear_slope[x, 4:8])
  corcoef <- cor(years, abundance, use = "complete.obs", method = c("pearson"))
  
  df_peryear_slope$pearson[x] <- corcoef
  
}


#plot flat pearson
df_heatmap <- dcast(df_peryear_flat, Sampling_Site ~ Species_Code, value.var = "pearson")
rownames(df_heatmap) <- df_heatmap$Sampling_Site
df_heatmap <- df_heatmap[,-1]
df_heatmap <- df_heatmap[order_island, order_species_flat]

df_heatmap <- as.matrix(df_heatmap)
df_heatmap <- t(df_heatmap)

pheatmap(df_heatmap, border_color=NA,color = colorRampPalette(c("darkorange", "white", "darkorange"))(40),
         cluster_cols = FALSE, cluster_rows = FALSE, na_col = "black")

#plot slope pearson
df_heatmap <- dcast(df_peryear_slope, Sampling_Site ~ Species_Code, value.var = "pearson")
rownames(df_heatmap) <- df_heatmap$Sampling_Site
df_heatmap <- df_heatmap[,-1]
df_heatmap <- df_heatmap[order_island, order_species_slope]

df_heatmap <- as.matrix(df_heatmap)
df_heatmap <- t(df_heatmap)

pheatmap(df_heatmap, border_color=NA,color = colorRampPalette(c("darkorange", "white", "darkorange"))(40),
         cluster_cols = FALSE, cluster_rows = FALSE, na_col = "black")


#---------Appendix Tab A.3 - Indicator species analysis---------------------------------------------------

df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv")

#df <- subset(df, Habitat == "flat")

#_____________________Flat vs slope, years, zones _______________________________

df_indicator <- dcast(df, Sample_Code + Sample_Year + Sampling_Site + Zone + Habitat + Depth_slope + Substrate + Inhabitated + Sampling_Time ~ Species_Code, 
                      fun.aggregate = sum, value.var = "RE")
df_species <- df_indicator[,-c(2:9)]
rownames(df_species) <- df_species$Sample_Code
df_species <- df_species[,-1]

year <- df_indicator$Sample_Year
island <- df_indicator$Sampling_Site
habitat <- df_indicator$Habitat
zone <- df_indicator$Zone
depth <- df_indicator$Depth_slope
subs <- df_indicator$Substrate
time <- df_indicator$Sampling_Time
inhab <- df_indicator$Inhabitated

#Indicator species analysis using IndVal.g - ABUNDANCE
#Indicator species analysis for habitat
#default func = "IndVal.g" --> it takes into account inequal sample size (".g")
indval <- multipatt(df_species, habitat, control = how(nperm=999))
summary(indval)
indval <- multipatt(df_species, zone, control = how(nperm=999))
summary(indval)
indval <- multipatt(df_species, year, control = how(nperm=999))
summary(indval)
indval <- multipatt(df_species, subs, control = how(nperm=999))
summary(indval)
indval <- multipatt(df_species, time, control = how(nperm=999))
summary(indval)
indval <- multipatt(df_species, inhab, control = how(nperm=999))
summary(indval)


#______________________Separated habitat - per year, per zones____________________________________

slope <- subset(df, Habitat == "slope")
flat <- subset(df, Habitat == "flat")

slope <- subset(slope, subset = Sampling_Site %not% "PJ3")

#"unmelt" the dataframe with reshape and
#to isolate values (N), Species code (x) and Sample codes (y) with col and row names
slope_indicator <- dcast(slope, Sample_Code + Sample_Year + Sampling_Site + Zone + Substrate + Depth_slope ~ Species_Code, 
                         fun.aggregate = sum, value.var = "RE")
slope_species <- slope_indicator[,-c(2:6)]
rownames(slope_species) <- slope_species$Sample_Code
slope_species <- slope_species[,-1]
year_slope <- slope_indicator$Sample_Year
subs_slope <- slope_indicator$Substrate
zone_slope <- slope_indicator$Zone
depth_slope <- slope_indicator$Depth_slope

flat_indicator <- dcast(flat, Sample_Code + Sample_Year + Sampling_Site + Zone + Substrate + Depth_slope ~ Species_Code, 
                        fun.aggregate = sum, value.var = "RE")
flat_species <- flat_indicator[,-c(2:6)]
rownames(flat_species) <- flat_species$Sample_Code
flat_species <- flat_species[,-1]
year_flat <- flat_indicator$Sample_Year
subs_flat <- flat_indicator$Substrate
zone_flat <- flat_indicator$Zone

#Indicator species analysis for years
#max.order = 1
indval <- multipatt(flat_species, year_flat, max.order = 1, control = how(nperm=999))
summary(indval,  indvalcomp=TRUE)
indval <- multipatt(flat_species, zone_flat,max.order = 1, control = how(nperm=999))
summary(indval)
indval <- multipatt(flat_species, subs_flat,max.order = 1, control = how(nperm=999))
summary(indval)


indval <- multipatt(slope_species, year_slope,max.order = 1, control = how(nperm=999))
summary(indval, indvalcomp=TRUE)
indval <- multipatt(slope_species, zone_slope,max.order = 1, control = how(nperm=999))
summary(indval)
indval <- multipatt(slope_species, subs_slope,max.order = 1, control = how(nperm=999))
summary(indval)
indval <- multipatt(slope_species, depth_slope, control = how(nperm=999))
summary(indval)



#---------Appendix Tab A.4 - ANOVA analysis-----------------------------------------

df <- read.csv("~/Downloads/Dataset/20210217_LBFdataset_12islands_17mostabundantspecies.csv")
df$Sample_Year <- as.character(as.numeric(df$Sample_Year))
df$Zone <- as.character(as.numeric(df$Zone))

slope <- subset(df, Habitat == "slope")
flat <- subset(df, Habitat == "flat")

#create dataframe for ANOVA results
Species_Code <- unique(df$Species_Code)
Anova_df <- as.data.frame(Species_Code)
Anova_df$Habitat <- NA
Anova_df$Sample_Year <- NA
Anova_df$Zone <- NA
Anova_df$Substrate <- NA
Anova_df$Depth_slope <- NA
Anova_df$Inhabitated <- NA
Anova_df$Sampling_Time <- NA
Anova_df$Island <- NA
Anova_df$Distance_Makassar_km <- NA

Anova_df_flat <- Anova_df[,-c(2, 6)] 
Anova_df_slope <- Anova_df[,-2] 

Column <- c("Habitat", "Sample_Year", "Zone", "Substrate", "Depth_slope", "Inhabitated", "Sampling_Time", "Island", "Distance_Makassar_km",
            "Species_Code", "RE", "Sample_Code")
df0 <- df[, Column]
x <- 0

#with full dataset
for (i in 2:10) {
  
  x <- x + 1
  df_mean <- df0 %>% group_by(df0[Column[[x]]], Species_Code, Sample_Code) %>% summarise(Mean = mean(RE))
  colnames(df_mean) <- c("variable", "Species_Code", "Sample_Code", "Mean")
  df_mean_spread <- dcast(df_mean, variable + Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "Mean")
  
  Anova_df[1,i] <- anova(lm(Amphisorus_SpL ~ variable, data = df_mean_spread))[1,5]
  Anova_df[2,i] <- anova(lm(Amphisorus_SpS ~ variable, data = df_mean_spread))[1,5]
  Anova_df[3,i] <- anova(lm(Amphistegina_lessonii ~ variable, data = df_mean_spread))[1,5]
  Anova_df[4,i] <- anova(lm(Amphistegina_lobifera ~ variable, data = df_mean_spread))[1,5]
  Anova_df[5,i] <- anova(lm(Amphistegina_radiata ~ variable, data = df_mean_spread))[1,5]
  Anova_df[6,i] <- anova(lm(Baculogypsinoides_spinosus ~ variable, data = df_mean_spread))[1,5]
  Anova_df[7,i] <- anova(lm(Calcarina_mayori ~ variable, data = df_mean_spread))[1,5]
  Anova_df[8,i] <- anova(lm(Calcarina_hispida ~ variable, data = df_mean_spread))[1,5]
  Anova_df[9,i] <- anova(lm(Calcarina_spengleri ~ variable, data = df_mean_spread))[1,5]
  Anova_df[10,i] <- anova(lm(Elphidium_spp ~ variable, data = df_mean_spread))[1,5]
  Anova_df[11,i] <- anova(lm(Heterostegina_depressa ~ variable, data = df_mean_spread))[1,5]
  Anova_df[12,i] <- anova(lm(Neorotalia_calcar ~ variable, data = df_mean_spread))[1,5]
  Anova_df[13,i] <- anova(lm(Neorotalia_gaimardi ~ variable, data = df_mean_spread))[1,5]
  Anova_df[14,i] <- anova(lm(Operculina_ammonoides ~ variable, data = df_mean_spread))[1,5]
  Anova_df[15,i] <- anova(lm(Peneroplis_planatus ~ variable, data = df_mean_spread))[1,5]
  Anova_df[16,i] <- anova(lm(Sorites_orbiculus ~ variable, data = df_mean_spread))[1,5]
  Anova_df[17,i] <- anova(lm(Sphaerogypsina_globulus ~ variable, data = df_mean_spread))[1,5]
  
} 


#for reef slope
Column1 <- c("Sample_Year", "Zone", "Substrate", "Depth_slope", "Inhabitated", "Sampling_Time", "Island", "Distance_Makassar_km",
            "Species_Code", "RE", "Sample_Code")
df1 <- slope[, Column1]
x <- 0

for (i in 2:9) {
  
  x <- x + 1

  df_mean <- df1 %>% group_by(df1[Column1[[x]]], Species_Code, Sample_Code) %>% summarise(Mean = mean(RE))
  colnames(df_mean) <- c("variable", "Species_Code", "Sample_Code", "Mean")
  df_mean_spread <- dcast(df_mean, variable + Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "Mean")
  
  Anova_df_slope[1,i] <- anova(lm(Amphisorus_SpL ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[2,i] <- anova(lm(Amphisorus_SpS ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[3,i] <- anova(lm(Amphistegina_lessonii ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[4,i] <- anova(lm(Amphistegina_lobifera ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[5,i] <- anova(lm(Amphistegina_radiata ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[6,i] <- anova(lm(Baculogypsinoides_spinosus ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[7,i] <- anova(lm(Calcarina_mayori ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[8,i] <- anova(lm(Calcarina_hispida ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[9,i] <- anova(lm(Calcarina_spengleri ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[10,i] <- anova(lm(Elphidium_spp ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[11,i] <- anova(lm(Heterostegina_depressa ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[12,i] <- anova(lm(Neorotalia_calcar ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[13,i] <- anova(lm(Neorotalia_gaimardi ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[14,i] <- anova(lm(Operculina_ammonoides ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[15,i] <- anova(lm(Peneroplis_planatus ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[16,i] <- anova(lm(Sorites_orbiculus ~ variable, data = df_mean_spread))[1,5]
  Anova_df_slope[17,i] <- anova(lm(Sphaerogypsina_globulus ~ variable, data = df_mean_spread))[1,5]
  
} 

#for reef flat
Column2 <- c("Sample_Year", "Zone", "Substrate", "Inhabitated", "Sampling_Time", "Island", "Distance_Makassar_km",
             "Species_Code", "RE", "Sample_Code")
df2 <- flat[, Column2]
x <- 0

for (i in 2:8) {
  
  x <- x + 1
  
  df_mean <- df2 %>% group_by(df2[Column2[[x]]], Species_Code, Sample_Code) %>% summarise(Mean = mean(RE))
  colnames(df_mean) <- c("variable", "Species_Code", "Sample_Code", "Mean")
  df_mean_spread <- dcast(df_mean, variable + Sample_Code ~ Species_Code, fun.aggregate = sum, value.var = "Mean")
  
  Anova_df_flat[1,i] <- anova(lm(Amphisorus_SpL ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[2,i] <- anova(lm(Amphisorus_SpS ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[3,i] <- anova(lm(Amphistegina_lessonii ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[4,i] <- anova(lm(Amphistegina_lobifera ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[5,i] <- anova(lm(Amphistegina_radiata ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[6,i] <- anova(lm(Baculogypsinoides_spinosus ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[7,i] <- anova(lm(Calcarina_mayori ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[8,i] <- anova(lm(Calcarina_hispida ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[9,i] <- anova(lm(Calcarina_spengleri ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[10,i] <- anova(lm(Elphidium_spp ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[11,i] <- anova(lm(Heterostegina_depressa ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[12,i] <- anova(lm(Neorotalia_calcar ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[13,i] <- anova(lm(Neorotalia_gaimardi ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[14,i] <- anova(lm(Operculina_ammonoides ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[15,i] <- anova(lm(Peneroplis_planatus ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[16,i] <- anova(lm(Sorites_orbiculus ~ variable, data = df_mean_spread))[1,5]
  Anova_df_flat[17,i] <- anova(lm(Sphaerogypsina_globulus ~ variable, data = df_mean_spread))[1,5]
  
} 

#round pvalues to 4 digits
Anova_df_flat <- cbind(Anova_df_flat[,1], round(Anova_df_flat[,2:8], digits = 4))
Anova_df_slope <- cbind(Anova_df_slope[,1], round(Anova_df_slope[,2:9], digits = 4))
Anova_df <- cbind(Anova_df[,1], round(Anova_df[,2:10], digits = 4))

write.csv(Anova_df_flat, "~/Downloads/Dataset/20210301-anova-results-flat.csv")
write.csv(Anova_df_slope, "~/Downloads/Dataset/20210301-anova-results-slope.csv")
write.csv(Anova_df, "~/Downloads/Dataset/20210301-anova-results-bothhabitats.csv")





