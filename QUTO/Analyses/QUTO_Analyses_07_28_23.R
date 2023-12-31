##############
# libraries #
##############

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)
library(tidyverse)
############################
#     Load Data Files      #
############################
#load in genepop file as a genind object 
#change working directory 
setwd("C:/Users/laguiniga/Documents/GitHub/QUTO_Genetic_Analyses")

#load in data file
QUTO_genind <- read.genepop("QUTO/Data_Files/QUTO_byisland.gen",
                            ncode = 3)
#allele frequency category lists 
all_cat_list <- c("global","glob_v_com","glob_com","glob_lowfr","glob_rare")

#
source("C:/Users/laguiniga/Documents/GitHub/QUTO_Genetic_Analyses/QUTO/Data_Files/Fa_sample_funcs.R")

#
dup_reps <- c(0:9)

#####################
#     Analyses      #
#####################
##reorganize genind file to have garden/wild pops
#separate garden individuals and rename
QUTO_garden_genind <- seppop(QUTO_genind)[[1]]
levels(QUTO_garden_genind@pop) <- "Garden"

#separate into wild pops 
QUTO_wild_genind <- repool(seppop(QUTO_genind)[2:5])

#repool into garden/wild genind
QUTO_garden_wild_genind <- repool(QUTO_garden_genind,
                                  QUTO_wild_genind)

##start genetic analyses
#create genetic summary of the genind file 
QUTO_garden_sum <- summary(seppop(QUTO_garden_wild_genind)[[1]])
#save mean for final output table 
QUTO_garden_hexp_mean <- mean(QUTO_garden_sum[[7]])
#allele numbers by pop 
QUTO_garden_nall <- QUTO_garden_sum$pop.n.all
#individual numbers
QUTO_garden_ind <- QUTO_garden_sum[[1]]
#save allelic richness for comparison
QUTO_garden_allrich_list <- allel.rich(seppop(QUTO_garden_wild_genind)[[1]])$all.richness
QUTO_garden_allrich_mean <- colMeans(QUTO_garden_allrich_list)	

##calculate summary statistics for wild populations 
#create a data frame to store results in 
QUTO_wild_sumstats <- matrix(nrow = length(levels(QUTO_wild_genind@pop)),
                             ncol = 4)

#loop to calculate sum stats by wild populations 
for(wp in 1:length(levels(QUTO_wild_genind@pop))){
  
  
  QUTO_wild_sumstats[wp,1] <- summary(seppop(QUTO_wild_genind)[[wp]])[[1]]
  QUTO_wild_sumstats[wp,2] <- summary(seppop(QUTO_wild_genind)[[wp]])[[4]]
  QUTO_wild_sumstats[wp,3] <- mean(allel.rich(seppop(QUTO_garden_wild_genind)[[wp]])$all.richness)
  QUTO_wild_sumstats[wp,4] <- mean(summary(seppop(QUTO_wild_genind)[[wp]])[[7]])
  
}

#name rows and columns 
rownames(QUTO_wild_sumstats) <- c("SR","C","SCr","SCl")
colnames(QUTO_wild_sumstats) <- c("N","NAll","AllRich","HExp")

#create garden data frame 
QUTO_garden_sumstats <- signif(cbind(QUTO_garden_ind, QUTO_garden_nall, QUTO_garden_allrich_mean, QUTO_garden_hexp_mean),3)
colnames(QUTO_garden_sumstats) <- c("N","NAll","AllRich","HExp")

#create dataframe
QUTO_garden_wild_sumstats <- rbind(QUTO_garden_sumstats, QUTO_wild_sumstats)

#write out data frame
write.csv(QUTO_garden_wild_sumstats, "C:/Users/laguiniga/Documents/GitHub/QUTO_Genetic_Analyses/QUTO/Analyses/Results/Sum_Stats/QUTO_sumstats_df.csv")


###################################
#     Ex situ representation      #
###################################
##reorganize data files for representation analyses
levels(QUTO_wild_genind@pop) <- rep("Wild", 4)

#now reorganize genind object
QUTO_garden_wild_onepop_genind <- repool(QUTO_garden_genind, QUTO_wild_genind)

##create allelic richness dataframe
#allelic richness data frame for statistical test
QUTO_allrich_df <- gather(as.data.frame(allel.rich(QUTO_garden_wild_onepop_genind)[[1]]))

QUTO_hexp_df <- cbind(summary(QUTO_garden_genind)[[7]],summary(QUTO_wild_genind)[[7]])
colnames(QUTO_hexp_df) <- c("Garden", "Wild")

QUTO_hexp_g_df <- gather(as.data.frame(QUTO_hexp_df))

#create a data frame to store means and pvalue
QUTO_exsitu_gendiv_comp_df <- matrix(nrow = 3, ncol = 2)

#save results
QUTO_exsitu_gendiv_comp_df[1,1] <- mean(QUTO_allrich_df[QUTO_allrich_df$key == "Garden",]$value)
QUTO_exsitu_gendiv_comp_df[2,1] <- mean(QUTO_allrich_df[QUTO_allrich_df$key == "Wild",]$value)
QUTO_exsitu_gendiv_comp_df[3,1] <- kruskal.test(QUTO_allrich_df[,2]~QUTO_allrich_df[,1])[[3]]

##hexp
QUTO_exsitu_gendiv_comp_df[1,2] <- mean(QUTO_hexp_g_df[QUTO_hexp_g_df$key == "Garden",]$value)
QUTO_exsitu_gendiv_comp_df[2,2] <- mean(QUTO_hexp_g_df[QUTO_hexp_g_df$key == "Wild",]$value)
QUTO_exsitu_gendiv_comp_df[3,2] <- kruskal.test(QUTO_hexp_g_df[,2]~QUTO_hexp_g_df[,1])[[3]]

rownames(QUTO_exsitu_gendiv_comp_df) <- c("Garden", "Wild", "p-value")
colnames(QUTO_exsitu_gendiv_comp_df) <- c("AllRich", "HExp")

##write out data frame 
write.csv(QUTO_exsitu_gendiv_comp_df, "C:/Users/laguiniga/Documents/GitHub/QUTO_Genetic_Analyses/QUTO/Analyses/Results/Garden_Wild_Comparison/QUTO_exsitu_gendiv_comp_df.csv")

colnames(QUTO_allrich_df) <- c("Population_Type", "Allelic_Richness")

###make a boxplot of the garden/wild comparisons
ggplot(QUTO_allrich_df, aes(x = Population_Type, y = Allelic_Richness, fill = Population_Type)) + 
  geom_boxplot() + scale_color_manual(values = c("lightblue", "hotpink"))

###hexp_Garden_Wild_Comp
##hexp
QUTO_exsitu_gendiv_comp_df[1,2] <- mean(QUTO_hexp_g_df[QUTO_hexp_g_df$key == "Garden",]$value)
QUTO_exsitu_gendiv_comp_df[2,2] <- mean(QUTO_hexp_g_df[QUTO_hexp_g_df$key == "Wild",]$value)
QUTO_exsitu_gendiv_comp_df[3,2] <- kruskal.test(QUTO_hexp_g_df[,2]~QUTO_hexp_g_df[,1])[[3]]

rownames(QUTO_exsitu_gendiv_comp_df) <- c("Garden", "Wild", "p-value")
colnames(QUTO_exsitu_gendiv_comp_df) <- c("AllRich", "HExp")

##write out data frame 
write.csv(QUTO_exsitu_gendiv_comp_df, "C:/Users/laguiniga/Documents/GitHub/QUTO_Genetic_Analyses/QUTO/Analyses/Results/Garden_Wild_Comparison/QUTO_exsitu_gendiv_comp_df.csv")

colnames(QUTO_hexp_df) <- c("Garden", "Wild")
QUTO_hexp_t_df <- gather(as.data.frame(QUTO_hexp_df))
colnames(QUTO_hexp_t_df) <- c("Population_Type", "Heterozygosity")

###ggplot_hexp_comp
ggplot(QUTO_hexp_t_df, aes(x = Population_Type, y = Heterozygosity, fill = Population_Type)) + 
  geom_boxplot() + scale_color_manual(values = c("lightblue", "green"))

###################################
#     Allelic Representation      #
###################################
#convert the wild genind object to a genpop object
QUTO_wild_genpop <- genind2genpop(seppop(QUTO_garden_wild_onepop_genind)[2]$Wild)

#create documents for comparison 
n_ind_W <- nrow(seppop(QUTO_garden_wild_onepop_genind)[[2]]@tab);  n_ind_G <- nrow(seppop(QUTO_garden_wild_onepop_genind)[[1]]@tab); 
QUTO_all_rep <- colSums(seppop(QUTO_garden_wild_onepop_genind)[[1]]@tab,na.rm=T)

#first calculate the frequency categories of alleles in the wild individuals   	
QUTO_all_cat <- get.allele.cat(QUTO_wild_genpop, 1, 1, n_ind_W, n_drop = 0, glob_only = TRUE)	

#remove regional alleles from QUTO all cat
QUTO_all_cat <- QUTO_all_cat[1:5]

#alleles existing 
QUTO_all_exist_df <- matrix(nrow = length(dup_reps), ncol = length(all_cat_list))
QUTO_garden_rep_df <- matrix(nrow = length(dup_reps), ncol = length(all_cat_list))
QUTO_garden_rep_per_df <- matrix(nrow = length(dup_reps), ncol = length(all_cat_list))

#loop to calculate allelic representation 
for(dup in dup_reps){
  for(cat in 1:length(QUTO_all_cat)){
    
    #first save all of the existing wild alleles
    QUTO_all_exist_df[dup+1,cat] <- round(sum(QUTO_all_rep[QUTO_all_cat[[cat]]] > dup))
    
    #now count the alleles represented in each category ex situ
    QUTO_garden_rep_df[dup+1,cat] <- round(sum(QUTO_all_rep[QUTO_all_cat[[cat]]] > dup)/length(QUTO_all_cat[[cat]]),4)
    
    #save in a presentation data frame 
    QUTO_garden_rep_df[dup+1,cat] <- paste0(signif((QUTO_garden_rep_df[dup+1,cat]*100),3), "% (", signif(QUTO_all_exist_df[dup+1,cat],3), ")")
    
  }
}


#alleles existing 
QUTO_all_exist_df <- matrix(nrow = length(dup_reps), ncol = length(all_cat_list))
QUTO_wild_cap_df <- matrix(nrow = length(dup_reps), ncol = length(all_cat_list))
QUTO_all_rep <- matrix(nrow = length(dup_reps), ncol = length(all_cat_list))
#calculating alleles that exist by allelic category
for(a in 1:length(list_allele_cat)) QUTO_all_exist_df[,a] <- round(sum(QUTO_all_rep[QUTO_all_cat[[a]]] > 0))

#now determine how many wild alleles were captured per category 
for(a in 1:length(list_allele_cat)) QUTO_wild_cap_df[,a] <- round(sum(QUTO_all_rep[QUTO_all_cat[[a]]] > 0)/length(sp_allele_cat[[a]]),4)

#code to store as one data frame 
for(a in 1:length(list_allele_cat))QUTO_all_rep[,a] <- paste0(signif((QUTO_wild_cap_df[,a]*100),3), "% (", signif(sp_all_exist_df[,a],3), ")")

colnames(QUTO_all_rep) <- list_allele_cat
write.csv(sp_allpop_gendiv_sumstat_df, "C:/Users/LAguiniga/Documents/QUTO/Analyses/Results/Sum_Stats/QUTO_sp_allele_cap_df.csv")