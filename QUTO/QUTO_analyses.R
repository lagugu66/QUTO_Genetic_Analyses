#####################
#     Libraries     #
#####################

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)
library(diveRsity)

#####################
#     Analyses      #
#####################
#load in genepop file as a genind object 
QUTO_genind <- read.genepop("C:/Users/LAguiniga/Documents/QUTO/QUTO_allpop.gen",
                            ncode = 3)

#allele frequency category lists 
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare",
                   "reg_rare","loc_com_d1","loc_com_d2","loc_rare")

#reorganize genind file 
QUTO_garden_genind <- seppop(QUTO_genind)[[1]]
levels(QUTO_garden_genind@pop) <- "Garden"

QUTO_wild_genind <- repool(seppop(QUTO_genind)[2:19])
levels(QUTO_wild_genind@pop) <- rep("Wild",18)

#smash back together
QUTO_garden_wild_genind <- repool(QUTO_garden_genind,
                                  QUTO_wild_genind)

##start genetic analyses
#create genetic summary of the genind file 
sp_sum <- summary(QUTO_garden_wild_genind)
#create poppr file 
sp_poppr <- poppr(QUTO_garden_wild_genind)
#save mean for final output table 
sp_hexp_mean <- sp_poppr[1:length(levels(QUTO_garden_wild_genind@pop)),10]
#allele numbers by pop 
sp_nall <- sp_sum$pop.n.all
#individual numbers
sp_ind <- sp_poppr[1:length(levels(QUTO_garden_wild_genind@pop)), 2:3]
#save allelic richness for comparison
sp_allrich_list <- allel.rich(QUTO_garden_wild_genind)$all.richness
sp_allrich_mean <- colMeans(allel.rich(QUTO_garden_wild_genind)$all.richness)	

#create data frame 
sp_allpop_gendiv_sumstat_df <- signif(cbind(sp_ind, sp_nall, sp_allrich_mean, sp_hexp_mean),3)
#create data frame 
sp_allpop_gendiv_sumstat_df <- signif(cbind(sp_ind, sp_nall, sp_allrich_mean, sp_hexp_mean),3)

#name rows 
rownames(sp_allpop_gendiv_sumstat_df) <- levels(QUTO_garden_wild_genind@pop)
colnames(sp_allpop_gendiv_sumstat_df) <- c("Ind","MLG", "NAll", "All_Rich", "Hexp")

#write out data frame
write.csv(sp_allpop_gendiv_sumstat_df, "C:/Users/LAguiniga/Documents/QUTO/Analyses/Results/Sum_Stats/QUTO_gendiv_sumstats_df.csv")

sp_allrich_df <- gather(as.data.frame(sp_allrich_list))
sp_allrich_p_value <- kruskal.test(sp_allrich_df[,2]~sp_allrich_df[,1])[3]

#make a boxplot
boxplot(sp_allrich_df[,2]~sp_allrich_df[,1], ylim = c(0,18), col = c("green","purple"))

#hexp
sp_hexp_df <- rbind(summary(QUTO_garden_genind)[[7]], summary(QUTO_wild_genind)[[7]])
sp_hexp_df <- as.data.frame(sp_hexp_df)
rownames(sp_hexp_df) <- c("Garden", "Wild")
sp_hexp_df <- t(sp_hexp_df)
sp_hexp_df <- gather(as.data.frame(sp_hexp_df))


sp_hexp_p_value <- (sp_hexp_df[,2]~sp_hexp_df[,1])[3]
sp_hexp_p_value <- kruskal.test(sp_hexp_df[,2]~sp_hexp_df[,1])[3]
boxplot(sp_hexp_df[,2]~sp_hexp_df[,1], ylim = c(0,1), col = c("blue","pink"))
#########
#     Allelic Representation
#convert the wild genind object to a genpop object
sp_wild_genpop <- genind2genpop(seppop(QUTO_garden_wild_genind)[2]$Wild)

#create documents for comparison 
n_ind_W <- nrow(QUTO_wild_genind@tab);  n_ind_G <- nrow(QUTO_garden_genind@tab); 
sp_alleles_cap <- colSums(seppop(QUTO_garden_wild_genind)[[1]]@tab,na.rm=T)

#first calculate the frequency categories of alleles in the wild individuals   	
sp_allele_cat <- get.allele.cat(sp_wild_genpop, 1, 1, n_ind_W, n_drop = 0, glob_only = TRUE)	

##
sp_all_exist_df <- matrix(nrow = 1, ncol = 9)
sp_wild_cap_df <- matrix(nrow = 1, ncol = 9)
sp_allele_cap <- matrix(nrow = 1, ncol = 9)
#calculating alleles that exist by allelic category
for(a in 1:length(list_allele_cat)) sp_all_exist_df[,a] <- round(sum(sp_alleles_cap[sp_allele_cat[[a]]] > 0))
    
#now determine how many wild alleles were captured per category 
for(a in 1:length(list_allele_cat)) sp_wild_cap_df[,a] <- round(sum(sp_alleles_cap[sp_allele_cat[[a]]] > 0)/length(sp_allele_cat[[a]]),4)
   
#code to store as one data frame 
for(a in 1:length(list_allele_cat))sp_allele_cap[,a] <- paste0(signif((sp_wild_cap_df[,a]*100),3), "% (", signif(sp_all_exist_df[,a],3), ")")

colnames(sp_allele_cap) <- list_allele_cat
write.csv(sp_allpop_gendiv_sumstat_df, "C:/Users/LAguiniga/Documents/QUTO/Analyses/Results/Sum_Stats/QUTO_sp_allele_cap_df.csv")