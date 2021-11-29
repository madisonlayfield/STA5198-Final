#########################
## TWAS-FUSION RESULTS ##
#########################

#Run TWAS-Fusion in terminal

#Import results
data_files <- list.files("/Users/MacBook/Desktop/t1d_results/", pattern=".dat") 
for(i in 1:length(data_files)) {                             
  assign(paste0("t1d","-", i),                                 
         read.csv(paste0("/Users/MacBook/Desktop/t1d_results/",
                         data_files[i]), sep=""))
}
t1d <- do.call(rbind, mget(ls(pattern = "t1d")))

#Find Significant using FDR correction
t1d$p_FDR <- p.adjust(t1d$TWAS.P, method="fdr")
fusion_sig <- t1d[t1d$p_FDR < 0.05 & complete.cases(t1d$p_FDR),]


###################
## FOCUS RESULTS ##
###################
library(data.table)

#Add chr, bp, beta, se, and p-values to LDSC sumstats
sumstats <- read.delim("/gpfs/home/mil18b/T1D/t1d.sumstats", sep=" ")
data <- read.delim("/gpfs/home/mil18b/T1D/34012112-GCST90014023-EFO_0001359-Build38.f.tsv.gz", sep=" ")
sumstats <- merge(sumstats, data, by.x = c("SNP"), by.y = c("SNP"), all.y=FALSE)
sumstats <- sumstats[,-c(9,10,13)]
colnames(sumstats) <- c("SNP", "Z", "A1", "A2", "N", "P", "CHR", "BP", "beta", "se")
write.table(sumstats, "t1d_sumstats.txt", sep=" ", row.names=FALSE, quote=FALSE)

#Run FOCUS in terminal

#Import results
data_files <- list.files("/Users/MacBook/Desktop/focus_results/", pattern=".tsv") 
for(i in 1:length(data_files)) {                             
  assign(paste0("t1d.chr", i),                                 
         read.csv(paste0("/Users/MacBook/Desktop/focus_results/",
                         data_files[i]), sep=""))
}
focus_t1d <- do.call(rbind, mget(ls(pattern = "t1d.chr")))

#Get genes of interest
focus_t1d <- focus_t1d[focus_t1d$mol_name %in% fusion_sig$ID,]

#Get likely causal genes (in credible set)
focus_sig <- focus_t1d[focus_t1d$in_cred_set == 1,]
