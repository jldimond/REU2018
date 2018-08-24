##################################################################################

setwd("~/Documents/Anthopleura_genome")

#open nanopore methylation data
Aele_meth_freq <-read.delim("methylation_frequency_thresh10.tsv", sep = "\t", header=T)

#get average meth frequency
Aele_meth_avg <- aggregate(methylated_frequency ~ chromosome, Aele_meth_freq, mean)

#get sum of CpG motifs
Aele_CpG_motifs <- aggregate(num_motifs_in_group ~ chromosome, Aele_meth_freq, sum)

#merge meth frequency and sum motifs
Meth_avg_sum_motifs <-merge(Aele_meth_avg,Aele_CpG_motifs,by="chromosome")

#change directory 
setwd("~/Documents/Coral-CpG-master/analyses/Aele")

#open file with seq length and CpG counts (length is V4, CpG count is V5)
#with CpG counts from the transcriptome sequences, we can estimate coverage of the nanopore data
Aele_comb <-read.delim("comb", header=F)

#just get seq length and CpG counts
Aele_comb2 <- Aele_comb[,c(1,4,5)]

#open CpG O/E data
Aele_cpgoe <-read.delim("ID_CpG", header=F)

#remove extra space in seq id
Aele_cpgoe2 <- as.data.frame(apply(Aele_cpgoe,2,function(x)gsub('\\s+', '',x)))

#change CpG data from factor to numeric
Aele_cpgoe2[,2] <-as.numeric(as.character(Aele_cpgoe2[,2]))

#merge nanopore data with seq length data
temp_merge <-merge(Meth_avg_sum_motifs,Aele_comb2,by.x="chromosome",by.y="V1")

#merge above data with CpG O/E
All_data <-merge(temp_merge,Aele_cpgoe2,by.x="chromosome",by.y="V1")

#rename columns
colnames(All_data)[4] <- "seq_length"
colnames(All_data)[5] <- "num_CpGs"
colnames(All_data)[6] <- "CpG_OE"

#now determine CpG coverage of nanopore data
#this divides total number of motifs (just another way of saying CpGs) from nanopore data
#by total number of CpGs in the sequence 
cpg_cov <- All_data[,3]/All_data[,5]

#bind cpg_cov to dataframe
All_data <- cbind(All_data, cpg_cov)

#now set threshold of 20% CpG coverage to data
All_data_thresh <- All_data[All_data[,7] >= 0.2,]

#now set threshold of 2 for CpG O/E 
All_data_thresh <- All_data_thresh[All_data_thresh[,6] <= 2,]

#smooth scatter plot of nanopore meth data vs. CpG O/E data
plot(All_data_thresh[,2], All_data_thresh[,6])

#linear regression of nanopore meth data vs. CpG O/E data
model <- lm(All_data_thresh[,2] ~ All_data_thresh[,6])
abline(model, col = "green")

# load ggpolot and ggpubr packages

library(ggplot2)
library(ggpubr)

# Scatter plot with density panels using ggplot and ggpubr
sp <- ggscatter(All_data_thresh, x = "methylated_frequency", y = "CpG_OE",
                color = "purple", size = 0.4, alpha = 0.1,
                add = "loess", add.params = list(color = "blue"),
                cor.coef = TRUE, 
                cor.coeff.args = list(method = "pearson", 
                                      label.x.npc = "right", label.y.npc = "top"),
                ylab = "CpG O/E", xlab = "Methylated frequency")+
  border()                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(All_data_thresh, "methylated_frequency")
yplot <- ggdensity(All_data_thresh, "CpG_OE")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE) 


