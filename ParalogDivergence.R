## 
##
## Script name: MBE_allfigures
##
## Purpose of script: 1) Generating plots from all analysis
##                    2) Statistical testing  
##
## Author: Sunil Kenchanmane Raju
## contact: kenchanmane@gmail.com
##
## Notes: Script combined from earlier versions (v2021)
##   
## _______________________________________________

setwd("C:/Users/kench/Dropbox/ANALYSIS/MBE_V2") #Laptop
setwd("C:/Users/kenchanm/Dropbox/ANALYSIS/MBE_V2")


df1 <- read.csv("mutation_single_gene_level.csv") # Mutation data from Grey et al 2021 
df2 <- read.csv("../GeneDuplication_V3/data/Athaliana/methylpy/results/Athaliana_MetClassified_genes_V6.csv")[,1:2] # methylation data from Raju et al. 2023


single_gene_level <- merge (df1, df2, by.x='gene', by.y="Feature")[,-2]  # -2 to remove the row numbers from df2

#write.csv(single_gene_level, file="mutation_single_gene_level_v6.csv") # v6 is after TE filtering

library(ggplot2)
library(gghalves)
library(RColorBrewer)

single_gene_level$order <- ifelse(single_gene_level$Methylation=="gbM",2,NA)
single_gene_level$order <- ifelse(single_gene_level$Methylation=="teM",3,single_gene_level$order)
single_gene_level$order <- ifelse(single_gene_level$Methylation=="unM",1,single_gene_level$order)
#single_gene_level$order <- ifelse(single_gene_level$Methylation=="Unclassified",4,single_gene_level$order)

single_gene_level <- single_gene_level[single_gene_level$Methylation != "Unclassified",]
single_gene_level$Methylation <- factor (single_gene_level$Methylation, levels = c("unM", "gbM", "teM"))

single_gene_level <- single_gene_level[single_gene_level$Duplication != "unclassified",]
single_gene_level <- single_gene_level[single_gene_level$Duplication != "singletons",]

single_gene_level$Duplication <- factor (single_gene_level$Duplication, levels = c("wgd", "tandem", "proximal", "transposed", "dispersed"))

#### Figures ####

# Fig 1B
# Genic methylation classification and polymorphisms in Arabidopsis accessions
# Figure is updated in Inkscape to change y-axis into percentages 

# violin plot

p <- ggplot(single_gene_level[!is.na(single_gene_level$Methylation),]) +
  theme_bw() +
  labs(title = "", subtitle = "")+
  geom_violin(aes(x=reorder(Methylation,order),y=NSNPs,fill=Methylation))+
  scale_fill_manual(values = c("unM"="Cornflowerblue", "gbM" = "#00b067", "teM"="#FC8D62"))+
  theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=12)) +
  theme(axis.title=element_text(size=14,face="bold")) +
  labs(x="",y="Polymorphism in A. thaliana accessions")
  
#halfplot
              
p <- ggplot(single_gene_level[!is.na(single_gene_level$Methylation),]) +
  theme_bw() +
  labs(title = "", subtitle = "")+
  geom_half_boxplot(aes(x=reorder(Methylation,order),y=NSNPs,fill=Methylation),
                    side="l",errorbar.draw=TRUE,outlier.color=NA, outlier.shape = NA) +
  scale_fill_manual(values = c("unM"="Cornflowerblue", "gbM" = "#00b067", "teM"="#FC8D62")) +
  coord_cartesian(ylim = quantile(single_gene_level$NSNPs, c(0.001, 0.9999), na.rm=TRUE)) +
  geom_half_point(aes(x=reorder(Methylation,order),y=NSNPs,color=Methylation),
                  side="r",size=0.25) +
  scale_color_manual(values = c("unM"="Cornflowerblue", "gbM" = "#00b067", "teM"="#FC8D62")) +
  theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=12)) +
  theme(axis.title=element_text(size=14,face="bold")) +
  labs(x="",y="Polymorphism in A. thaliana accessions")
p

# Fig 1C, 1D and S1

#Correlation between NSNPs and frequency of methylation classification in the natural population

library("ggpubr")
single_gene_levela <- single_gene_level[single_gene_level$Methylation != "Unclassified",]
single_gene_levelb <- single_gene_levela[!is.na(single_gene_levela$Methylation),]

# gbM - Fig 1C
ggscatter(single_gene_levelb, x ="p_gbM" , y = "NSNPs",
          cor.coeff.args = list(method = "Spearman", label.x.npc = 0.05, label.y.npc = 0.95),
          add = "reg.line", conf.int = TRUE, size=0.5, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Frequency of gbM in 928 A. thaliana accessions", ylab = "SNP polymorphism")  
#Change y each time

# teM - Fig 1D
ggscatter(single_gene_levelb, x ="p_teM" , y = "NSNPs",
          cor.coeff.args = list(method = "Spearman", label.x.npc = 0.05, label.y.npc = 0.95),
          add = "reg.line", conf.int = TRUE, size=0.5, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Frequency of teM in 928 A. thaliana accessions", ylab = "SNP polymorphism")  
#Change y each time

# unM - Fig S1
ggscatter(single_gene_levelb, x ="p_unM" , y = "NSNPs",
          cor.coeff.args = list(method = "Spearman", label.x.npc = 0.05, label.y.npc = 0.95),
          add = "reg.line", conf.int = TRUE, size=0.5, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Frequency of teM in 928 A. thaliana accessions", ylab = "SNP polymorphism")  
#Change y each time

# Histone supplement figure 
# Change of the name of the histone modification each time
# might have to change ylim to play with the outliers

p <- ggplot(single_gene_level[!is.na(single_gene_level$Methylation),]) +
  theme_bw() +
  #labs(title = "Genic methylation - H3K4me1", subtitle = "")+
  geom_violin(aes(x=reorder(Methylation,order),y=H3K56ac, fill=Methylation),
              outlier.color=NA, size =0.5) +
  theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=12)) +
  theme(axis.text.y=element_text(color="black",hjust=1,size=12)) +
  labs(x="",y="H3K56ac")
p + scale_fill_manual(values = c("unM"="Cornflowerblue", "gbM" = "#00b067", "teM"="#FC8D62")) + ylim(-1.5,7)

#Figure 2A
#GenicMethylation-TajimasD

p <- ggplot(single_gene_level[!is.na(single_gene_level$Methylation),]) +
  theme_bw() +
  #labs(title = "Genic methylation - Tajima's D", subtitle = "")+
  scale_fill_manual(values = c("unM"="Cornflowerblue", "gbM" = "#00b067", "teM"="#FC8D62")) +
  geom_half_boxplot(aes(x=reorder(Methylation,order),y=TajimaD,fill=Methylation),
                    side="l",errorbar.draw=TRUE,outlier.color=NA) +
  geom_half_point(aes(x=reorder(Methylation,order), y=TajimaD,color=Methylation),
                  side="r",size=0.5, alpha=0.2) +
  scale_color_manual(values = c("unM"="Cornflowerblue", "gbM" = "#00b067", "teM"="#FC8D62")) +
  theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=12)) +
  theme(axis.text.y=element_text(color="black",hjust=1,size=12)) +
 
  labs(x="",y="Tajima's D")
p

#Figure S2
#GenicMethylation-DnDs

p <- ggplot(single_gene_level[!is.na(single_gene_level$Methylation),]) +
  theme_bw() +
  #labs(title = "Genic methylation - Tajima's D", subtitle = "")+
  scale_fill_manual(values = c("unM"="Cornflowerblue", "gbM" = "#00b067", "teM"="#FC8D62")) +
  geom_boxplot(aes(x=reorder(Methylation,order),y=DnDs,fill=Methylation)) +
  theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=12)) +
  theme(axis.text.y=element_text(color="black",hjust=1,size=12)) +
  labs(x="",y="DnDs")
p + ylim(0, 3)

#### Figure 3A ####
# Duplicate pairs - NSNPs

df3 <- read.csv("../GeneDuplication_V3/figures_V6/Athaliana/Athaliana_Duplicate_pair_met.csv", header=TRUE)

df4 <- merge(df3, df1[, c(2,29,36,44,45,49:51)], by.x="Duplicate.1", by.y ="gene") #data for dup pairs
df5 <- merge(df4, df1[, c(2,29,36,44,45,49:51)], by.x="Duplicate.2", by.y ="gene")

df5a <- df5[df5$Duplicate.1_Methylation != "Unclassified",] #removing unclassified
df5a <- df5a[df5a$Duplicate.2_Methylation != "Unclassified",]

df5b <- df5a[c(2,1,3:21)]

#write.csv(df5b, "DuplicatePairs_Mutation_Data.csv")

data2 <- data.frame(
  name=c("unM-gbM", "unM-gbM","unM-teM","unM-teM","gbM-teM", "gbM-teM"),
  name2=c("Paralog1", "Paralog2"),
  value=c(6.872882859, 6.76726761, 8.141845168, 12.97429769, 8.760142288, 13.13597122), # values from gbm-tem, unm-tem, gbm-unm files
  CI=c(0.004263063, 0.004210303, 0.010626456, 0.017627378, 0.12621255, 0.144018837) # STD Error 
)

data2$name <- factor (data2$name, levels = c("unM-gbM", "unM-teM", "gbM-teM"))
library(ggplot2)
library(RColorBrewer)
p<- ggplot(data2, aes(x=name, y=value, fill=name2)) +
  theme_bw()+
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-CI, ymax=value+CI), width=.4, size=0.4, position=position_dodge(.9))+
  labs(x="",y="Polymorphism in A. thailana accessions")
p + scale_fill_brewer(palette = "Paired") + ylim(0,16)

#Duplicate pairs abs difference in NSNPS
library(dplyr)
df6 <- select(df5b, 1,2,7,11,18)
df6$abs <- abs(df6$NSNPs.x-df6$NSNPs.y)

df6$Classification <- factor(df6$Classification, levels = c("gbM-gbM", "unM-unM", "teM-teM", "gbM-unM", "gbM-teM", "unM-teM" ))
library(ggplot2)
library(gghalves)

p <- ggplot(df6) +
  theme_bw() +
  labs(title = "Transposed duplicates - NSNPs", subtitle = "")+
  geom_half_boxplot(aes(x=Classification,y=abs,fill=Classification),
                    side="l",errorbar.draw=TRUE,outlier.color=NA) +
  geom_half_point(aes(x=Classification,y=abs,color=Classification),
                  side="r",size=0.5) +
  theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=10)) +
  labs(x="",y="Polymorphism (Absolute Difference)")
p +ylim(0,25)


# Figure 3C
# Duplicate pairs - Tajima's D

data3 <- data.frame(
  name=c("unM-gbM", "unM-gbM","unM-teM","unM-teM","gbM-teM", "gbM-teM"),
  name2=c("Paralog1", "Paralog2"),
  value=c(-0.787067978, -0.817163335, -0.826661136, -1.044256103, -0.942665268, -0.995735717), # value calculated from mutations_single_gene_analysis.csv 
  CI=c(0.000694368, 0.000721975, 0.001661927, 0.001639029, 0.010428345, 0.013468702) # CI value obtained by prop.test (code below in statistical test section) 
)

data3$name <- factor (data3$name, levels = c("unM-gbM", "unM-teM", "gbM-teM"))
library(ggplot2)
library(RColorBrewer)
p<- ggplot(data3, aes(x=name, y=value, fill=name2)) +
  theme_bw()+
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-CI, ymax=value+CI), width=.4, size=0.4, position=position_dodge(.9)) 
p + scale_fill_brewer(palette = "Paired") + ylim(-1.3,0)


# Duplication - NSNPs

p <- ggplot(single_gene_level[!is.na(single_gene_level$Duplication),]) +
  theme_bw() +
  #labs(title = "Genic methylation - Tajima's D", subtitle = "")+
  scale_fill_manual(values = c("#003049", "#d62828", "#f77f00", "#fcbf49", "#eae2b7")) +
  geom_half_boxplot(aes(x=Duplication,y=NSNPs,fill=Duplication),
                    side="l",errorbar.draw=TRUE,outlier.color=NA, alpha=0.5) +
  #coord_cartesian(ylim = quantile(single_gene_level$TajimaD, c(0.001, 0.999))) +
  geom_half_point(aes(x=Duplication, y=NSNPs,color=Duplication),
                  side="r",size=0.5, alpha=0.1) +
  scale_color_manual(values = c("#003049", "#d62828", "#f77f00", "#fcbf49", "#eae2b7")) +
  theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=12)) +
  theme(axis.text.y=element_text(color="black",hjust=1,size=12)) +
  labs(x="",y="NSNPs")
p + ylim(0, 30)

# Duplication - Tajimas D

p <- ggplot(single_gene_level[!is.na(single_gene_level$Duplication),]) +
  theme_bw() +
  #labs(title = "Genic methylation - Tajima's D", subtitle = "")+
  scale_fill_manual(values = c("#003049", "#d62828", "#f77f00", "#fcbf49", "#eae2b7")) +
  geom_half_boxplot(aes(x=Duplication,y=TajimaD,fill=Duplication),
                    side="l",errorbar.draw=TRUE,outlier.color=NA, alpha=0.5) +
  #coord_cartesian(ylim = quantile(single_gene_level$TajimaD, c(0.001, 0.999))) +
  geom_half_point(aes(x=Duplication, y=TajimaD,color=Duplication),
                  side="r",size=0.5, alpha=0.1) +
  scale_color_manual(values = c("#003049", "#d62828", "#f77f00", "#fcbf49", "#eae2b7")) +
  theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=12)) +
  theme(axis.text.y=element_text(color="black",hjust=1,size=12)) +
  labs(x="",y="Tajima's D")
p + ylim(-2, 1.5)

# Duplication - Dn/Ds

p <- ggplot(single_gene_level[!is.na(single_gene_level$Duplication),]) +
  theme_bw() +
  #labs(title = "Genic methylation - Tajima's D", subtitle = "")+
  scale_fill_manual(values = c("#003049", "#d62828", "#f77f00", "#fcbf49", "#eae2b7")) +
  geom_boxplot(aes(x=Duplication,y=DnDs,fill=Duplication),
                    side="l",errorbar.draw=TRUE,outlier.color=NA, alpha=0.5)
  p + ylim(0, 3.5)

#### Transposed duplicates ####
  
  setwd("C:/Users/kench/Dropbox/ANALYSIS/GeneDuplication_V3")
  
  library(dplyr)
  
  df1 <- read.table("data/Athaliana/dupgen/results-unique/Athaliana.transposed.pairs-unique", header=TRUE,sep="\t")
  df2 <- df1[,c(1,3)]
  
  df3 <- read.csv("data/Athaliana/methylpy/results/Athaliana_MetClassified_genes_V6.csv",header=TRUE)
  df4 <- df3[,c(2,31)]
  
  df5 <- merge(df2, df4, by.x="Transposed", by.y="Feature")
  colnames(df5)[3] <- "Transposed_met" 
  df6 <- merge(df5, df4, by.x="Parental", by.y="Feature")
  colnames(df6)[4] <- "Parental_met"
  
  df6 <- select(df6, 2,1,3,4)
  
  #write.csv(df6, "figures_V6/Athaliana/Transposed_direction_methylation.csv")
  
  df7 <- read.csv("../MBE_V2/mutation_single_gene_level_v6.csv", check.names = TRUE) # from the output of single_gene_analysis.R
  df8 <- merge(df6, df7[, c(2,45)], by.x="Transposed", by.y ="gene")
  colnames(df8)[5] <- "Transposed_NSNPs"
  df9 <- merge(df8, df7[, c(2,45)], by.x="Parental", by.y ="gene")
  colnames(df9)[6] <- "Parental_NSNPs"
  df9 <- select(df9, 2,1,3,4,5,6)
  
  
  df9 <- df9[df9$Transposed_met != "Unclassified",]
  df9 <- df9[df9$Parental_met != "Unclassified",]
  
  df9 <- df9[!is.na(df9$Transposed_met),]
  df9 <- df9[!is.na(df9$Parental_met),]
  
  df9$classification <- NA
  df9$classification <- paste(df9$Transposed_met, "_", df9$Parental_met, sep="")
  df9$abs <- abs(df9$Transposed_NSNPs-df9$Parental_NSNPs)
  
  # Have to do this better but for now, will do it manually - 
  # classification.x and mut_per_base.x is for transposed and .y is for parental
  
  write.csv(df9, "../MBE_V2/Transposed_direction_methylation_NSNPs.csv")
  
  df9 <- read.csv("../Transposed_direction_methylation_NSNPs.csv")
  
  library(ggplot2)
  library(gghalves)
  
  df9$classification <- factor (df9$classification, 
                                 levels = c("gbM_gbM", "unM_unM", "teM_teM","gbM_unM", "unM_gbM", "unM_teM", "gbM_teM", "teM_unM", "teM_gbM"))
  
  p <- ggplot(df9) +
    theme_bw() +
    labs(title = "Transposed duplicates - Polymorphism", subtitle = "")+
    geom_half_boxplot(aes(x=classification,y=abs,fill=classification),
                      side="l",errorbar.draw=TRUE,outlier.color=NA) +
    geom_half_point(aes(x=classification,y=abs,color=classification),
                    side="r",size=0.5) +
    theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=10)) +
    labs(x="Classification",y="Percent Polymorphism (Absolute Difference)")
  p
  
  library(ggpubr)
  
  df11 <- df9[df9$classification!= "gbM_gbM", ]
  df11 <- df11[df11$classification!= "unM_unM", ]
  df11 <- df11[df11$classification!= "teM_teM", ]
  
  df11$classification <- factor (df11$classification, 
                                 levels = c("unM_gbM", "unM_teM", "gbM_teM","gbM_unM", "teM_unM", "teM_gbM"))
  ggpaired(df11, cond1 = "Transposed_NSNPs", cond2 = "Parental_NSNPs", facet.by = "classification",
           palette = "npg", color = "condition", point.size = 1, line.size =0.005, line.color = "grey")
  
# GC content all genes

data16 <- data.frame(
  name=c("unM", "gbM","teM"),
  #name2=c( "C:T", "non-C:T"),
  value=c( 0.448341232,0.439786656, 0.428396055), # value calculated from mutations_single_gene_analysis.csv 
  CI=c(0.000298478,0.000336158,0.001372073) # CI value obtained by prop.test (code below in statistical test section) 
)

data16$name <- factor (data16$name, levels = c("unM", "gbM", "teM"))
library(ggplot2)
library(RColorBrewer)
p<- ggplot(data16, aes(x=name, y=value, fill=name)) +
  theme_bw()+
  geom_bar(stat="identity", 
           position="stack")+
  geom_errorbar(aes(ymin=value-CI, ymax=value+CI), width=0.4, size=0.4, position=position_dodge()) 
p + scale_fill_manual(values = c("unM"="Cornflowerblue", "gbM" = "#00b067", "teM"="#FC8D62")) + ylim(0, 0.5)

#### Statistical Tests ####

#ANOVA tests following https://www.scribbr.com/statistics/anova-in-r/ 

install.packages(c("ggplot2", "ggpubr", "tidyverse", "broom", "AICcmodavg"))
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)

#### ChIP data analysis ####

df10 <- read.csv("PDS5_enrich_03-12-23.csv", check.names = F)[,1:2] #From latest dataset from Grey
df11 <- merge(single_gene_level, df10, by="gene")

df3 <- read.csv("../GeneDuplication_V3/figures_V6/Athaliana/Athaliana_Duplicate_pair_met.csv", header=TRUE)

df4 <- merge(df3, df11[, c(1,65)], by.x="Duplicate.1", by.y ="gene") #data for dup pairs
df5 <- merge(df4, df11[, c(1,65)], by.x="Duplicate.2", by.y ="gene")

df5a <- df5[df5$Duplicate.1_Methylation != "Unclassified",] #removing unclassified
df5a <- df5a[df5a$Duplicate.2_Methylation != "Unclassified",]

df5b <- df5a[c(2,1,3:9)]

write.csv(df5b, "DuplicatePairs_PDS5_03-12-23.csv")

# Manually creating dataframe with values for plotting PDS5 in duplicate pairs from df5b
# t-test significance values are in the manually created file 

# data14 <- data.frame(
#   name=c("gbM-unM", "gbM-unM","unM-teM","unM-teM","gbM-teM", "gbM-teM"),
#   name2=c("Paralog1", "Paralog2"),
#   value=c( 2.301247336, 2.069071976, 1.440033779, 0.088270441,  1.567684995, -0.313787141), 
#   CI=c( 0.036646584, 0.029603799, 0.056664143, 0.099254883,  0.282242585, 0.332743916) 
# )
# 
# data14$name <- factor (data14$name, levels = c("gbM-unM", "unM-teM", "gbM-teM"))
# 
# p<- ggplot(data14, aes(x=name, y=value, fill=name2)) +
#   theme_bw()+
#   geom_bar(stat="identity", 
#            position=position_dodge()) +
#   geom_errorbar(aes(ymin=value-CI, ymax=value+CI), width=.4, size=0.4, position=position_dodge(.9))+
#   labs(x="",y="PDS5-ChIP")+
#   theme(axis.text.y=element_text(color="black",hjust=1,size=12))
# p + scale_fill_brewer(palette = "Paired") 
# 
# 
#### H3K4me1-Duplicate Pairs ####

df10 <- read.csv("PDS5_enrich_03-12-23.csv", check.names = F)[,c(1,3)]
df11 <- merge(single_gene_level, df10, by="gene")

df3 <- read.csv("../GeneDuplication_V3/figures_V6/Athaliana/Athaliana_Duplicate_pair_met.csv", header=TRUE)

df4 <- merge(df3, df11[, c(1,65)], by.x="Duplicate.1", by.y ="gene") #data for dup pairs
df5 <- merge(df4, df11[, c(1,65)], by.x="Duplicate.2", by.y ="gene")

df5a <- df5[df5$Duplicate.1_Methylation != "Unclassified",] #removing unclassified
df5a <- df5a[df5a$Duplicate.2_Methylation != "Unclassified",]

df5b <- df5a[c(2,1,3:9)]

write.csv(df5b, "DuplicatePairs_H3K4me1_03-13-23.csv")

# Manually creating dataframe with values for plotting H3K4me1 in duplicate pairs
# t-test significance values are in the manually created file 

# 
# data14 <- data.frame(
#   name=c("gbM-unM", "gbM-unM","unM-teM","unM-teM","gbM-teM", "gbM-teM"),
#   name2=c("Paralog1", "Paralog2"),
#   value=c( 0.522170088, 0.133711265, -0.924217272, -0.390878546,  0.185870599, -0.88680682),
#   CI=c(0.020922673, 0.020251669, 0.022927231, 0.028879286,  0.137126602, 0.094229551)
# )
# 
# data14$name <- factor (data14$name, levels = c("gbM-unM", "unM-teM", "gbM-teM"))
# 
# p<- ggplot(data14, aes(x=name, y=value, fill=name2)) +
#   theme_bw()+
#   geom_bar(stat="identity",
#            position=position_dodge()) +
#   geom_errorbar(aes(ymin=value-CI, ymax=value+CI), width=.4, size=0.4, position=position_dodge(.9))+
#   labs(x="",y="H3K4me1")+
#   theme(axis.text.y=element_text(color="black",hjust=1,size=12))
# p + scale_fill_brewer(palette = "Paired")


#### H3K4me1 & PDS5 boxplots ####

# Manually creating a csv file to plot enrichment scores for duplicate pairs from df5b

df20 <- read.csv("h3k4me1_pairs_for_boxplot_03-13-23.csv", header=T)

df20$classification <- factor(df20$classification, levels = c("gbM-unM", "unM-teM", "gbM-teM"))

 p <- ggplot(df20, aes(x=classification, y=value, fill=paralog)) +
  theme_bw()+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype=2)+
  scale_fill_brewer(palette = "Paired")
p + ylim(-9,5)

#PDS5

df21 <- read.csv("pds5_pairs_for_boxplot_03-12-23.csv", header=T) # from new data

df21$classification <- factor(df21$classification, levels = c("gbM-unM", "unM-teM", "gbM-teM"))

ggplot(df21, aes(x=classification, y=value, fill=paralog)) +
  theme_bw()+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype=2)+
  scale_fill_brewer(palette = "Paired")

#Transposed pairs boxplot
install.packages("ggpubr")

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

library(ggpubr)

df21 <- read.csv("../Transposed_pairs_NSNPs_boxplot.csv", header=T)

df21$classification <- factor(df21$classification, levels = c("unM -- > gbM", "teM -- > gbM", "gbM -- > unM", "teM -- > unM", "unM -- > teM", "gbM -- > teM",
                                                              "gbM -- > gbM", "teM -- > teM", "unM -- > unM"))

compare_means(value ~ paralog, data = df21, group.by = "classification") #statistical test

 p <- ggplot(df21, aes(x=classification, y=value, fill=paralog)) +
  theme_bw()+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  geom_vline(xintercept = c(2.5, 4.5, 6.5), linetype=2)+
   labs(x="",y="Polymorphism in A thalina") +
   theme(axis.text.y=element_text(color="black",hjust=1,size=12)) #Adding stats on tot he plots. requires ggpubr
 
p + stat_compare_means(label =  "p.signif", label.x = 1.5) 
 
#### Multivariate regression ####
setwd("C:/Users/kench/Dropbox/ANALYSIS/MBE_V2/mutationbias")

df1 <- read.csv("../test_met_level_correlation2.csv")
df1$ABS_CG <- as.numeric(df1$ABS_CG)

df2 <- read.csv("../test_met_level_correlation2_gbm_gbm.csv")
df2$ABS_CG <- as.numeric(df2$ABS_CG)

df3 <- read.csv("../test_met_level_correlation2_unm_unm.csv")
df3$ABS_CG <- as.numeric(df3$ABS_CG)

df4 <- read.csv("../test_met_level_correlation2_tem_tem.csv")
df4$ABS_CG <- as.numeric(df4$ABS_CG)

df5 <- read.csv("../test_met_level_correlation2_diverged.csv")
df5$ABS_CG <- as.numeric(df5$ABS_CG)

df6 <- read.csv("../test_met_level_correlation2_NOT_diverged.csv")
df6$ABS_CG <- as.numeric(df6$ABS_CG)



install.packages("ggpubr")

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

library(ggpubr)

cor.test(x = df1$ABS_CG, y = df1$ABS_NSNPs, method= "spearman")

## multiple regression ##
fit1 <- lm(ABS_NSNPs ~ ABS_CG + ABS_CHG + ABS_CHH + ABS_CG * ABS_CHG + ABS_CG * ABS_CHH + ABS_CHG * ABS_CHH, data = df1)
fit1a <- lm(ABS_NSNPs ~ ABS_CG + ABS_CH + ABS_CG * ABS_CH, data = df1)

anova(fit1a)
summary(fit1)
sqrt(summary(fit1)$r.squared)

fit2 <- lm(ABS_NSNPs ~ ABS_CG + ABS_CHG + ABS_CHH + ABS_CG * ABS_CHG + ABS_CG * ABS_CHH + ABS_CHG * ABS_CHH, data = df2)
fit2a <- lm(ABS_NSNPs ~ ABS_CG + ABS_CH + ABS_CG * ABS_CH, data = df2)

anova(fit2)
summary(fit2)
sqrt(summary(fit2)$r.squared)

fit3 <- lm(ABS_NSNPs ~ ABS_CG + ABS_CHG + ABS_CHH + ABS_CG * ABS_CHG + ABS_CG * ABS_CHH + ABS_CHG * ABS_CHH, data = df3)
fit3a <- lm(ABS_NSNPs ~ ABS_CG + ABS_CH + ABS_CG * ABS_CH, data = df3)

anova(fit3)
summary(fit3)
sqrt(summary(fit3)$r.squared)

fit4 <- lm(ABS_NSNPs ~ ABS_CG + ABS_CHG + ABS_CHH + ABS_CG * ABS_CHG + ABS_CG * ABS_CHH + ABS_CHG * ABS_CHH, data = df4)
fit4a <- lm(ABS_NSNPs ~ ABS_CG + ABS_CH + ABS_CG * ABS_CH, data = df4)

anova(fit4)
summary(fit4)
sqrt(summary(fit4)$r.squared)

fit5 <- lm(ABS_NSNPs ~ ABS_CG + ABS_CHG + ABS_CHH + ABS_CG * ABS_CHG + ABS_CG * ABS_CHH + ABS_CHG * ABS_CHH, data = df5)
fit5a <- lm(ABS_NSNPs ~ ABS_CG + ABS_CH + ABS_CG * ABS_CH, data = df5)

library(car)
fit.type3 <-Anova(lm(ABS_NSNPs ~ ABS_CG + ABS_CHG + ABS_CHH + ABS_CG * ABS_CHG + ABS_CG * ABS_CHH + ABS_CHG * ABS_CHH, data = df5), type ='III')

summary(fit.type3)
pander(fit.type3)

anova(fit5a)
summary(fit5)
sqrt(summary(fit5)$r.squared)

fit6 <- lm(ABS_NSNPs ~ ABS_CG + ABS_CHG + ABS_CHH + ABS_CG * ABS_CHG + ABS_CG * ABS_CHH + ABS_CHG * ABS_CHH, data = df6)
fit6a <- lm(ABS_NSNPs ~ ABS_CG + ABS_CH + ABS_CG * ABS_CH, data = df6)

anova(fit6a)
summary(fit6)
sqrt(summary(fit6)$r.squared)



