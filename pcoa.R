library(RColorBrewer)
library(reshape2)
library(vegan)
library(permute)
library(ellipse)
library(data.table)
library(readxl)
library(writexl)
library(dplyr)
library(tidyverse)
library(cluster)
library(ggplot2)
library(nlme)                    
library(lme4)                    
library(lattice)
library ("plyr")
library(broom)
  


x<-read.csv("9CRC.species.csv",row.names = 1)

#计算距离
pcoa_data_numeric <- x[, -c(1:6)]
pcoa_data_numeric[] <- lapply(pcoa_data_numeric, function(x) as.numeric(as.character(x)))
# Step 1: Calculate Bray-Curtis Distance Matrix
bray_dist <- vegdist(pcoa_data_numeric, method = "bray")
df_dist.matrix<-as.matrix(bray_dist)
write.csv(df_dist.matrix,"10CRC_dist/bray_dist.csv")


##PCOA分析
#使用ape包中的hclust（）函数进行层次聚类，可选择方法有single、complete、median、mcquitty、average 、centroid 、ward
df_hc1 <- hclust(bray_dist,method="average")
#可视化
tiff(filename = "10CRC_dist/cluster-scale.tif",width = 18,height = 100,units ="cm",compression="lzw",bg="white",res=600)
plot(as.dendrogram(df_hc1),type="rectangle",horiz=T)
dev.off ()
# Step 2: Perform PCoA (optional for visualization)
pcoa_results <- cmdscale(bray_dist, eig = TRUE, k = (nrow(pcoa_data_numeric)-1)) # 'k'是{1, 2, ..  n - 1}
pcoa_eig <- pcoa_results$eig
barplot(pcoa_eig)
pcoa_exp <- pcoa_results$eig/sum(pcoa_results$eig)  #各轴解释量#
site <- pcoa_results$point   #样本坐标值#
species <- wascores(pcoa_results$points[,1:2], pcoa_data_numeric)
write.csv(pcoa_exp , '10CRC_dist/pcoa_exp.csv', col.names = NA, quote = FALSE)
write.csv(site, '10CRC_dist/pcoa_sample.csv',  col.names = NA, quote = FALSE)
write.csv(species, '10CRC_dist/pcoa_species.csv',  col.names = NA, quote = FALSE)

# Step 3: PERMANOVA using adonis2，可循环每一个变量
permanova_results <- adonis2(bray_dist ~ x$Cohort, method = "bray")
# View the PERMANOVA results
print(permanova_results)

#PCOA绘图
r<-site[,1:3]
r[,3]<-rownames(r)
colnames(r)[3]<-"ID"
colnames(r)[1]<-"V1"
colnames(r)[2]<-"V2"
x <- x %>% mutate(cohort_group = paste(Cohort, Group, sep = "_"))
meta<-cbind(row.names(x),x$cohort_group) #cohort cohort_group

colnames(meta)<-c("ID","Group")
r2<-merge(meta,r,by="ID")
r2$V1<-as.numeric(r2$V1)
r2$V2<-as.numeric(r2$V2)
#r2$Group <- factor(r2$Group, levels=c("Control","CRC"),ordered=TRUE)

#置信椭圆
calc_ellipse <- function(df, level = 0.95) {
  do.call(rbind, lapply(split(df, df$Group), function(data) {
    ell <- ellipse(cov(data[,c("V1", "V2")]), centre = colMeans(data[,c("V1", "V2")]), level = level)
    data.frame(x = ell[,1], y = ell[,2], Group = unique(data$Group))
  }))
}

# Calculate ellipse data
ellipse_data <- calc_ellipse(r2, level = 0.7)

tiff(filename = "10CRC_dist/pca-Cohort_GROUP.tif",width = 25,height = 15,units ="cm",compression="lzw",bg="white",res=600)
ggplot(r2, aes(V1, V2)) +
  geom_polygon(data = ellipse_data, aes(x, y, fill = Group), alpha = 0.3) +
  geom_point(size = 6, shape = 21, colour = "black", aes(fill = Group)) +
  #scale_fill_manual(values = c("#92AAD6", "#FFB754", "#B75454")) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", face = "bold", size = 20),
    axis.line = element_line(colour = "black", size = 1.5),
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(colour = "black", face = "bold", size = 20)
  ) +
  labs(x = 'PCoA1: 12.9%', y = 'PCoA2: 8.25%')
dev.off ()



