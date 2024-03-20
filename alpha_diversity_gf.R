setwd("~/Desktop/GF_16s_1")
getwd()
path <- "/Users/charlotte.kenneally/Desktop/GF_16s_1"

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
library(phyloseq)

## read in files -- specify headers of row names present and that column 1 only contains the row names
otu <-read.csv("OTU600_2.csv", header = TRUE, row.names = 1)
map <-read.csv("samdat600_2.csv", header = TRUE, row.names = 1) 
tax <-read.csv("tax600_GF_Dec1_23.csv", header = TRUE, row.names = colnames(otu))

## reread in files now that i have removed the 1C samples 
otu <-read.csv("OTU600_Jan8_no1Csamples.csv", header = TRUE, row.names = 1)
map <-read.csv("samdat600_Jan8_no1Csamples.csv", header = TRUE, row.names = 1)
tax <-read.csv("tax600_GF_Dec1_23.csv", header = TRUE, row.names = colnames(otu))
tax <- as.matrix(tax[2:8])

ps3 <- phyloseq(otu_table(otu, taxa_are_rows = FALSE), 
                tax_table(tax), sample_data(map)) 
class(map) 
class(dt_ps3)
### make data frame from ps object 
dt_ps3 <- data.table(psmelt(ps3))
### check that it is a data frame
class(dt_sub_goodphyla) ## true 

rowSums(otu_table(ps3)) %>%
  sort()
ps3_cut3000 <- rarefy_even_depth(ps3, rngseed = 777, sample.size = 3000, replace = FALSE)
rowSums(otu_table(ps3_cut3000)) %>%
  sort()

### make smaller data frames for each one i want to display on graphs at once 
donor <- subset_samples(ps3_cut3000, sample_type == "donor")
mouse_stools <- subset_samples(ps3_cut3000, sample_type == "stool")
round1 <- subset_samples(ps3_cut3000, round == "1")
round2 <- subset_samples(ps3_cut3000, round == "2")
round3 <- subset_samples(ps3_cut3000, round == "3")
rec_pooled1 <- subset_samples(ps3_cut3000, donor_type == "pooled1")
cr_rec <- subset_samples(ps3_cut3000, donor == "Cr")
mr1ko_rec <- subset_samples(ps3_cut3000, donor == "mr1ko")
bl6_rec <- subset_samples(ps3_cut3000, donor == "BL6-1")
controls <- subset_samples(ps3_cut3000, sample_type == "ctrl")

### make it a data.frame
dt_donor <- data.table(psmelt(donor))
dt_mouse_stools <- data.table(psmelt(mouse_stools))
dt_round1 <- data.table(psmelt(round1))
dt_round2 <- data.table(psmelt(round2))
dt_round3 <- data.table(psmelt(round3))
dt_rec_pooled1 <-data.table(psmelt(rec_pooled1))
dt_cr_rec <- data.table(psmelt(cr_rec))
dt_mr1ko_rec <- data.table(psmelt(mr1ko_rec))
dt_bl6_rec <- data.table(psmelt(bl6_rec))
dt_controls <- data.table(psmelt(controls))

### alpha diversity 
obs_richness <- estimate_richness(ps3_cut3000, measures = "Observed")
shannon_richness <- estimate_richness(ps3_cut3000, measures = "Shannon")
simpson_richness <- estimate_richness(ps3_cut3000, measures = "Simpson")
inv_richness <- estimate_richness(ps3_cut3000, measures = "InvSimpson")

### observed richness all
plot_richness(ps3_cut3000, x = "cage", measures = "Observed", color = "cage")+
  geom_boxplot(alpha = 0.6)+
  theme()

### simpson all 
plot_richness(ps3_cut3000, x = "cage", measures = "Simpson", color = "cage")+
  geom_boxplot(alpha = 0.6)+
  theme()

## shannon boxplot with every sample 
plot_richness(ps3_cut3000, x = "cage", measures = "Shannon", color = "cage")+
  geom_boxplot(alpha = 0.6)+
  theme() 
### inverse
plot_richness(ps3_cut3000, x = "cage", measures = "InvSimpson", color = "cage")+
  geom_boxplot(alpha = 0.6)+
  theme() 

## richness of donors alone
plot_richness(donor, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6)+
  labs(title = "Shannon Diversity all donors")+
  theme()
plot_richness(donor, x = "cage", measures = "InvSimpson", color = "cage")+ 
  geom_boxplot(alpha = 0.6)+
  labs(title = "Inverse Simpson Diversity all donors")+
  theme()
plot_richness(donor, x = "cage", measures = "Observed", color = "cage")+ 
  geom_boxplot(alpha = 0.6)+
  labs(title = "Observed Richness all donors")+
  theme()

### recipients of pooled 1 
plot_richness(rec_pooled1, x = "cage", measures = "Observed", color = "cage")+
  geom_boxplot()+
  labs(title = "Observed Richness of Recipients of Pooled 1")+
  theme_classic() 
plot_richness(rec_pooled1, x = "cage", measures = "Observed", color = "cage")+
  geom_violin()+
  labs(title = "Observed richness of recipients of pooled 1")+
  theme_classic()
### ABUNDANCE BAR PLOTS
# barplot abundance pool 1
dt_rec_pooled1 %>%
  ggplot(aes(cage, Abundance, color = Phylum))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")
# barplot abundance all stools
dt_mouse_stools %>%
  ggplot(aes(cage, Abundance, color = Phylum))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")
# barplot abundance all donors
dt_donor %>%
  ggplot(aes(cage, Abundance, color = Phylum))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")

## barplot abundance of stools grouped by round 1 and 2 only 
a = dt_round1 %>%
  drop_na(Phylum) %>%
  ggplot(aes(cage, Abundance, color = Phylum))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("Round 1")
b = dt_round2 %>% 
  drop_na(Phylum) %>%
  ggplot(aes(cage, Abundance, color = Phylum))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("Round 2")
ggarrange(a,b, nrow = 1, common.legend = TRUE, legend = "right")

# dt_mouse_stools %>%
  drop_na(Phylum) %>%
  ggplot(aes(round, Abundance, color = Phylum))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("Round 1")+
  facet_wrap(~cage)

### do that again and try to replace x-axis with the donor names 
c = dt_round1 %>%
  drop_na(Phylum) %>%
  ggplot(aes(cage, Abundance, color = Phylum))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  scale_x_discrete(labels=c("PDxCre(KO)", "Mr1KO", "Pooled 1"))+
  theme_classic()+
  ggtitle("Round 1")
d = dt_round2 %>%
  drop_na(Phylum) %>%
  ggplot(aes(cage, Abundance, color = Phylum))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  scale_x_discrete(labels=c("PDxCre(KO)", "C57BL6", "Mr1KO", "Pooled 1"))+
  theme_classic()+
  labs(title="Round 2")
ggarrange(c,d, nrow = 1, common.legend = TRUE, legend = "right")


## with all the stool minus pooled mix 2
dt_rd1_2_3B %>%
  drop_na(Phylum) %>%
  ggplot(aes(cage, Abundance))+ ##color = Phylum
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  scale_x_discrete(name = "Donor", limits = c("1A", "2A", "1B", "2C", "2B", "1D", "2D", "3B"), 
                   labels=c("1A"="PDxCre(KO)", "2A"="PDxCre(KO)", "1B"="Mr1KO", 
                            "2C"="Mr1KO", "2B"="C57BL6", "1D"="Pooled 1", "2D"="Pooled 1", "3B"="Pooled 1"))+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF","#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
    ## "#5991FFFF","#8CB2FFFF","#BFD4FFFF","#E6EEFFFF","#F7FAFFFF","#FFFFCCFF","#FFFF99FF","#FFFF00FF","#FFCC00FF","#FF9900FF","#FF6600FF","#FF0000FF"))+
    #"#8E0152FF","#C51B7DFF","#DE77AEFF","#F1B6DAFF","#FDE0EFFF","#F7F7F7FF","#E6F5D0FF","#B8E186FF","#7FBC41FF","#4D9221FF","#276419FF","#417839FF"))+ 
                               ## "#A7D3D4FF", "#F8D564FF", "#EDD6D5FF", "#EA879CFF", "#FF3D7FFF", "#7FC7AFFF", "#E9E29CFF",
                              ##  "#E0E4CCFF", "#F38630FF", "#FA6900FF", "#BF9BDDFF"))+
  ## paletteer::scale_fill_paletteer_d("wesanderson::FantasticFox1")+
  ## scale_fill_brewer(palette = "PiYG")+
  labs(title = "Relative Abundance")

display.brewer.all()

install.packages("paletteer")
library(paletteer)

### received cr:NIH
plot_richness(cr_rec, x = "cage", measures = "Observed", color = "cage")+
  geom_boxplot()+
  labs(title = "Observed Richness of Mice who received Cr:NIH")+
  theme_classic()

### received mr1ko
plot_richness(mr1ko_rec, x = "cage", measures = "Observed", color = "cage")+
  geom_boxplot()+
  labs(title = "Observed Richness of Mice who received Mr1KO")+
  theme_classic()

### received C57BL6
plot_richness(bl6_rec, x = "cage", measures = "Observed", color = "cage")+
  geom_boxplot()+
  labs(title = "Observed Richness of Mice who received C57BL6")+
  theme_classic()

### richness stool alone 
plot_richness(mouse_stools, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6)+
  theme_classic()+
  labs(title = "Shannon Diversity all mouse stool samples")
plot_richness(mouse_stools, x = "cage", measures = "Simpson", color = "cage")+ 
  geom_boxplot(alpha = 0.6)+
  theme()
plot_richness(mouse_stools, x = "cage", measures = "Observed", color = "cage")+ 
  geom_boxplot(alpha = 0.6)+
  theme()
#plot_richness(rec_pooled1, x = "cage", measures = "Observed", color = "cage")+
  geom_violin(rec_pooled1, "Observed", add = "boxplot", fill = "cage")+
  labs(title = "Observed richness of recipients of pooled 1")+
  theme_classic()

#dt_mouse_stools%>%
  ggplot(aes(cage, "Observed", color = "cage"))+
  geom_violin(stat = "ydensity", position = "dodge")+
  geom_boxplot(aes(stat = "identity", position = "fill"))+
  theme_classic()

#dt_mouse_stools%>%
  ggplot(aes(cage, "Observed", color = "cage"))+
  geom_violin(aes(stat = "ydensity", position = "dodge", fill = "#a6cee3", "#b2df8a", "#fdbf6f"))+
  theme_classic()

  
### stats
  # compare pooled 1 groups diversity (1D, 2D, 3B)
  # also compare the single donors with same donor between 1 and 2
library("ggpubr")
a_my_comparisons <- list(c("1D", "2D", "3B"), c("1A", "2A"), c("1B", "2C"), c("2D", "3B"), c("1D", "3B"))
symnum = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
plot_richness(mouse_stools, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6)+
  theme_classic()+
  scale_x_discrete(labels=c("1A"="PDxCre(KO)", "1B"="Mr1KO", "1D"="Pooled 1", "2A"="PDxCre(KO)", "2B"="C57BL6", "2C"="Mr1KO", "2D"="Pooled 1", "3A"="Pooled 2", "3B"="Pooled 1", "3D"="Pooled 2"))+
  ## scale_x_discrete(limits = c("PDxCre(KO)", "Mr1KO", "Pooled 1", "PDxCre(KO)", "C57BL6", "Mr1KO", "Pooled 1", "Pooled 2", "Pooled 1", "Pooled 2"))+
  labs(title = "Shannon Diversity all mouse stool samples")+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum)

## again but with 3B only and do not include 3A or 3D bc pooled 2 
## make PS object with just rounds 1, rounds2, and round 3 cage 3B 
rd1_2_3B <- subset_samples(ps3_cut3000, sample_type == "stool")
dt_rd1_2_3B <- data.table(psmelt(rd1_2_3B))

rd1_2_3B <- subset_samples(rd1_2_3B, cage != "3A")
rd1_2_3B <- subset_samples(rd1_2_3B, cage != "3D")

## plot with significance between donors 
plot_richness(rd1_2_3B, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6)+
  theme_classic()+
  scale_x_discrete(name = "Donor", limits = c("1A", "2A", "1B", "2C", "2B", "1D", "2D", "3B"), labels=c("1A"="PDxCre(KO)", "2A"="PDxCre(KO)", "1B"="Mr1KO", 
                  "2C"="Mr1KO", "2B"="C57BL6", "1D"="Pooled 1", "2D"="Pooled 1", "3B"="Pooled 1"))+
  ## scale_x_discrete(limits = c("PDxCre(KO)", "Mr1KO", "Pooled 1", "PDxCre(KO)", "C57BL6", "Mr1KO", "Pooled 1", "Pooled 2", "Pooled 1", "Pooled 2"))+
  labs(title = "Shannon Diversity")+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum)
  
### all richness round1 together

a = plot_richness(round1, x = "cage", measures = "Observed", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("ASV Richness") +
  scale_x_discrete(labels = c("1A"="PDxCre(KO)", "1B"="Mr1KO", "1D"="Pooled 1"))
## a + scale_x_discrete(name = c("1A"="PDxCre(KO)", "1B"="Mr1KO", "1D"="Pooled 1"))
b = plot_richness(round1, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Shannon") 
c = plot_richness(round1, x = "cage", measures = "Simpson", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Simpson") 

ggarrange(a,b,c, nrow = 1, common.legend = TRUE)

## with stats
round1_alpha_comparisions <- list(c("1A", "1D"), c("1B", "1D"))
a = plot_richness(round1, x = "cage", measures = "Observed", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("ASV Richness")+
  stat_compare_means(method = "wilcox.test", comparisons = round1_alpha_comparisions, label = "p.signif", symnum.args = symnum)
b = plot_richness(round1, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Shannon")+
  stat_compare_means(method = "wilcox.test", comparisons = round1_alpha_comparisions, label = "p.signif", symnum.args = symnum)
c = plot_richness(round1, x = "cage", measures = "Simpson", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Simpson")+
  stat_compare_means(method = "wilcox.test", comparisons = round1_alpha_comparisions, label = "p.signif", symnum.args = symnum)

ggarrange(a,b,c, nrow = 1, common.legend = TRUE)
  
### all richness round2 together 
a = plot_richness(round2, x = "cage", measures = "Observed", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("ASV Richness")
b = plot_richness(round2, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Shannon") 
c = plot_richness(round2, x = "cage", measures = "Simpson", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Simpson") 

ggarrange(a,b,c, nrow = 1, common.legend = TRUE)
# WITH STATS
round2_alpha_comparisions <- list(c("2A", "2D"), c("2B", "2D"), c("2C", "2D"))
a = plot_richness(round2, x = "cage", measures = "Observed", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("ASV Richness")+
  stat_compare_means(method = "wilcox.test", comparisons = round2_alpha_comparisions, label = "p.signif", symnum.args = symnum)
b = plot_richness(round2, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Shannon")+
  stat_compare_means(method = "wilcox.test", comparisons = round2_alpha_comparisions, label = "p.signif", symnum.args = symnum)
c = plot_richness(round2, x = "cage", measures = "Simpson", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Simpson")+
  stat_compare_means(method = "wilcox.test", comparisons = round2_alpha_comparisions, label = "p.signif", symnum.args = symnum)

ggarrange(a,b,c, nrow = 1, common.legend = TRUE)
  
### round 3
a = plot_richness(round3, x = "cage", measures = "Observed", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("ASV Richness")
b = plot_richness(round3, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Shannon") 
c = plot_richness(round3, x = "cage", measures = "Simpson", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Simpson") 

ggarrange(a,b,c, nrow = 1, common.legend = TRUE)
 ## not doing any tests for sigf bc its all pooled, not single compared to multi in this round 

### richness donors?? 
a = plot_richness(donor, x = "cage", measures = "Observed", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("ASV Richness")
b = plot_richness(donor, x = "cage", measures = "Shannon", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Shannon") 
c = plot_richness(donor, x = "cage", measures = "Simpson", color = "cage")+ 
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  theme_classic()+
  ylab("Simpson") 

ggarrange(a,b,c, nrow = 1, common.legend = TRUE)



