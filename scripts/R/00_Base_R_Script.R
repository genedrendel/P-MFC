##script created for analysis of metagenomic (16S amplicon) data
##script created by J L WOOD 03.12.2016
##we will:
##make a phyloseq object
##learn how to subset a phyloseq object
##perform rarfaction curves
##assess alpha diversity
##assess beta diveristy
##creat exploratory barcharts
##script updated 3.04.2017 by JLWOOD
##script updated 4.05.2020 by JLWOOD - CLR normalisation added 

##PHYLOSEQ OBJECTS_____________________________________________________________________________
##get started: change working dir to top folder of your metagenomic project

##for alpha diversity and NMDS/PCOs we'll rarefy the data
##read in our OTU table and treatment file
##format should be samples as columns and OTUs as rows
##on first run you will need to open the .txt and remove the '#' from the first row
##otherwise R will not read it. also convert 'OTU ID' to 'OTU_ID'
##you can also make a excel file and read in as a .csv

otu_table <- subset(as.data.frame(read.table("08otu_table/raw_readmap.txt", sep="\t", header=TRUE,row.names = "OTU_ID")), select = -taxonomy)


##look at the total reads per sample and decide on a rarefaction depth
rare<-colSums(otu_table)
rare
plot(rare)
##rarefy to depth 4000(we'll be doing this later, we'll lose one sample)

##next we need a tree - we use it for the ordinations
##this will generated in the UNIX script and placed in the 09R_files folder for us
##install the package ape

library(ape)

TREE <- read.tree("09R_files/tree_seqs.phy")

####now a taxonomy table
##make this file in excel from HPC outputs
##headings in the .csv file should be: 'OTU_ID' "Domain' 'Phylum' 'Class' 'Order' 'Family' 'Genus' 'Species'

taxmat <- as.matrix(read.csv("09R_files/tax_table.csv", row.names=1, header=TRUE))

##row names for taxmat and otu table MUST be the same
rownames(taxmat)
rownames(otu_table)

##Load and create a metadata data frame
##run make this in excel:first row should be your sample names, subsequent rows are treatment aspects
##IMPORTANT : you must have more than one column in your treatment file 
##(remember the first column will be pulled out to be used as labels so more than 2 coloums is needed)

treat <- as.data.frame(read.csv("09R_files/mapping_file.csv", row.names=1, header=TRUE))

##OK time to make our phyloseq object
#to install the phyloseq package you will ned to follow instructions on
#https://bioconductor.org/packages/release/bioc/html/phyloseq.html

library(phyloseq)
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
TREAT = sample_data(treat)

OBJ1 = phyloseq(OTU,TAX,TREAT,TREE)

##you can use the following to confrim correct labels/treatments ect have been assigned

sample_data(OBJ1)

##you may wan to remove mitochondria and chloroplasts. you can use this script

library(magrittr)
OBJ1 <- OBJ1 %>%
  subset_taxa(
    Domain == "Bacteria" &
    Family  != "mitochondria" &
    Class   != "Chloroplast"
  )

##the following walkthrough detail many pyloseq preprocessing options
##https://joey711.github.io/phyloseq/preprocess.html

##here is a few typical way of accessing/subseting phyloseq data:
sample_data(OBJ1) #look at sample data table
otu_table(OBJ1) #look at otu table
tax_table(OBJ1) #look at taxonomy table

##subsetting:
OBJ1_rhizo <- subset_samples(OBJ1, Location == "Rhizosphere")
sample_data(OBJ1_rhizo)
OBJ1_sphing = subset_taxa(OBJ1, Order=="Sphingomonadales")

##normalisation opton 1: rarefy data
set.seed(8385) #there is an element of randomness in rarfying. this eliminates that
OBJ1_r <- rarefy_even_depth(OBJ1, sample.size = 3500,  rngseed = TRUE, replace = FALSE, trimOTUs = TRUE)

#normalisation opton 2: transform to TSS
OBJ1_ts = transform_sample_counts(OBJ1, function(OTU) OTU/sum(OTU) )

##ALPHA DIVERSITY_____________________________________________________________________________
##BOXPLOTS, ANOVA, TUKEY TEST

library(ggplot2)
library(microbiome)
##package microbiome needs to be loaded from bioconductor: 
##https://www.bioconductor.org/packages/release/bioc/html/microbiome.html

a <- plot_richness(OBJ1_r, x = "Group")# measures=c("Simpson","Chao1", "Shannon"))#+ theme_bw()
a

##we need to know if data meets normailtiy and homogeneity of varience assumptions
##so we can decide how to test for statistical differences

##make a new object (Div) containing sample information and alpha diversity estimates

r <- estimate_richness(OBJ1_r)#exract the alpha diveristy values
e <- evenness(OBJ1_r, index = "all", zeroes = TRUE)#exract eveness values
treat <- sample_data(OBJ1_r)#extract the OBJ1_r treatment file
Div <- cbind(treat,r,e)

##repeat the below tests for each diveristy metric you are interested in

help(shapiro.test)
shapiro.test(Div$Chao1)
##test for normality
##if p < 0.05 reject null hyp that "samples are normally distributed"
##if not normal dist. use kruskal test
##if data is from normal distribution use anova
##result: p > 0.05 for Shannon and Chao1, but not Simpson


help(bartlett.test)
bartlett.test(Simpson~Group, Div)
##test for homegeneity of varience
##if p < 0.05 reject null hyp that "varience is the same"
##if varience is different use kruskal test
##if varience is same use anova
## result: p > 0.05 for Chao1, Shannon and Simpson


##example ANOVA with post hoc test #parametric
ANOVA1<-aov(Div$Shannon ~ Div$Location * Div$Cd_dose)
summary(ANOVA1)

ANOVA1<-aov(Div$Shannon ~ Div$Group)
summary(ANOVA1)

TUKEY<- TukeyHSD(ANOVA1,'Div$Group', conf.level=0.95)
TUKEY


##WORK IN PROGRESS - AUTOMATING A WAY TO PUT THE LETTERS ON THE BOXPLOTS
# library
library(multcompView)
# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){

   # Extract labels and factor levels from Tukey post-hoc 
     Tukey.levels <- TUKEY[[variable]][,4]
     Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
     
     #I need to put the labels in the same order as in the boxplot :
     Tukey.labels$treatment=rownames(Tukey.labels)
     Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
     return(Tukey.labels)
     }

LABELS <- generate_label_df(TUKEY, "Div$Group")
LABELS

my_colors<-c("green4","red3","darkorange1","green4","red3","darkorange1") #make sure you colours appear in the order that corresponds to your factor levels

# Draw the basic boxplot
a <- boxplot(Div$Shannon ~ Div$Group , ylim=c(min(Div$Shannon), 1.1*max(Div$Shannon)) , col=my_colors[as.numeric(LABELS[,1])] , ylab="Shannon Diversity" , main="")
# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )
text( c(1:nlevels(Div$Group)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )

##End WORK IN PROGRESS



##example Kruskal wallis with post hoc test
kruskal.test(Shannon ~ Group, Div) 

library(FSA)
PT = dunnTest(Shannon ~ Group, data=Div, method="bh")    
PT

##BOXPLOTS 
library(ggplot2)

Div
colnames(Div)#check colnames() to decide what factor you want on your x-axis > we'll use 'Group'
levels(Div$Group) #check the categories in your x-axis facor and the order the will appear > we may need to re-order
lim <- c("Bulk_0","Bulk_20","Bulk_100","Rhizo_0","Rhizo_20","Rhizo_100")#order the axis how you want them
#choose some colours/ colour vectors
fill <- "gold1" #box-plot fill colour
line <- "goldenrod2" #boxplot line colour

#choose indices you are interested in
#-ed out line is to be filled in based on stats

s1 <- ggplot(Div, aes(Group,Shannon)) +  geom_boxplot()
s1 <- s1  + theme_bw()+  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
s1 <- s1 + stat_summary(geom = 'text', label = c("a","ab","a","b","b","ab"), fun.y = max, vjust = -1)
s1 <- s1 + scale_x_discrete(name = "Title", limits= lim, labels = lim)  
s1 <- s1 + scale_y_continuous(name = "Shannon diversty",  breaks = seq(5, 6.5, 0.5), limits=c(5, 6.5))
s1 <- s1 + geom_boxplot(fill = fill, colour = line, alpha = 0.5)
s1


c1 <- ggplot(Div, aes(Group,Chao1)) +  geom_boxplot() 
c1 <- c1 + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
#c1 <- c1 + stat_summary(geom = 'text', label = c("a","a","a","b","b","b"), fun.y = max, vjust = -1)
c1 <- c1 + scale_x_discrete(name = "Title", limits= lim, labels = lab)  
c1 <- c1 + scale_y_continuous(name = "Chao1",  breaks = seq(0, 2000, 500), limits=c(0, 2000))
c1 <- c1 + geom_boxplot(fill = fill, colour = line, alpha = 0.5)
c1


p1 <- ggplot(Div, aes(Group, simpson)) +  geom_boxplot() 
p1 <- p1 + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) 
#p1 <- p1 + stat_summary(geom = 'text', label = c("a","ab","a","b","b","ab"), fun.y = max, vjust = -1)
p1 <- p1 + scale_x_discrete(name = "Title", limits= lim, labels = lab)  
p1 <- p1 + scale_y_continuous(name = "Simpson",  breaks = seq(0,0.5,0.1), limits=c(0,0.5))
p1 <- p1 + geom_boxplot(fill = fill, colour = line, alpha = 0.5)
p1


#make a panel image
library(Rmisc)
multiplot(s1,c1,p1, cols=3)


##BETA DIVERSITY_____________________________________________________________________________
##ORDINATIONS DENDROGRAMS ANOSIM PERMANOVA

##some ordinations first
library("ggplot2")

##NMDS:
NMDS1 <- ordinate(OBJ1_r, "NMDS", distance="unifrac", weighted=TRUE, parallel=TRUE)
NMDS1 #use to check that stress is  < 0.2

##PCoA:
PCoA1 <- ordinate(OBJ1_r, "PCoA", distance="unifrac", weighted=TRUE, parallel=TRUE)
PCoA1$values #use to check eigen values if you wish

p1<-plot_ordination(OBJ1_r, NMDS1, color="Group", shape="Location", label=NULL)
p1

#work out the order of plot labels wiht in the facor we are using (i.e 'Group'):
levels(p1$data$Group)


colvec<-c("green4","red3","darkorange1","green4","red3","darkorange1") #make sure you colours appear in the order that corresponds to your factor levels
black <- rep("black",6) #colour vector for shape outlines

##now start again:
help(ordinate)

NMDS1 <- ordinate(OBJ1_r, "PCoA", distance="unifrac", weighted=TRUE, parallel=TRUE)
p1<-plot_ordination(OBJ1_r, NMDS1 , color="Group", shape="Location", label=NULL) + theme_bw()
print(p1)
p1 <- p1 +  geom_point(aes(colour=factor(Group), fill = factor(Group), shape= factor(Location)), size = 3.5)+ ggtitle(NULL)+ theme(legend.position="none") 
p1 <- p1 + scale_shape_manual(values = c(21,24))
p1 <- p1 + scale_colour_manual(values = black)
p1 <- p1+ scale_fill_manual(values = colvec, breaks = c("Bulk_0", "Bulk_20", "Bulk_100","Rhizo_0", "Rhizo_20", "Rhizo_100"))
p1

NMDS2 <- ordinate(OBJ1_r, "PCoA", distance="unifrac", weighted=FALSE, parallel=TRUE)
p2 <- plot_ordination(OBJ1_r, NMDS2 , color="Group", shape="Location", label=NULL) + theme_bw()
p2 <- p2 +  geom_point(aes(colour=factor(Group), fill = factor(Group), shape= factor(Location)), size = 3.5)+ ggtitle(NULL)+ theme(legend.position="none") 
p2<- p2+ scale_shape_manual(values = c(21,24))
p2<- p2 + scale_colour_manual(values = black)
p2<- p2 + scale_fill_manual(values = colvec, breaks = c("Bulk_0", "Bulk_20", "Bulk_100","Rhizo_0", "Rhizo_20", "Rhizo_100"))
p2

#make a panel image
library(Rmisc)
multiplot(p1,p2, cols=2)


##NOW TO THE STATISTICS - ANOSIM AND PERMANOVA

 
library(vegan)
#pull the variables you want to test:
Location <- get_variable(OBJ1_r, "Location")
Cd_dose <- get_variable(OBJ1_r, "Cd_dose")

##ANOSIMS(good for one-way beta diversity analysis)
##https://sites.google.com/site/mb3gustame/hypothesis-tests/anosim
##Paper explining ANOSIM:
##Non-parametric multivariate analyses of changes in community structure
##Australian Journal of Ecology. (1993) K. R. CLARKE

Group_ano <- anosim(distance(OBJ1_r, "wunifrac"), Location)
Group_ano

Group_ano <- anosim(distance(OBJ1_r, "unifrac"), Location)
Group_ano

##PERMANOVA (function adonis())is similar to anosim but can handle more complex designs
##https://sites.google.com/site/mb3gustame/hypothesis-tests/manova/npmanova

##make distance matricies (these are weighted and unweighted unifrac. can also use 'bray')
OBJ1_r_wu <- distance(OBJ1_r, "wunifrac")
OBJ1_r_u <- distance(OBJ1_r, "unifrac")

Group_ado = adonis(OBJ1_r_wu ~ Location * Cd_dose)
Group_ado

Group_ado = adonis(OBJ1_r_u ~ Location * Cd_dose)
Group_ado

##HEIRACHICAL CLUSTERING____________________________________________________________________
library(vegan)

##make disimialrity matrix
OBJ1_r_b <- distance(OBJ1_r, "bray")


# This is the actual hierarchical clustering call, specifying average-link clustering

OBJ1.hclust <- hclust(OBJ1_r_b, method="average")
plot(OBJ1.hclust)

# Convert hclust into a dendrogram and plot
library(dendextend)
dend <- as.dendrogram(OBJ1.hclust)
plot(dend)


##for each colour scheme work out how many colour we need and
##in what order the colour vecor will need to be

length(levels(get_variable(OBJ1_r, "Group")))
unique(levels(get_variable(OBJ1_r, "Group")))

#make  a vector that asscoiates colurs you want with labels
colorscale<- c("green4","red3", "darkorange1","green4", "red3","darkorange1")
colours_to_use <- colorscale[get_variable(OBJ1_r, "Group")] 

#put them in the right order adn add them to the dend object
colours_to_use <- colours_to_use[order.dendrogram(dend)]
labels_colors(dend) <- colours_to_use

#labels to use:
labs_to_use <- get_variable(OBJ1_r, "Location")
labs_to_use <- labs_to_use[order.dendrogram(dend)]
labels(dend) <- as.vector(labs_to_use)

plot(dend)

# Color in function of the cluster
par(mar=c(1,1,1,7))
dend %>%
  set("labels_col", "black")  %>% #label colours
  set("labels_cex", 0.6)  %>% #label text size
  set("branches_k_color", value = rep(c("#138d75","goldenrod1", "darkorange" ),2), k = 6) %>% #colur k major splits
  set("branches_lwd", 1.8)%>% #branch line weight
  set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 1) %>%  # node point size
  set("leaves_col", labels_colors(dend)) %>% # node point color
  plot(horiz=TRUE, axes=TRUE)
abline(v = 350, lty = 2)



##BARCHARTS________________________________________________________________________________

##PHYLOGENETIC STRUCTURE
##BARCHARTS
#set colour pallet - we'll make a large colur vecors by stringing together
#RColorBrewer colour vectors
#https://www.r-graph-gallery.com/38-rcolorbrewers-palettes/

library(RColorBrewer)
colvec1 <- brewer.pal(11, "RdYlBu")
colvec4 <- brewer.pal(10, "PuOr")
colvec3 <- brewer.pal(10, "BrBG")
colvec2 <- brewer.pal(10, "PiYG")
colvec5 <- brewer.pal(10, "PiYG")
colvec6 <- brewer.pal(10, "BrBG")
colvec <- c(colvec1, colvec2,colvec3,colvec4,colvec5,colvec6) #this is a large colourful vecotor

library(microbiome)

OBJ1_rr <- transform(OBJ1_r, "compositional")#transform to relative abundance -microbiom pkg
OBJ1_rr <- transform_sample_counts(OBJ1_r, function(x) x / sum(x) )#alternative way transform to relative abundance - phyloseq pkg

OBJ1_ord <- tax_glom(OBJ1_rr,taxrank = "Order")#concatenate OTUs at a higher phylogenetic level
OBJ1_ord #check how make OTUs you have - too amny makes an untidy fig legend
#you can either concatenate at a higher tax rank or remove really rare otus thus:
OBJ1_ord2 <- filter_taxa(OBJ1_ord, function(x) mean(x) > 0.0005, TRUE)
OBJ1_ord2 #much more manageable number of OTUs

##PLOT1 - all samples
p1 <- plot_bar(OBJ1_ord2, "Sample", fill="Class")
p1 <- p1 + scale_fill_manual(values=colvec)
p1 <- p1 + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1

##sometimes changes are subtle and hard to see
##one optinon is to plot only OTUs that are differetnially abundant (a lesson for later)


##PLOT2 - all samples faceted by phylum
OBJ1_ord2 <- filter_taxa(OBJ1_ord, function(x) mean(x) > 0.005, TRUE)
p2 <- plot_bar(OBJ1_ord2, "Sample", fill="Class", facet_grid=~Phylum)
p2 <- p2 + scale_fill_manual(values=colvec)
p2 <- p2 + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2


lim <- c("Bulk_0","Bulk_20","Bulk_100","Rhizo_0","Rhizo_20","Rhizo_100")#order the axis how you want them
lim_Cd_dose <- c("0_mg","20_mg","100_mg")#order the axis how you want them

##PLOT3 - Cd-dose groups faceted by Location
OBJ1_ord2 <- filter_taxa(OBJ1_ord, function(x) mean(x) > 0.005, TRUE)
p3 <- plot_bar(OBJ1_ord2, "Cd_dose", fill="Order", facet_grid=~Location)
p3 <- p3 + scale_x_discrete(name = "Title", limits= lim_Cd_dose)  
p3 <- p3 + scale_fill_manual(values=colvec)
p3 <- p3 + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3
#note in this plot abundance goes to 5 b/c there are 5 reps 
#(this is also why each colour has 5 segments) 
#in 0mg rhizo ther are only 4 reps so abundances apper lower

##you can pull out specific OTUs if you have a list of thier names in a vector object
##OBJ1_DE = prune_taxa(VectorOfNames, OBJ1_rr)

#PLOT4 - custom barchart - for when you know what you want to plot
#this method will remove the segmentation in the above plots

OBJ1_rr <- transform(OBJ1_r, "compositional")#transform to relative abundance
OBJ1_ord <- tax_glom(OBJ1_rr,taxrank = "Order")#concatenate to tax level we want to plot
OBJ1_ord2 <- filter_taxa(OBJ1_ord, function(x) mean(x) > 0.005, TRUE)#remove really rare stuff

#using our colour vector again:
library(RColorBrewer)
colvec1 <- brewer.pal(11, "RdYlBu")
colvec4 <- brewer.pal(10, "PuOr")
colvec3 <- brewer.pal(10, "BrBG")
colvec2 <- brewer.pal(10, "PiYG")
colvec5 <- brewer.pal(10, "PiYG")
colvec6 <- brewer.pal(10, "BrBG")
colvec <- c(colvec1, colvec2,colvec3,colvec4,colvec5,colvec6)

melt<-psmelt(OBJ1_ord2)#extract the information from the phyloseq object
head(melt) #check column headings
melt<-melt[sort.list(melt[,10]), ] #reorder 'melt' by the column we will be plotting by (in this cae 'order si in col 10) 

##ordering the x axis labels:
melt$Cd_dose<- as.character(melt$Cd_dose)
#Then turn it back into an ordered factor
melt$Cd_dose <- factor(melt$Cd_dose, levels=unique(melt$Cd_dose))
melt$Cd_dose <- factor(melt$Cd_dose, levels=c("0_mg","20_mg","100_mg"))

t1<-ggplot(melt, aes(Cd_dose, Abundance/5, fill=Order))+ geom_bar(stat = "identity",)
t1 <- t1 + facet_grid(Location~.) + theme(axis.text.x = element_text(angle = 45, size = 9,vjust = 0.7))+ scale_fill_manual(values = colvec)
t1 <- t1+ labs(title = NULL, x=expression(Cadmium~dose~(mg~kg^{-1}~soil)), y= "Relative abundance")
t1 <- t1 + scale_x_discrete(labels= c(rep("0 mg",5),rep("20 mg",5), rep("100 mg",5))) 
t1 

#trends in Sphingomonadales and Sphingobacterales are way clearer
#note we've devided abundance by 5 but the treatment missing a sample will still be skewed
#better to plot the individual sampels if you can
#can do this by adding a custom column in your treamtent file (we have added 'Fig1' and use it in PLOT5)


#PLOT5 - plot of differentially abundant OTUs NB must understand how DESeq model design works FIRST
#https://bioconductor.org/packages/release/bioc/html/DESeq2.html #for downloading the package
#https://www.genomatix.de/online_help/help_regionminer/DESeq2.pdf
#https://lashlock.github.io/compbio/R_presentation.html
#https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf
#https://joey711.github.io/phyloseq-extensions/DESeq2.html

library("DESeq2")

q1 <- subset_samples(OBJ1, Fig1 != "0_mg4")#use unrarefied data for DE so manually remove the sample that we know is undersampled

diagdds = phyloseq_to_deseq2(q1 , ~ Location + Cd_dose) #model for your DE analysis READ VINGITTE FIRST
diagdds$Cd_dose <- relevel(diagdds$Cd_dose, ref="0_mg") ##reset reference level
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
resultsNames(diagdds)

sigtab = results(diagdds, contrast=c("Cd_dose","100_mg","0_mg"))
sigtab_otu = cbind(as(sigtab, "data.frame"), as(tax_table(q1)[rownames(sigtab), ], "matrix"))


dim(sigtab_otu)
sigtab_sigg <- subset(sigtab_otu, sigtab_otu$padj < 0.05) 
dim(sigtab_sigg) #149 OTus
OTU <- rownames(sigtab_sigg)
OBJ1_rr  = transform_sample_counts(OBJ1_r, function(x) x / sum(x) )
OBJ1_DE = prune_taxa(OTU, OBJ1_rr)
OBJ1_DEo <- tax_glom(OBJ1_DE,taxrank = "Order")
OBJ1_DE2 <- subset_taxa(OBJ1_DEo,rowMax(otu_table(OBJ1_DEo)) > 0.005)

library(RColorBrewer)
colvec2 <- brewer.pal(11, "RdYlBu")
colvec3 <- brewer.pal(10, "PuOr")
colvec3 <- brewer.pal(10, "BrBG")
colvec1 <- brewer.pal(10, "PiYG")
colvec5 <- brewer.pal(10, "PiYG")
colvec6 <- brewer.pal(10, "BrBG")
colvec <- c(colvec1, colvec2,colvec3,colvec4,colvec5,colvec6)

melt<-psmelt(OBJ1_DE2)
head(melt)
melt<-melt[sort.list(melt[,11]), ]

##ordering the x axis:

melt$Fig1<- as.character(melt$Fig1)
#Then turn it back into an ordered factor
melt$Fig1 <- factor(melt$Fig1, levels=unique(melt$Fig1))
melt$Fig1 <- factor(melt$Fig1, levels=c("0_mg1","0_mg2","0_mg3","0_mg4","0_mg5","20_mg1","20_mg2","20_mg3","20_mg4","20_mg5","100_mg1","100_mg2","100_mg3","100_mg4","100_mg5"))

t1<-ggplot(melt, aes(Fig1, Abundance/5, fill=Order))+ geom_bar(stat = "identity",)
t1 <- t1 + facet_grid(Location~.) + theme(axis.text.x = element_text(angle = 45, size = 9,vjust = 0.7))+ scale_fill_manual(values = colvec)
t1 <- t1+ labs(title = NULL, x=expression(Cadmium~dose~(mg~kg^{-1}~soil)), y= "Relative abundance")
t1 <- t1 + scale_x_discrete(labels= c(rep("0 mg",5),rep("20 mg",5), rep("100 mg",5))) 
t1 


