# Intro and Starting Notes ------------------------------------------------
##script created for analysis of metagenomic (16S amplicon) data
##script created by J L WOOD 03.12.2016
##we will:
##make a phyloseq object
##learn how to subset a phyloseq object
##perform rarefaction curves
##assess alpha diversity
##assess beta diveristy
##create exploratory barcharts
##script updated 3.04.2017 by JLWOOD
##script updated 4.05.2020 by JLWOOD - CLR normalisation added
##Script updated and modified by Gene Drendel 2020 - Total Sum Scaling and P-MFC project specific changes
#23.11.2020 Script now uploaded to P-MFC github repo, changes to be tracked and maintained there

##PHYLOSEQ OBJECTS_____________________________________________________________________________
##get started: change working dir to top folder of your metagenomic project
##for alpha diversity and NMDS/PCOs we'll rarefy the data
##read in our OTU table and treatment file
##format should be samples as columns and OTUs as rows
##on first run you will need to open the .txt and remove the '#' from the first row
##otherwise R will not read it. also convert 'OTU ID' to 'OTU_ID'
##you can also make a excel file and read in as a .csv

## Working Directory -------------------------------------------------------
#Set your working directory, location for data tables, tree, etc
setwd("~/Documents/University/Analysis/PMFC_18/2020 rerun outputs/Format for phyloseq")

## Import ------------------------------------------------------------------
#Import .csv as OTU table
#Replaced import with the below: Not sure if any functional difference once imported like this? But the other one didn't work, assuming mostly just a txt vs csv thing
otu_table <- as.data.frame(read.csv("raw_readmap.csv", header=TRUE,row.names = "OTU_ID"))

##look at the total reads per sample and decide on a rarefaction depth
rare<-colSums(otu_table)
rare
plot(rare)
##rarefy to depth 4000(we'll be doing this later, we'll lose one sample)

##next we need a tree - we use it for the ordinations
##this will generated in the UNIX script and placed in the 09R_files folder for us

library(ape)

TREE <- read.tree("rooted_tree.nwk")

####now a taxonomy table
##make this file in excel from HPC outputs
##headings in the .csv file should be: 'OTU_ID' "Domain' 'Phylum' 'Class' 'Order' 'Family' 'Genus' 'Species'

taxmat <- as.matrix(read.csv("tax_table.csv", row.names=1, header=TRUE))

##row names for taxmat and otu table MUST be the same
rownames(taxmat)
rownames(otu_table)

##Load and create a metadata data frame
##run make this in excel:first row should be your sample names, subsequent rows are treatment aspects
##IMPORTANT : you must have more than one column in your treatment file
##(remember the first column will be pulled out to be used as labels so more than 2 coloums is needed)

treat <- as.data.frame(read.csv("mapping_file.csv", row.names=1, header=TRUE))

##OK time to make our phyloseq object
#to install the phyloseq package you will ned to follow instructions on
#https://bioconductor.org/packages/release/bioc/html/phyloseq.html

library(phyloseq)
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
TREAT = sample_data(treat)
OBJ1 = phyloseq(OTU,TAX,TREAT,TREE)

##you can use the following to confirm correct labels/treatments ect have been assigned

sample_data(OBJ1)


## Remove unwanted taxa ----------------------------------------------------
##you may wan to remove mitochondria and chloroplasts. you can use this script

library(magrittr)
OBJ1 <- OBJ1 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
    Family  != "mitochondria" &
    Class   != "Chloroplast"
  )


##the following walkthrough detail many phyloseq preprocessing options
##https://joey711.github.io/phyloseq/preprocess.html

##here is a few typical way of accessing/subseting phyloseq data:
sample_data(OBJ1) #look at sample data table
otu_table(OBJ1) #look at otu table
tax_table(OBJ1) #look at taxonomy table

#Install and run Shiny Phyloseq for interactive exploration of data (however can only import as .biom, .rdata or netwick trees)
install.packages("shiny")
shiny::runGitHub("shiny-phyloseq","joey711")


# Subsetting --------------------------------------------------------------
##Subsetting:
#First thing first GET RID OF THE CONTROL SAMPLE AND ONLY SUBSET FROM THAT OBJECT SO THAT IT STOPS FUCKING WITH ORDINATION DISTANCES
#Similarly, will use this to remove probleamtic samples , e.g those with low read depth that entirely skew ordination patterns by being entirely alone
#To do this, created additional column of treatment file indicating Y (yes) and N (no) for all data minus the control
#Lets call this OBJ1_exp as an indication of it being the proper "experimental" dataset


## Experimental Sample Subset (and TSS it) ----------------------------------------------
#Main Experimental Data Subset - run this one EVERY TIME
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
#TSS Transform on ALL Experimental data as pre-treatment for ordinations (rather than doing TSS for each individual subset)
#Overall TSS - then run the individual subsets, BUT ensure these are off of the TSS set, not the regular set in the subsetting section
OBJ1_exp_tss = transform_sample_counts(OBJ1_exp, function(OTU) OTU/sum(OTU) )
OBJ1_exp_tss

#Subsets on TSS data

## Treatment Subsets -------------------------------------------------------
#Treatments (by connection) - for trying to get all the info on treatment combos in one (e.g keep in mind for ordinations = using open and closed symbols + colour
#to distinguish the significant inoculum+connection effect that seems to be happening based on the permanova
OBJ_Unin_Conn_tss <- subset_samples(OBJ1_exp_tss, Treatment == "Uninoculated Connected")
OBJ_Unin_Unconn_tss <- subset_samples(OBJ1_exp_tss, Treatment == "Uninoculated Unconnected")
OBJ_Geo_Conn_tss <- subset_samples(OBJ1_exp_tss, Treatment == "Geobacter Connected")
OBJ_Geo_Unconn_tss <- subset_samples(OBJ1_exp_tss, Treatment == "Geobacter Unconnected")
OBJ_Mont_Conn_tss <- subset_samples(OBJ1_exp_tss, Treatment == "Montebello Connected")
OBJ_Mont_Unconn_tss <- subset_samples(OBJ1_exp_tss, Treatment == "Montebello Unconnected")
OBJ_Pseudo_Conn_tss <- subset_samples(OBJ1_exp_tss, Treatment == "Pseudomonas Connected")
OBJ_Pseudo_Unconn_tss <- subset_samples(OBJ1_exp_tss, Treatment == "Pseudomonas Unconnected")

#Locations
OBJ1_anode_tss <- subset_samples(OBJ1_exp_tss, Location == "Anode")
OBJ1_cathode_tss <- subset_samples(OBJ1_exp_tss, Location == "Cathode")

#Connection
OBJ1_connected_tss <- subset_samples(OBJ1_exp_tss, Connection == "Connected")
OBJ1_unconnected_tss <- subset_samples(OBJ1_exp_tss, Connection == "Unconnected")

#Week/Time
OBJ_W0_tss <- subset_samples(OBJ1_exp_tss, Week == "Zero")
OBJ_W6_tss <- subset_samples(OBJ1_exp_tss, Week == "Six")
OBJ_W8_tss <- subset_samples(OBJ1_exp_tss, Week == "Eight")
OBJ_W14_tss <- subset_samples(OBJ1_exp_tss, Week == "Fourteen")

#Treatment
OBJ_Unin_tss <- subset_samples(OBJ1_exp_tss, Inoculum == "Uninoculated")
OBJ_Geo_tss <- subset_samples(OBJ1_exp_tss, Inoculum == "Geobacter")
OBJ_Mont_tss <- subset_samples(OBJ1_exp_tss, Inoculum == "Montebello")
OBJ_Pseudo_tss <- subset_samples(OBJ1_exp_tss, Inoculum == "Pseudomonas")

## Non-transformed subsets ------------------------------------
#Subsets for NON transformed dataset
#Retained for if used/needed for non-beta

#Locations
OBJ1_anode <- subset_samples(OBJ1_exp, Location == "Anode")
OBJ1_cathode <- subset_samples(OBJ1_exp, Location == "Cathode")

#Connection
OBJ1_connected <- subset_samples(OBJ1_exp, Connection == "Connected")
OBJ1_unconnected <- subset_samples(OBJ1_exp, Connection == "Unconnected")

#Week/Time
OBJ_W0 <- subset_samples(OBJ1_exp, Week == "Zero")
OBJ_W6 <- subset_samples(OBJ1_exp, Week == "Six")
OBJ_W8 <- subset_samples(OBJ1_exp, Week == "Eight")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")

sample_data(OBJ1_anode)


## Subset/Agglomerate taxa -------------------------------------------------
#e.g subset for taxa too (ATTENTION RE:still need to change for tss vs norm if goin to do this for ordinations vs boxplots etc)
#Experimental
OBJ1_Desulf = subset_taxa(OBJ1_exp_tss, Order=="Desulfuromonadales")
OBJ1_Geobacteraceae = subset_taxa(OBJ1_exp_tss, Family=="Geobacteraceae")
OBJ1_Geobacter = subset_taxa(OBJ1_exp_tss, Genus=="Geobacter")

#Agglomerating at certain taxonomic levels
#Will be basing off of TSS'd Data objects, currently mainly interested in end point stats, so W14
OBJ_W14_tss_ORD <- tax_glom(OBJ_W14_tss,taxrank = "Order")
OBJ_W14_tss_FAM <- tax_glom(OBJ_W14_tss,taxrank = "Family")
OBJ_W14_tss_GEN <- tax_glom(OBJ_W14_tss,taxrank = "Genus")


# Rarefy / TSS ---------------------------------------------
##normalisation opton 1: rarefy data
set.seed(8385) #there is an element of randomness in rarfying. this eliminates that
OBJ1_r <- rarefy_even_depth(OBJ1, sample.size = 3500,  rngseed = TRUE, replace = FALSE, trimOTUs = TRUE)

#TSS Total Sum Scaling , do these before beta analyses, at the absolute least i.e use for ordinations
#Perfectly valid to use for alpha diversity too though apparently!
#e.g Diversity measures generally work on proportions so that's not a problem, but may be best to multiply by 100 and round to full numbers, mostly thinking about normality tests here)
#normalisation opton 2: transform to TSS (Total Sum Scaling)
#Overall
OBJ1_ts = transform_sample_counts(OBJ1, function(OTU) OTU/sum(OTU) )
OBJ1_ts

#Can doublecheck these again for transformed or normalised data
sample_data(OBJ1_ts) #look at sample data table
otu_table(OBJ1_ts) #look at otu table
tax_table(OBJ1_ts) #look at taxonomy table

#AFTER MEETING W JEN
#THIS SCRIPT TO CUT OUT/SUBSET SAMPLES WIHTOUT THE PESKY OUTLIER SAMPLES THAT SKEW THE ORDINATIONS
OBJ1_s <- subset_samples(OBJ1_pmta2, PortStar != "PORT") #assuming for e.g Portstar = Sanmple_Name and "PORT" = "W6UG3A", or Location and "Root"

#Normality tests don't seem to like the TSS data, have to convert to integers?
OBJ1_ts_rounded = transform_sample_counts(OBJ1_ts, function(OTU) round(OTU,digits=0) )
OBJ1_ts_rounded

#Re: TSS proportions and integers
#Microbiome analysis says can multiply by 1,000,000 for easier interpretation...let's test..?
OBJ1_ts_multiplied = transform_sample_counts(OBJ1_ts, function(OTU) OTU*1000000 )
otu_table(OBJ1_ts_multiplied)
#Josh Suggested trying by 100 insetad....may as well try a few variations 1000 etc... but question remains what should really be done
#e.g are shannons and simpsons impacted ONLY by ratios, or do they need whole numbers? Need to test some comparisons and see what other people do

# Alpha Diversity ---------------------------------------------------------
##ALPHA DIVERSITY_____________________________________________________________________________
##BOXPLOTS, ANOVA, TUKEY TEST

# Update to latest version as of 3/2/21
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("microbiome")

library(ggplot2)
library(microbiome)
##package microbiome needs to be loaded from bioconductor:
##https://www.bioconductor.org/packages/release/bioc/html/microbiome.html

#RAW
a <- plot_richness(OBJ1, x = "Location")# measures=c("Simpson","Chao1", "Shannon"))#+ theme_bw()
a

#TSS Rounded
a <- plot_richness(OBJ1_ts_rounded, x = "Location")# measures=c("Simpson","Chao1", "Shannon"))#+ theme_bw()
a

#TSS Multiplied
a <- plot_richness(OBJ1_ts_multiplied, x = "Location")# measures=c("Simpson","Chao1", "Shannon"))#+ theme_bw()
a


## Normality Tests ---------------------------------------------------------


##we need to know if data meets normailtiy and homogeneity of varience assumptions
##so we can decide how to test for statistical differences

##make a new object (Div) containing sample information and alpha diversity estimates

r <- estimate_richness(OBJ1_ts_rounded)#extract the alpha diversity values
e <- evenness(OBJ1_ts_rounded, index = "all", zeroes = TRUE)#extract evenness values
treat <- sample_data(OBJ1_ts_rounded)#extract the OBJ1_r treatment file
Div <- cbind(treat,r,e)

##repeat the below tests for each diversity metric you are interested in

help(shapiro.test)
shapiro.test(Div$Chao1)
##test for normality
##if p < 0.05 reject null hyp that "samples are normally distributed"
##if not normal dist. use kruskal test
##if data is from normal distribution use anova
##result: p > 0.05 for Shannon and Chao1, but not Simpson


help(bartlett.test)
bartlett.test(Simpson~Location, Div)
##test for homegeneity of varience
##if p < 0.05 reject null hyp that "varience is the same"
##if varience is different use kruskal test
##if varience is same use anova
## result: p > 0.05 for Chao1, Shannon and Simpson


##example ANOVA with post hoc test #parametric
ANOVA1<-aov(Div$Shannon ~ Div$Location * Div$Connection)
summary(ANOVA1)

ANOVA1<-aov(Div$Shannon ~ Div$Location)
summary(ANOVA1)

TUKEY<- TukeyHSD(ANOVA1,'Div$Location', conf.level=0.95)
TUKEY


## WIP - Boxplot labels ----------------------------------------------------


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

LABELS <- generate_label_df(TUKEY, "Div$Location")
LABELS

my_colors<-c("green4","red3","darkorange1","green4","red3","darkorange1") #make sure you colours appear in the order that corresponds to your factor levels

# Draw the basic boxplot
a <- boxplot(Div$Shannon ~ Div$Location , ylim=c(min(Div$Shannon), 1.1*max(Div$Shannon)) , col=my_colors[as.numeric(LABELS[,1])] , ylab="Shannon Diversity" , main="")
# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )
text( c(1:nlevels(Div$Group)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )

##End WORK IN PROGRESS

##example Kruskal wallis with post hoc test
kruskal.test(Shannon ~ Location, Div)

library(FSA)
PT = dunnTest(Shannon ~ Location, data=Div, method="bh")
PT


## Boxplots ----------------------------------------------------------------

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


# Heatmaps ----------------------------------------------------------------
##PHYLOSEQ HEATMAPS (WIP)_____________________________________________
##INPROGRESS

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())

#Overall
gpt_all <- subset_taxa(OBJ1_exp_tss, Kingdom=="Bacteria")
gpt_all <- prune_taxa(names(sort(taxa_sums(gpt_all),TRUE)[1:30]), gpt_all)
plot_heatmap(gpt_all, sample.label="Location")
plot_heatmap(gpt_all, "NMDS", "unifrac", "Connection", "Family", sample.label="Location", weighted=FALSE)
plot_heatmap(gpt_all, "NMDS", "unifrac", "Connection", "Family", sample.label="Location", weighted=TRUE)
heatmap(otu_table(gpt))

gpt_prot <- subset_taxa(OBJ1_exp_tss, Phylum=="Proteobacteria")
gpt_prot <- prune_taxa(names(sort(taxa_sums(gpt_prot),TRUE)[1:30]), gpt_prot)
plot_heatmap(gpt_prot)
plot_heatmap(gpt_prot, "NMDS", "unifrac", "Connection", "Family", weighted=FALSE)
plot_heatmap(gpt_prot, "NMDS", "unifrac", "Connection", "Family", weighted=TRUE)
plot_heatmap(gpt_prot)

#Subsetting...test to try to get labels legible and see time series
gpt_sub_an <- subset_taxa(OBJ1_anode_tss, Kingdom=="Bacteria")
gpt_sub_an <- prune_taxa(names(sort(taxa_sums(gpt_sub_an),TRUE)[1:30]), gpt_sub_an)
plot_heatmap(gpt_sub_an, sample.label="Connection")
plot_heatmap(gpt_sub_an, method = "NMDS", distance = "unifrac", sample.label = "Connection", taxa.label = "Genus")
heatmap(otu_table(gpt_sub_an))

#Subsetting w14...test 
gpt_sub14 <- subset_taxa(OBJ_W14_tss, Kingdom=="Bacteria")
gpt_sub14 <- prune_taxa(names(sort(taxa_sums(gpt_sub14),TRUE)[1:30]), gpt_sub14)
plot_heatmap(gpt_sub14, method = "MDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.label = "Family", weighted=TRUE)
heatmap(otu_table(gpt_sub14))

#Try everything by time for those asvs
gpt_alltime <- subset_taxa(OBJ1_exp_tss, Kingdom=="Bacteria")
gpt_alltime <- prune_taxa(names(sort(taxa_sums(gpt_alltime),TRUE)[1:30]), gpt_alltime)
plot_heatmap(gpt_alltime, method = "MDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted=TRUE)
heatmap(otu_table(gpt_alltime))

#Now agglomerate to diff taxa levels instead of asv
OBJ_all_tss_FAM <- tax_glom(OBJ1_exp_tss,taxrank = "Family")
OBJ_all_tss_FAM <- subset_taxa(OBJ_all_tss_FAM, Kingdom=="Bacteria")
OBJ_all_tss_FAM <- prune_taxa(names(sort(taxa_sums(OBJ_all_tss_FAM),TRUE)[1:30]), OBJ_all_tss_FAM)
plot_heatmap(OBJ_all_tss_FAM, method = "MDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted=TRUE)
heatmap(otu_table(OBJ_all_tss_FAM))


## Heatmap for specialisation index -------------------------------------
#Seperating out most variable asvs across trial for paper figure
#Import heatmap subset .csv's (in this case: top30 ASVs sorted by abundance across locations and selected when above 0.8 specialisation index as caluclated)

#Prep libraries and metadata regardless of which subset is being run
library(phyloseq)
library(ape)
library(ggplot2)
theme_set(theme_bw())

#metadata applies for both top and hyper specialists
treat_spec <- as.data.frame(read.csv("mapping_file_spec_index.csv", row.names=1, header=TRUE))

## Top specialists ---------------------------------------------------------

#asv_id for top specialists by abundance across all samples
otu_table_spec <- as.data.frame(read.csv("readmap_spec_index.csv", header=TRUE,row.names = "OTU_ID"))
#taxonomy for top specialists by abundance across all samples
taxmat_spec <- as.matrix(read.csv("tax_spec_index.csv", row.names=1, header=TRUE))
#use this one isntead of the line above now as I have appended higher taxonomy to levels that were missing to make figures clearer
taxmat_spec <- as.matrix(read.csv("tax_spec_index_append_hi_tax.csv", row.names=1, header=TRUE))

OTU = otu_table(otu_table_spec, taxa_are_rows = TRUE)
TAX = tax_table(taxmat_spec)
TREAT = sample_data(treat_spec)
TREE <- read.tree("rooted_tree.nwk")
OBJ_SPEC = phyloseq(OTU,TAX,TREAT,TREE)

library(magrittr)
OBJ_SPEC <- OBJ_SPEC %>%
  subset_taxa(
    Kingdom == "Bacteria" &
    Family  != "mitochondria" &
    Class   != "Chloroplast"
  )

OBJ1_spec_ts = transform_sample_counts(OBJ_SPEC, function(OTU) OTU/sum(OTU) )
OBJ1_spec_ts

#Top 30 ASVs above 0.8 specialisation index and high abundance
High_spec_gp <- subset_taxa(OBJ1_spec_ts, Kingdom=="Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:30]), High_spec_gp)
plot_heatmap(High_spec_gp, sample.label="Location") # this will do the default phyloseq ordination based sorting with ASV ID labels
plot_heatmap(High_spec_gp, sample.label="Location", taxa.label = "Species") # replace ASV labels with species
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.label = "Connection", taxa.label = "Species") #these arguments will let you reorganise order and labelling
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Species")
heatmap(otu_table(High_spec_gp))
#lable legibility
ASVheatplot_ORD <- plot_heatmap(High_spec_gp, sample.label="Location", taxa.label = "Species") # this will do the default phyloseq ordination based sorting
ASVheatplot_ORD + theme(axis.text.x = element_text(size=9, angle=80, hjust=0.4), axis.text.y=element_text(size=11))

#Lets try glom to group the asv's together
#will give a  smaller heatmap but might be an interesting/quick way of seeing if trends hold
#e.g will definitely combine all the individual Geobacter asv's but I'm unsure of the shared taxa between other ones

#Species
High_spec_SPE <- tax_glom(OBJ1_spec_ts,taxrank = "Species")
High_spec_SPE  <- subset_taxa(High_spec_SPE, Kingdom=="Bacteria")
High_spec_SPE  <- prune_taxa(names(sort(taxa_sums(High_spec_SPE),TRUE)[1:30]), High_spec_SPE)
plot_heatmap(High_spec_SPE, sample.label="Location", taxa.label = "Species") # this will do the default phyloseq ordination based sorting
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Species", taxa.label = "Species", weighted=TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Species", taxa.label = "Species", weighted=TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Species", weighted=TRUE)
heatmap(otu_table(High_spec_SPE))
#lets take out the individual species map ordered by the ordination for more legible font
speciesheatplot_ORD <- plot_heatmap(High_spec_SPE, sample.label="Location", taxa.label = "Species") # this will do the default phyloseq ordination based sorting
speciesheatplot_ORD + theme(axis.text.x = element_text(size=9, angle=80, hjust=0.4), axis.text.y=element_text(size=11))

#Genus
High_spec_GEN <- tax_glom(OBJ1_spec_ts,taxrank = "Genus")
High_spec_GEN  <- subset_taxa(High_spec_GEN, Kingdom=="Bacteria")
High_spec_GEN  <- prune_taxa(names(sort(taxa_sums(High_spec_GEN),TRUE)[1:30]), High_spec_GEN)
plot_heatmap(High_spec_GEN, sample.label="Location", taxa.label = "Genus") # this will do the default phyloseq ordination based sorting
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Genus", taxa.label = "Genus", weighted=TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Genus", taxa.label = "Genus", weighted=TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Genus", taxa.label = "Genus", weighted=TRUE)
heatmap(otu_table(High_spec_GEN))

#Family
High_spec_FAM <- tax_glom(OBJ1_spec_ts,taxrank = "Family")
High_spec_FAM  <- subset_taxa(High_spec_FAM, Kingdom=="Bacteria")
High_spec_FAM  <- prune_taxa(names(sort(taxa_sums(High_spec_FAM),TRUE)[1:30]), High_spec_FAM)
plot_heatmap(High_spec_FAM, sample.label="Location", taxa.label = "Family")
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Family", taxa.label = "Family", weighted=TRUE)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Family", taxa.label = "Family", weighted=TRUE)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted=TRUE, cex=2.5)
heatmap(otu_table(High_spec_FAM))
#test for labelling font size
heatplot <- plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted=TRUE)
heatplot + theme(axis.text.x = element_text(size=10, angle=90, hjust=0.4), axis.text.y=element_text(size=11))

#Order
High_spec_ORD <- tax_glom(OBJ1_spec_ts,taxrank = "Order")
High_spec_ORD  <- subset_taxa(High_spec_ORD, Kingdom=="Bacteria")
High_spec_ORD  <- prune_taxa(names(sort(taxa_sums(High_spec_ORD),TRUE)[1:30]), High_spec_ORD)
plot_heatmap(High_spec_ORD, sample.label="Location", taxa.label = "Order")
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted=TRUE)
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order", weighted=TRUE)
heatmap(otu_table(High_spec_ORD))

#Class
High_spec_CLA <- tax_glom(OBJ1_spec_ts,taxrank = "Class")
High_spec_CLA  <- subset_taxa(High_spec_CLA, Kingdom=="Bacteria")
High_spec_CLA  <- prune_taxa(names(sort(taxa_sums(High_spec_CLA),TRUE)[1:30]), High_spec_CLA)
plot_heatmap(High_spec_CLA, sample.label="Location", taxa.label = "Class")
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted=TRUE)
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Class", taxa.label = "Class", weighted=TRUE)
heatmap(otu_table(High_spec_CLA))

#Phyla
High_spec_PHY <- tax_glom(OBJ1_spec_ts,taxrank = "Phylum")
High_spec_PHY  <- subset_taxa(High_spec_PHY, Kingdom=="Bacteria")
High_spec_PHY  <- prune_taxa(names(sort(taxa_sums(High_spec_PHY),TRUE)[1:30]), High_spec_PHY)
plot_heatmap(High_spec_PHY, sample.label="Location", taxa.label = "Phylum")
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted=TRUE)
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Phylum", taxa.label = "Phylum", weighted=TRUE)
heatmap(otu_table(High_spec_PHY))

# "Hyper" specialists -----------------------------------------------------
# These are ASVs chose using the specialist index worksheet (as above). This timelooking specifically for ASVs that appear solely in one location, hence "hyper" specialists. 
# These are not sorted or picked by abundance, but instead by getting the # of samples that each ASV appearred in, and sorting for highest in one location, and lowest in the other two
# Excel allows the sorting for those 3 columns to best match the ascending/descending cimbination, and then ASVs with a sample count appearannce of 6 or greater were selected (initally was oging to go with # of 8, but the roots had zero at this and only 1 at 7)

#asv_id for "hyper" specialists
otu_table_hyper_spec <- as.data.frame(read.csv("readmap_hyper_spec_index.csv", header=TRUE,row.names = "OTU_ID"))
#taxonomy for "hyper" specialists
taxmat_hyper_spec <- as.matrix(read.csv("tax_hyper_spec_index.csv", row.names=1, header=TRUE))

OTU = otu_table(otu_table_hyper_spec, taxa_are_rows = TRUE)
TAX = tax_table(taxmat_hyper_spec)
TREAT = sample_data(treat_spec)
TREE <- read.tree("rooted_tree.nwk")
OBJ_HYPER_SPEC = phyloseq(OTU,TAX,TREAT,TREE)

# Note: doesn't seem to like making heatmaps from this particular TSS set....
# possibly due to the lack of overlap in presence of asvs (e.g entirely zeros for any two locations for each asv)
# So leaving on raw count....to be decided what to do here, maybe go back to rarefy?

#All ASVs
# Running on raw ASVs for now because it does not agree with taxa sorting otherwise, presumable cannot sort by taxa using distiance because there are too many zeros for a difference to be calculated
High_spec_gp <- subset_taxa(OBJ_HYPER_SPEC, Kingdom=="Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:140]), High_spec_gp)
plot_heatmap(High_spec_gp, sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
hyperheatplot <- plot_heatmap(High_spec_gp, sample.order = "Location", sample.label="Location", taxa.order = "Family", taxa.label = "Family")
hyperheatplot + theme(axis.text.x = element_text(size=8, angle=90, hjust=0.4), axis.text.y=element_text(size=7, angle = 0))
heatmap(otu_table(High_spec_gp))

#Top ASVs
High_spec_gp <- subset_taxa(OBJ_HYPER_SPEC, Kingdom=="Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:50]), High_spec_gp)
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
heatmap(otu_table(High_spec_gp))

#Lets try glom to group the asv's together
#will give a  smaller heatmap but might be an interesting/quick way of seeing if trends hold
#e.g will definitely combine all the individual Geobacter asv's but I'm unsure of the shared taxa between other ones

#Species
High_spec_SPE <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Species")
High_spec_SPE  <- subset_taxa(High_spec_SPE, Kingdom=="Bacteria")
High_spec_SPE  <- prune_taxa(names(sort(taxa_sums(High_spec_SPE),TRUE)[1:140]), High_spec_SPE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Species", taxa.label = "Species", weighted=TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Species", taxa.label = "Species", weighted=TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Genus", weighted=TRUE)
heatmap(otu_table(High_spec_SPE))

#Genus
High_spec_GEN <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Genus")
High_spec_GEN  <- subset_taxa(High_spec_GEN, Kingdom=="Bacteria")
High_spec_GEN  <- prune_taxa(names(sort(taxa_sums(High_spec_GEN),TRUE)[1:140]), High_spec_GEN)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Genus", taxa.label = "Genus", weighted=TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Genus", taxa.label = "Genus", weighted=TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Genus", taxa.label = "Genus", weighted=TRUE)
heatmap(otu_table(High_spec_GEN))

#Family
High_spec_FAM <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Family")
High_spec_FAM  <- subset_taxa(High_spec_FAM, Kingdom=="Bacteria")
High_spec_FAM  <- prune_taxa(names(sort(taxa_sums(High_spec_FAM),TRUE)[1:140]), High_spec_FAM)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted=TRUE)
heatmap(otu_table(High_spec_FAM))

#Order
High_spec_ORD <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Order")
High_spec_ORD  <- subset_taxa(High_spec_ORD, Kingdom=="Bacteria")
High_spec_ORD  <- prune_taxa(names(sort(taxa_sums(High_spec_ORD),TRUE)[1:140]), High_spec_ORD)
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order", weighted=TRUE)
hyperheatplot_ORD <- plot_heatmap(High_spec_ORD, sample.order = "Location", sample.label="Location", taxa.order = "Order", taxa.label = "Order")
hyperheatplot_ORD + theme(axis.text.x = element_text(size=8, angle=80, hjust=0.4), axis.text.y=element_text(size=8.5, angle = 0))
heatmap(otu_table(High_spec_gp))

heatmap(otu_table(High_spec_ORD))

#Class
High_spec_CLA <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Class")
High_spec_CLA  <- subset_taxa(High_spec_CLA, Kingdom=="Bacteria")
High_spec_CLA  <- prune_taxa(names(sort(taxa_sums(High_spec_CLA),TRUE)[1:30]), High_spec_CLA)
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Class", taxa.label = "Class", weighted=TRUE)
heatmap(otu_table(High_spec_CLA))

#Phyla
High_spec_PHY <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Phylum")
High_spec_PHY  <- subset_taxa(High_spec_PHY, Kingdom=="Bacteria")
High_spec_PHY  <- prune_taxa(names(sort(taxa_sums(High_spec_PHY),TRUE)[1:30]), High_spec_PHY)
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Phylum", taxa.label = "Phylum", weighted=TRUE)
heatmap(otu_table(High_spec_PHY))

# Beta diversity ----------------------------------------------------------



##BETA DIVERSITY_____________________________________________________________________________
##ORDINATIONS DENDROGRAMS ANOSIM PERMANOVA


## Ordinations -------------------------------------------------------------
##some ordinations first
library("ggplot2")
library("RColorBrewer")

##NMDS:
NMDS1 <- ordinate(OBJ1_exp_tss, "NMDS", distance="unifrac", weighted=TRUE, parallel=TRUE)
NMDS1 #use to check that stress is  < 0.2

#Subsetted timepoints WEIGHTED
NMDS_W0w <- ordinate(OBJ_W0_tss, "NMDS", distance="unifrac", weighted=TRUE, parallel=TRUE)
NMDS_W0w
NMDS_W6w <- ordinate(OBJ_W6_tss, "NMDS", distance="unifrac", weighted=TRUE, parallel=TRUE)
NMDS_W6w
NMDS_W8w <- ordinate(OBJ_W8_tss, "NMDS", distance="unifrac", weighted=TRUE, parallel=TRUE)
NMDS_W8w
NMDS_W14w <- ordinate(OBJ_W14_tss, "NMDS", distance="unifrac", weighted=TRUE, parallel=TRUE)
NMDS_W14w
#subsetted timepoints UNWEIGHTED
NMDS_W0u <- ordinate(OBJ_W0_tss, "NMDS", distance="unifrac", weighted=FALSE, parallel=TRUE)
NMDS_W0u
NMDS_W6u <- ordinate(OBJ_W6_tss, "NMDS", distance="unifrac", weighted=FALSE, parallel=TRUE)
NMDS_W6u
NMDS_W8u <- ordinate(OBJ_W8_tss, "NMDS", distance="unifrac", weighted=FALSE, parallel=TRUE)
NMDS_W8u
NMDS_W14u <- ordinate(OBJ_W14_tss, "NMDS", distance="unifrac", weighted=FALSE, parallel=TRUE)
NMDS_W14u

#Subsetted as only geobacter relatived taxa just out of curiosoty....
NMDS_Desulf <- ordinate(OBJ1_exp_tss, "NMDS", distance="unifrac", weighted=TRUE, parallel=TRUE)
NMDS_W0


## A Colour Palette Interlude ------------------------------------------------

#Section to play with themes and colours (now integrated below for plots, must move move/delete this later)
#- quick and dirty point and text size fix too (not sure if this will encounter problems down the line for shifting colvecs and mapp file ordering?)
p1u + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) + scale_color_manual(values = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
p1u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3) + scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#000000"))
p1u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3) + scale_color_brewer(palette = "Dark2")
p1u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3) + scale_color_brewer(palette = "Set1")
p1u + theme_dark() + theme(text = element_text(size = 14)) + geom_point(size = 3)

#### The palette with black that above colours in grey_theme plot were chosen from:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##Modifications to colour palette to account for showing all treatments in one ordination
#Dark and light colours of each to indicate connected/unconnected circuit
#e.g Geobacter connected: Dark Blue, Geobacter unconnected: Light Blue
#New colour palette for full dataset in one ordination is below:

#Geobacter Connected: Blue:"#177BB5"
#Geobacter Unconnected: Light Blue:"#56B4E9"
#Montebello Connected: Orange:"#BF8300"
#Montebello Unconnected: Light Orange:"#E09900"
#Pseudomonas Connected: Green:"#008F47"
#Pseudomonas Unconnected: Light Green:"#00B85C"
#Uninoculated Connected: Black ="#141414"
#Uninoculated Unconnected: Grey ="#7A7A7A"
my_colvec <- c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A")

#### To use for fills, add
bp + scale_fill_manual(values = cbp2)
### To use for line and point colors, add
sp + scale_colour_manual(values=cbp2)
#Themes
# https://ggplot2.tidyverse.org/reference/ggtheme.html
#Colours
#### https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/ , https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html


## Back to ordinations now -----------------------------------------------------
#Plotting the subsetted ordinations to look at individual time points, treatments, locations
#Week0
#Unweighted
p1u<-plot_ordination(OBJ_W0_tss, NMDS_W0u, color="Connection", shape="Location", label=NULL)
p1u<-plot_ordination(OBJ_W0_tss, NMDS_W0u, color="Inoculum", shape="Location", label=NULL)
p1u<-plot_ordination(OBJ_W0_tss, NMDS_W0u, color="Treatment", shape="Location", label=NULL)
p1u
p1u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p1w<-plot_ordination(OBJ_W0_tss, NMDS_W0w, color="Connection", shape="Location", label=NULL)
p1w<-plot_ordination(OBJ_W0_tss, NMDS_W0w, color="Inoculum", shape="Location", label=NULL)
p1w<-plot_ordination(OBJ_W0_tss, NMDS_W0w, color="Treatment", shape="Location", label=NULL)
p1w
p1w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Same plot with sample names added if needing to single out problematic sames etc
p1<-plot_ordination(OBJ_W0_tss, NMDS_W0, color="Inoculum", shape="Location", label="Sample_Name")
p1

#Week6
#unweighted
p2u<-plot_ordination(OBJ_W6_tss, NMDS_W6u, color="Connection", shape="Location", label=NULL)
p2u<-plot_ordination(OBJ_W6_tss, NMDS_W6u, color="Inoculum", shape="Location", label=NULL)
p2u<-plot_ordination(OBJ_W6_tss, NMDS_W6u, color="Treatment", shape="Location", label=NULL)
p2u
p2u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p2w<-plot_ordination(OBJ_W6_tss, NMDS_W6w, color="Connection", shape="Location", label=NULL)
p2w<-plot_ordination(OBJ_W6_tss, NMDS_W6w, color="Inoculum", shape="Location", label=NULL)
p2w<-plot_ordination(OBJ_W6_tss, NMDS_W6w, color="Treatment", shape="Location", label=NULL)
p2w
p2w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

p2<-plot_ordination(OBJ_W6_tss, NMDS_W6, color="Inoculum", shape="Location", label="Sample_Name")
p2

#Week8
#Unweighted
p3u<-plot_ordination(OBJ_W8_tss, NMDS_W8u, color="Connection", shape="Location", label=NULL)
p3u<-plot_ordination(OBJ_W8_tss, NMDS_W8u, color="Inoculum", shape="Location", label=NULL)
p3u<-plot_ordination(OBJ_W8_tss, NMDS_W8u, color="Treatment", shape="Location", label=NULL)
p3u
p3u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p3w<-plot_ordination(OBJ_W8_tss, NMDS_W8w, color="Connection", shape="Location", label=NULL)
p3w<-plot_ordination(OBJ_W8_tss, NMDS_W8w, color="Inoculum", shape="Location", label=NULL)
p3w<-plot_ordination(OBJ_W8_tss, NMDS_W8w, color="Treatment", shape="Location", label=NULL)
p3w
p3w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Week14
#Unweighted
p4u<-plot_ordination(OBJ_W14_tss, NMDS_W14u, color="Connection", shape="Location", label=NULL)
p4u<-plot_ordination(OBJ_W14_tss, NMDS_W14u, color="Inoculum", shape="Location", label=NULL)
p4u<-plot_ordination(OBJ_W14_tss, NMDS_W14u, color="Treatment", shape="Location", label=NULL)
p4u
p4u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#weighted
p4w<-plot_ordination(OBJ_W14_tss, NMDS_W14w, color="Connection", shape="Location", label=NULL)
p4w<-plot_ordination(OBJ_W14_tss, NMDS_W14w, color="Inoculum", shape="Location", label=NULL)
p4w<-plot_ordination(OBJ_W14_tss, NMDS_W14w, color="Treatment", shape="Location", label=NULL)
p4w
p4w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Same plot with sample names added if needing to single out problematic sames etc
p4<-plot_ordination(OBJ_W14_tss, NMDS_W14, color="Connection", shape="Location", label="Sample_Name")
p4

##PCoA:
PCoA1 <- ordinate(OBJ1_ts, "PCoA", distance="unifrac", weighted=TRUE, parallel=TRUE)
PCoA1$values #use to check eigen values if you wish

p1<-plot_ordination(OBJ1_ts, NMDS1, color="Treatment", shape="Location", label=NULL)
p1

p2<-plot_ordination(OBJ1_ts, NMDS1, color="Week", shape="Location", label=NULL)
p2

#work out the order of plot labels wiht in the facor we are using (i.e 'Group'):
levels(p1$data$Week)

colvec<-c("green4","red3","darkorange1","green4","red3","darkorange1") #make sure you colours appear in the order that corresponds to your factor levels
black <- rep("black",6) #colour vector for shape outlines

##now start again:
help(ordinate)

NMDS1 <- ordinate(OBJ1_ts, "PCoA", distance="unifrac", weighted=TRUE, parallel=TRUE)
p1<-plot_ordination(OBJ1_ts, NMDS1 , color="Week", shape="Location", label=NULL) + theme_bw()
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

## ANOSIM AND PERMANOVA ----------------------------

library(vegan)
#Setup objects for testing significance, effect size, correlation.
#First pull the variables you want to test into objects:
#By location at end
Location <- get_variable(OBJ_W14_tss, "Location")
Location <- get_variable(OBJ_W14_tss_GEN, "Location")
#By connection at end
Connection <- get_variable(OBJ_W14_tss, "Connection")
Connection <- get_variable(OBJ_W14_tss_FAM, "Connection")
#By Inoculum at end
Inoculum <- get_variable(OBJ_W14_tss, "Inoculum")
Inoculum <- get_variable(OBJ_W14_tss_ORD, "Inoculum")

#Quick test of W0 for these too for discussions sake if comes up as question
LocationW0 <- get_variable(OBJ_W0_tss, "Location")
#By connection at end
ConnectionW0 <- get_variable(OBJ_W0_tss, "Connection")
#By Inoculum at end
InoculumW0 <- get_variable(OBJ_W0_tss, "Inoculum")

##ANOSIMS(good for one-way beta diversity analysis)
##https://sites.google.com/site/mb3gustame/hypothesis-tests/anosim
##Paper explining ANOSIM:
##Non-parametric multivariate analyses of changes in community structure
##Australian Journal of Ecology. (1993) K. R. CLARKE

#Weighted
W14_Loc_ano_w <- anosim(distance(OBJ_W14_tss, "wunifrac"), Location)
W14_Loc_ano_w

#Unweighted
W14_Loc_ano_u <- anosim(distance(OBJ_W14_tss, "unifrac"), Location)
W14_Loc_ano_u

W14_Conn_ano <- anosim(distance(OBJ_W14_tss, "wunifrac"), Connection)
W14_Conn_ano

W14_Inoc_ano <- anosim(distance(OBJ_W14_tss, "wunifrac"), Inoculum)
W14_Inoc_ano

##PERMANOVA (function adonis())is similar to anosim but can handle more complex designs
##https://sites.google.com/site/mb3gustame/hypothesis-tests/manova/npmanova

##make distance matricies (these are weighted and unweighted unifrac. can also use 'bray')
OBJ1_W14perm_wu <- distance(OBJ_W14_tss, "wunifrac")
OBJ1_W14perm_u <- distance(OBJ_W14_tss, "unifrac")

#Taxa levels Weighted
OBJ1_W14perm_wu_G <- distance(OBJ_W14_tss_GEN, "wunifrac")
OBJ1_W14perm_wu_F <- distance(OBJ_W14_tss_FAM, "wunifrac")
OBJ1_W14perm_wu_O <- distance(OBJ_W14_tss_ORD, "wunifrac")

##Quick diversity test of W0 too (see above)
OBJ1_W0perm_wu <- distance(OBJ_W0_tss, "wunifrac")
OBJ1_W0perm_u <- distance(OBJ_W0_tss, "unifrac")

#Two ways, location by connection
W14_Group_ado = adonis(OBJ1_W14perm_wu ~ Location * Connection, permutations = 9999)
W14_Group_ado

#Three way comparison, should be able to just add in inoculum (order of factors shouldn't matter)
#Also added in permutations
#adonis Analysis of Variance, Perumatatioal
#NOTE must make factor for inoculum as above for location etc (three way version), *now fixed*
W14_Group_ado_w = adonis(OBJ1_W14perm_wu ~ Location * Connection * Inoculum, permutations = 9999)
W14_Group_ado_w

W14_Group_ado_u = adonis(OBJ1_W14perm_u ~ Location * Connection * Inoculum, permutations = 9999)
W14_Group_ado_u

#Quick test of W0 stats for discussions sake (see above)
W0_Group_ado_w = adonis(OBJ1_W0perm_wu ~ LocationW0 * ConnectionW0 * InoculumW0, permutations = 9999)
W0_Group_ado_w

W0_Group_ado_u = adonis(OBJ1_W0perm_u ~ LocationW0 * ConnectionW0 * InoculumW0, permutations = 9999)
W0_Group_ado_u

#Taxa Levels Weighted
#Order
W14_Group_ado_w_O = adonis(OBJ1_W14perm_wu_O ~ Location * Connection * Inoculum, permutations = 9999)
W14_Group_ado_w_O
#Family
W14_Group_ado_w_F = adonis(OBJ1_W14perm_wu_F ~ Location * Connection * Inoculum, permutations = 9999)
W14_Group_ado_w_F
#Genus
W14_Group_ado_w_G = adonis(OBJ1_W14perm_wu_G ~ Location * Connection * Inoculum, permutations = 9999)
W14_Group_ado_w_G

# DESeq2 -------------------------------------------------------------------

#html and bioclite may just be leftovers from old install processs, but keeping here in case they end up being needed
# BiocManager should handle it all though
install.packages("htmltools")
library(htmltools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

library(ggplot2)
library("DESeq2")

#Reminder if not already done: remove unwanted samples (the sampling control set, not relevant for logchange etc)
#i.e don't use the TSS transformed data
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")

#Phyloseq to DESEQ
diagdds = phyloseq_to_deseq2(OBJ1_exp, ~ Connection + Location) #Re: order of factors here, this would be testing for the effect of location, controlling for connection
#e.g for the above ~ Treatment + Location , the order of these matters, first one is what is controlled for and second one after the + is one is tested need to double check which is which
#run this on everything instead of just the specialists subset, because context of the full dataset is important for this differential abundance 

diagdds$Location<- relevel(diagdds$Location, ref="Root") # sets the reference point, baseline or control to be compared against
#look at the p adjusted value from the output NOT the regular p 

#Subset Week 14
diagdds <- diagdds[ , diagdds$Week == "Fourteen" ]
#check that subset was done as expected
as.data.frame( colData(diagdds) )

#Run model and factors
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
res # print out results

#for printing out results, alpha is a filter, dont run this if you want to get everything and sort them out yourself
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]

# WIP FOR SUBSETTING: Subset deseq results BY ASV to bypass need for manual creation of specialist OTU and TAX tables..
# Should be possible using the base R subset command (as opposed to subset_taxa) and specifiying your list of individual ASVs
# Functions fine for results themselves, but still working on fixing agglomerated taxa for the subset...may now be fixed!
#BPrepare W14 subsets
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")

#Bind taxononmy to results
res = cbind(as(res, "data.frame"), as(tax_table(OBJ_W14)[rownames(res), ], "matrix"))
res

#Aglommerate to certain taxa intead of ASV #I THINK THIS IS THE STUMBLING BLOCK, is overwriting deseq with the entire physeq objsect when we just want it to append taxa at the family level to match
res <- tax_glom(OBJ_W14,taxrank = "Family")
res
write.csv(as.data.frame(res), 
          file="DESeq2_specialist_subset_TEST.csv")

#This one will only ever work at ASV level, but may come in handy for any other work we do on the specilaist or "hyper" specilast subsets, lets us pull them out immediately
specialist_res_subset <- subset((res), rownames((res)) %in% c('0238e0e03ffd3faa629954545d336e61', '052f174cab37be599600e58c78283fa2', '0592d48e1457e0fafb90b4158fa522c2', '084e19e29dc9ebe03f4401003d10bff6', '26be318a519d7fe51cc6cf5d2378d1c8', '275c04dcabb809d184f0b6838763e20e', '2e7210652ae3b77f31c42743c148406b', '321da2e457e5a9b769814b25476c4b10', '33d9e0d38932e8ce14c359de566b7d08', '34c559c02664a1ac5ece941ca9000309', '4b8d75e30b64a18cf561c75cdc17043f', '4bb91812872f443514a3977be8c58774', '4ce53584fbaa2aa3650f10bbe615c714', '57e66fdc78f87cd026455c6394730932', '5fca9caffa56a57fcc31d7ccba92d008', '770af6feae23fae0ab3f1282a0ccbf18', '79eb38e43351aa3b12eae197935b81fb', '893c52ddb9dd678876c58f35b7ecbad6', '89514854a16cbb3269c2e9e94a05e9d7', '908664fbed1b4350a0be7c1dc38094a8', '9690acde73fcd49a81835468c4b92397', '9f9b3c564446dde56bbc0ef69261369f', 'a79e5838a5e0b50040e7f9aecf923028', 'b7f611ae7c6166d62354f04ca25391d8', 'ba37b62f122aca2aaff8ad84244df273', 'bd5b2a1bc31a73a4f37561a50e8e238c', 'c0664e6873e08cb34f7333015aabb07c', 'c36dba3bdc497773e638b8c441a9a0eb', 'c752096e70b99e9b3feeb36fb2beb2e1', 'd8a16afe0d36d2502e504377df9e7e44'))
write.csv(as.data.frame(specialist_res_subset), 
          file="DESeq2_specialist_subset_TEST.csv")
## END SUBSETTING WIP


#Different Comparison Direction Sheets
sigtabA = results(diagdds, contrast=c("Location","Root","Anode"))
sigtabB = results(diagdds, contrast=c("Location","Root","Cathode"))
sigtabC = results(diagdds, contrast=c("Location","Anode","Cathode"))
sigtabA # this appears as the same default res output (e.g strong responding Geobacter asv has same values, negative change, i.e it is log fold change AWAY from root, TO Anode)
sigtabB # following the above logic this should display change from Root to Cathode
sigtabC # and finally, this should display change from Anode to Cathode

#Bind results sheets with OTU Taxa
sigtab_otuA = cbind(as(sigtabA, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtabA), ], "matrix"))
sigtab_otuA

sigtab_otuB = cbind(as(sigtabB, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtabB), ], "matrix"))
sigtab_otuB

sigtab_otuC = cbind(as(sigtabC, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtabC), ], "matrix"))
sigtab_otuC

#Export .csv
write.csv(as.data.frame(sigtab_otuA), 
          file="DESeq2_resultsA.csv")
write.csv(as.data.frame(sigtab_otuB), 
          file="DESeq2_resultsB.csv")
write.csv(as.data.frame(sigtab_otuC), 
          file="DESeq2_resultsC.csv")

#convert to matrix
sigtab_otu = cbind(as(sigtab, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtab), ], "matrix"))

write.csv(as.data.frame(res), 
          file="raw_results.csv")

#only print out results with adjusted P higher than:
resSig <- subset(res, padj < 0.1)
resSig
write.csv(as.data.frame(resSig), 
          file="sig_results.csv")

resultsNames(diagdds)

# Figure of results
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
 
# Phylum order
x = tapply(sigtab_otu$log2FoldChange, sigtab_otu$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_otu$Phylum = factor(as.character(sigtab_otu$Phylum), levels=names(x))
# Class order
x = tapply(sigtab_otu$log2FoldChange, sigtab_otu$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_otu$Phylum = factor(as.character(sigtab_otu$Phylum), levels=names(x))
 
#Visulise
ggplot(sigtab_otu, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

####Jen's DESEQ2 script to compare, figure out adapting commands
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

# Heirarchical Clustering -------------------------------------------------



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


# Barcharts ---------------------------------------------------------------

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

# Deseq example ------------------------------------------------------------

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

