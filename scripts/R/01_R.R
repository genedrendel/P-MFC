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


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")

#For installing Phyloseq and associated packages, use Bioconductor / Biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

BiocManager::valid()

# Working Directory -------------------------------------------------------
#Set your working directory, location for data tables, tree, etc
setwd("~/Documents/University/Analysis/PMFC_18/2020 rerun outputs/Format for phyloseq")

#Quick import all to skip the below
library(phyloseq)
library(ape)
library(magrittr)
library(ggplot2)

otu_table <- as.data.frame(read.csv("raw_readmap.csv", header = TRUE,row.names = "OTU_ID"))
taxmat <- as.matrix(read.csv("tax_table.csv", row.names = 1, header = TRUE))
treat <- as.data.frame(read.csv("mapping_file.csv", row.names = 1, header = TRUE))
TREE <- read.tree("rooted_tree.nwk")
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
TREAT = sample_data(treat)
OBJ1 = phyloseq(OTU,TAX,TREAT,TREE)

#Calculate stats for dataset (before any subsetting)
OBJ1
sample_sums(OBJ1)
sort(sample_sums(OBJ1))
hist(sample_sums(OBJ1), main = "Histogram: Read Counts", xlab = "Total Reads", 
     border = "blue", col = "green", las = 1, breaks = 24)
TREAT$total_reads <- sample_sums(OBJ1)
TREAT$total_reads
print(mean(TREAT$total_reads))
print(var(TREAT$total_reads))
print(sd(TREAT$total_reads))


#Calculate stats for dataset (subsetted for just those being used in paper #1: wk0 and wk14)
OBJ1_paper <- subset_samples(OBJ1, Paper_Subset == "Included")
OBJ1_paper <- subset_samples(OBJ1, Treatment_Half_Trim == "Retain")

sample_sums(OBJ1_paper)
sort(sample_sums(OBJ1_paper))
hist(sample_sums(OBJ1_paper), main = "Histogram: Read Counts", xlab = "Total Reads", 
     border = "blue", col = "green", las = 1, breaks = 24)
total_reads <- sample_sums(OBJ1_paper)
total_reads
print(mean(total_reads))
print(var(total_reads))
print(sd(total_reads))


#continue processing for data anlysis
OBJ1 <- subset_taxa(OBJ1, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
OBJ1 <- OBJ1 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
    Order  != "Rickettsiales" &
    Family  != "mitochondria" &
    Class   != "Chloroplast"
  )
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")




#filter
OBJ1_exp = filter_taxa(OBJ1_exp, function(x) sum(x > 2) > (0.01*length(x)), TRUE)

######################################### test
#prevalence FILTERING??
# Subset to the remaining phyla
# Create table, number of features for each phyla
table(tax_table(OBJ1_exp)[, "Phylum"], exclude = NULL)
OBJ1_exp <- subset_taxa(OBJ1_exp, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(OBJ1_exp),
               MARGIN = ifelse(taxa_are_rows(OBJ1_exp), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(OBJ1_exp),
                    tax_table(OBJ1_exp))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Subset to the remaining phyla, line at 1% prevelance
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(OBJ1_exp, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(OBJ1_exp),color = Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.03, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position = "none")

# Define prevalence threshold as 1% of total samples
prevalenceThreshold = 0.01 * nsamples(OBJ1_exp)
prevalenceThreshold

keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
OBJ1_exp = prune_taxa(keepTaxa, OBJ1_exp)


# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(OBJ1_exp),
               MARGIN = ifelse(taxa_are_rows(OBJ1_exp), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(OBJ1_exp),
                    tax_table(OBJ1_exp))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})



############

OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
OBJ1_exp_tss = transform_sample_counts(OBJ1_exp, function(OTU) OTU/sum(OTU) )
#Half cut that retains Pseudomonas (use these ones moving forward from now!)
OBJ_Overall_TRIM <- subset_samples(OBJ1_exp, Treatment_Half_Trim == "Retain")
OBJ_Overall_TRIM_tss <- subset_samples(OBJ1_exp_tss, Treatment_Half_Trim == "Retain")
OBJ_W0 <- subset_samples(OBJ1_exp, Week == "Zero")
OBJ_W0_TRIM <- subset_samples(OBJ_W0, Treatment_Half_Trim == "Retain")
OBJ_W0_tss <- subset_samples(OBJ1_exp_tss, Week == "Zero")
OBJ_W0_TRIM_tss <- subset_samples(OBJ_W0_tss, Treatment_Half_Trim == "Retain")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")
OBJ_W14_TRIM <- subset_samples(OBJ_W14, Treatment_Half_Trim == "Retain")
OBJ_W14_tss <- subset_samples(OBJ1_exp_tss, Week == "Fourteen")
OBJ_W14_TRIM_tss <- subset_samples(OBJ_W14_tss, Treatment_Half_Trim == "Retain")

#Main Functional object  outputs from this:
OBJ_Overall_TRIM
OBJ_Overall_TRIM_tss
OBJ_W0_TRIM
OBJ_W0_TRIM_tss
OBJ_W14_TRIM
OBJ_W14_TRIM_tss

####TEMP SECTION FOR TESTING / IMPORVING TRANFOMATION (TSS vs CSS vs Anderson log)
install.packages("remotes")
remotes::install_github("vmikk/metagMisc")

library(metagMisc)
#metaGmisc TSS total sum scaling (this can also use any other method available within decostand) https://rdrr.io/rforge/vegan/man/decostand.html
OBJ1_exp_tss = phyloseq_standardize_otu_abundance(OBJ1_exp, method = "total")
#Anderson log, a modified version of log (x+1)
OBJ1_exp_anderson = physeq_transform_anderson_log(OBJ1_exp)
#CSS cumulative sum scaling
OBJ1_exp_css = phyloseq_transform_css(OBJ1_exp, norm = TRUE, log = TRUE)

OBJ_Overall_TRIM <- subset_samples(OBJ1_exp, Treatment_Half_Trim == "Retain")
OBJ_Overall_TRIM_tss <- subset_samples(OBJ1_exp_tss, Treatment_Half_Trim == "Retain")
OBJ_W0 <- subset_samples(OBJ1_exp, Week == "Zero")
OBJ_W0_TRIM <- subset_samples(OBJ_W0, Treatment_Half_Trim == "Retain")
OBJ_W0_tss <- subset_samples(OBJ1_exp_tss, Week == "Zero")
OBJ_W0_TRIM_tss <- subset_samples(OBJ_W0_tss, Treatment_Half_Trim == "Retain")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")
OBJ_W14_TRIM <- subset_samples(OBJ_W14, Treatment_Half_Trim == "Retain")
OBJ_W14_tss <- subset_samples(OBJ1_exp_tss, Week == "Fourteen")
OBJ_W14_TRIM_tss <- subset_samples(OBJ_W14_tss, Treatment_Half_Trim == "Retain")
#CSS
OBJ_Overall_TRIM_css <- subset_samples(OBJ1_exp_css, Treatment_Half_Trim == "Retain")
OBJ_W0_css <- subset_samples(OBJ1_exp_css, Week == "Zero")
OBJ_W0_TRIM_css <- subset_samples(OBJ_W0_css, Treatment_Half_Trim == "Retain")
OBJ_W14_css <- subset_samples(OBJ1_exp_css, Week == "Fourteen")
OBJ_W14_TRIM_css <- subset_samples(OBJ_W14_css, Treatment_Half_Trim == "Retain")
#Anderson log
OBJ_Overall_TRIM_anderson <- subset_samples(OBJ1_exp_anderson, Treatment_Half_Trim == "Retain")
OBJ_W0_anderson <- subset_samples(OBJ1_exp_anderson, Week == "Zero")
OBJ_W0_TRIM_anderson <- subset_samples(OBJ_W0_anderson, Treatment_Half_Trim == "Retain")
OBJ_W14_anderson <- subset_samples(OBJ1_exp_anderson, Week == "Fourteen")
OBJ_W14_TRIM_anderson <- subset_samples(OBJ_W14_anderson, Treatment_Half_Trim == "Retain")
########





####POTENTIALLY IGNORE ALL OF THIS...t be continued and investigated vs ANCOM-BC..
# spoke to Josh and he just uses TSS for everything....even picrust

#okay so if the story is that we need to use TSS + log x+1 for everything that was previously just TSS'd, let's do a test of the adnerson log
#from meta g misc pacakge with miscelaneous functions
##Using metagMiscpackage to do both TSS and anderson log
##metagmisc standardise function using "total" is the same as the old method of manually doing the TSS transform (as below)
#PHYLO_tss_manual = transform_sample_counts(Path_PHYLO, function(OTU) OTU/sum(OTU) )
#https://github.com/vmikk/metagMisc/
#https://rdrr.io/github/vmikk/metagMisc/
#anderson log
#Old version of the below commands
#Path_PHYLO_tss <- phyloseq_standardize_otu_abundance(Path_PHYLO, method = "total")
#Path_PHYLO_tss 
#Path_PHYLO_log <- phyloseq_standardize_otu_abundance(Path_PHYLO_tss, method = "log")
#Path_PHYLO_log
#devtools::install_github("vmikk/metagMisc")

install.packages("remotes")
remotes::install_github("vmikk/metagMisc")
library(metagMisc)
OBJ1_exp_tss = transform_sample_counts(OBJ1_exp, function(OTU) OTU/sum(OTU) )
OBJ1_exp_tss = phyloseq_standardize_otu_abundance(OBJ1_exp, method = "total")
OBJ1_exp_tss = physeq_transform_anderson_log(OBJ1_exp_tss)
OBJ1_exp_tss = physeq_transform_anderson_log(OBJ1_exp)
######################


#Deprecated, see above for half cut import workflow including overall, 0, 14, and raw vs tss
# CHOOSE ONE. -note after comparisons, unless anything else comes to light, should only need/use the half cut dataset from now on
#Plus one OR the other of the treatment cut versions (if ending up using them). New Treatment TSS'd subset
#Full cut , getting rid of both pseudo and montebello
OBJ_Overall_TRIMfull <- subset_samples(OBJ1_exp, Treatment_Trim == "Retain")
OBJ_W0_TRIMfull <- subset_samples(OBJ1_exp, Treatment_Trim == "Retain")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")
OBJ_W14_TRIMfull <- subset_samples(OBJ_W14, Treatment_Trim == "Retain")
# OR
#Half cut that retains Pseudomonas
OBJ_Overall_TRIMhalf <- subset_samples(OBJ1_exp, Treatment_Half_Trim == "Retain")
OBJ_W0_TRIMhalf <- subset_samples(OBJ1_exp, Treatment_Half_Trim == "Retain")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")
OBJ_W14_TRIMhalf <- subset_samples(OBJ_W14, Treatment_Half_Trim == "Retain")

### Import ------------------------------------------------------------------
#Import .csv as OTU table
#Replaced import with the below: Not sure if any functional difference once imported like this?
#But the other one didn't work, assuming mostly just a txt vs csv thing
otu_table <- as.data.frame(read.csv("raw_readmap.csv", header = TRUE,row.names = "OTU_ID"))

##look at the total reads per sample and decide on a rarefaction depth
rare <- colSums(otu_table)
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

taxmat <- as.matrix(read.csv("tax_table.csv", row.names = 1, header = TRUE))

##row names for taxmat and otu table MUST be the same
rownames(taxmat)
rownames(otu_table)

##Load and create a metadata data frame
##run make this in excel:first row should be your sample names, subsequent rows are treatment aspects
##IMPORTANT : you must have more than one column in your treatment file
##(remember the first column will be pulled out to be used as labels so more than 2 coloums is needed)

treat <- as.data.frame(read.csv("mapping_file.csv", row.names = 1, header = TRUE))

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

#### Remove unwanted taxa ----------------------------------------------------
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

#### Subsetting --------------------------------------------------------------

#First thing first GET RID OF THE CONTROL SAMPLE AND ONLY SUBSET FROM THAT OBJECT
#SO THAT IT STOPS FUCKING WITH ORDINATION DISTANCES
#Similarly, will use this to remove problematic samples , 
#e.g those with low read depth that entirely skew ordination patterns by being entirely alone
#To do this, created additional column of treatment file indicating Y (yes) and N (no) for all data minus the control
#Lets call this OBJ1_exp as an indication of it being the proper "experimental" dataset

#### Experimental Sample Subset (and TSS it) ----------------------------------------------
#Main Experimental Data Subset - run this one EVERY TIME
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
#TSS Transform on ALL Experimental data as pre-treatment for ordinations (rather than doing TSS for each individual subset)
#Overall TSS - then run the individual subsets, BUT ensure these are off of the TSS set, not the regular set in the subsetting section
OBJ1_exp_tss = transform_sample_counts(OBJ1_exp, function(OTU) OTU/sum(OTU) )
OBJ1_exp_tss

#Subsets on TSS data

#### Treatment Subsets -------------------------------------------------------
#Treatments (by connection) - for trying to get all the info on treatment combos in one
#(e.g keep in mind for ordinations = using open and closed symbols + colour
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

#New Treatment Cuts RAW COUNTS
#New Treatment TSS'd subsets

#Overall cut, full cut , getting rid of both pseudo and monebello
OBJ_Overall_TRIMfull <- subset_samples(OBJ1_exp, Treatment_Trim == "Retain")
#Just Week Zero
OBJ_W0_TRIMfull <- subset_samples(OBJ1_exp, Treatment_Trim == "Retain")
#Just Week 14
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")
OBJ_W14_TRIMfull <- subset_samples(OBJ_W14, Treatment_Trim == "Retain")
#OR
#Half trim that retains Pseudomonas
OBJ_Overall_TRIMhalf <- subset_samples(OBJ1_exp, Treatment_Half_Trim == "Retain")
#Just Week Zero
OBJ_W0 <- subset_samples(OBJ1_exp, Week == "Zero")
OBJ_W0_TRIMhalf <- subset_samples(OBJ1_exp, Treatment_Half_Trim == "Retain")
#Just Week 14
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")
OBJ_W14_TRIMhalf <- subset_samples(OBJ_W14, Treatment_Half_Trim == "Retain")

#New Treatment TSS'd subsets
#Overall - cut both pseudo and monebello
OBJ_Overall_TRIM_tss <- subset_samples(OBJ1_exp_tss, Treatment_Trim == "Retain")
#Half trim that retains Pseudomonas
OBJ_Overall_TRIM_tss <- subset_samples(OBJ1_exp_tss, Treatment_Half_Trim == "Retain")
#Just Week Zero
OBJ_W0_TRIM_tss <- subset_samples(OBJ1_W0_tss, Treatment_Trim == "Retain")
#Just Week 14
OBJ_W14_TRIM_tss <- subset_samples(OBJ_W14_tss, Treatment_Trim == "Retain")

#Wk14 Connected and Unconnected Subsets
OBJ1_W14_connected_tss <- subset_samples(OBJ_W14_tss, Connection == "Connected")
OBJ1_W14_unconnected_tss <- subset_samples(OBJ_W14_tss, Connection == "Unconnected")

#W14 locations (for including Root in comparisons)
OBJ1_14anode_tss <- subset_samples(OBJ_W14_tss, Location == "Anode")
OBJ1_14cathode_tss <- subset_samples(OBJ_W14_tss, Location == "Cathode")
OBJ1_14roots_tss <- subset_samples(OBJ_W14_tss, Location == "Root")

#Treatment
OBJ_Unin_tss <- subset_samples(OBJ1_exp_tss, Inoculum == "Uninoculated")
OBJ_Geo_tss <- subset_samples(OBJ1_exp_tss, Inoculum == "Geobacter")
OBJ_Mont_tss <- subset_samples(OBJ1_exp_tss, Inoculum == "Montebello")
OBJ_Pseudo_tss <- subset_samples(OBJ1_exp_tss, Inoculum == "Pseudomonas")

#### Non-transformed subsets ------------------------------------
#Subsets for NON transformed dataset
#Retained for other tests that use their own transformations , treatments etc

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

#W14 locations (for including Root in comparisons)
OBJ1_14anode <- subset_samples(OBJ_W14, Location == "Anode")
OBJ1_14cathode <- subset_samples(OBJ_W14, Location == "Cathode")
OBJ1_14roots <- subset_samples(OBJ_W14, Location == "Root")

sample_data(OBJ1_anode)

#### Agglomerate taxa -------------------------------------------------
#e.g subset for taxa too (ATTENTION RE:still need to change for tss vs norm if goin to do this for ordinations vs boxplots etc)
#Experimental
OBJ1_Desulf = subset_taxa(OBJ1_exp_tss, Order == "Desulfuromonadales")
OBJ1_Geobacteraceae = subset_taxa(OBJ1_exp_tss, Family == "Geobacteraceae")
OBJ1_Geobacter = subset_taxa(OBJ1_exp_tss, Genus == "Geobacter")

#Agglomerating at certain taxonomic levels
#Will be basing off of TSS'd Data objects, currently mainly interested in end point stats, so W14
OBJ_W14_tss_ORD <- tax_glom(OBJ_W14_tss,taxrank = "Order")
OBJ_W14_tss_FAM <- tax_glom(OBJ_W14_tss,taxrank = "Family")
OBJ_W14_tss_GEN <- tax_glom(OBJ_W14_tss,taxrank = "Genus")

#### Subset specific ASVs/OTUs -----------------------------------------------

## WIP
#This one will only ever work at ASV level of course, but can then go back and agglomerate the subset if you with, 
#and may come in handy for any other work we do on the specilaist or "hyper" specilast subsets, lets us pull them out immediately
#If you read the help for subset_taxa in phyloseq it states it is just a convenience wrapper for the base R subset function
#that allows for easy passing of a variable in the sample_data to subset the phyloseq object by. It also provides details of how it uses the base subset function.
#To subset your phyloseq object by a list of OTU ids (and there may be a better way of doing this)
#you can subset your otu_table using the base R subset function if the row names match your list of OTU ids, 
#then merge the resulting OTU table back with the other components on your phyloseq object, e.g:
#assuming you have a phyloseq object named 'physeq'
OBJ1_ASVsubset <- subset(otu_table(OBJ1), rownames(otu_table(OBJ1)) %in% c('0238e0e03ffd3faa629954545d336e61', '052f174cab37be599600e58c78283fa2', '0592d48e1457e0fafb90b4158fa522c2',
                                                                           '084e19e29dc9ebe03f4401003d10bff6', '26be318a519d7fe51cc6cf5d2378d1c8', '275c04dcabb809d184f0b6838763e20e',
                                                                           '2e7210652ae3b77f31c42743c148406b', '321da2e457e5a9b769814b25476c4b10', '33d9e0d38932e8ce14c359de566b7d08',
                                                                           '34c559c02664a1ac5ece941ca9000309', '4b8d75e30b64a18cf561c75cdc17043f', '4bb91812872f443514a3977be8c58774',
                                                                           '4ce53584fbaa2aa3650f10bbe615c714', '57e66fdc78f87cd026455c6394730932', '5fca9caffa56a57fcc31d7ccba92d008',
                                                                           '770af6feae23fae0ab3f1282a0ccbf18', '79eb38e43351aa3b12eae197935b81fb', '893c52ddb9dd678876c58f35b7ecbad6',
                                                                           '89514854a16cbb3269c2e9e94a05e9d7', '908664fbed1b4350a0be7c1dc38094a8', '9690acde73fcd49a81835468c4b92397',
                                                                           '9f9b3c564446dde56bbc0ef69261369f', 'a79e5838a5e0b50040e7f9aecf923028', 'b7f611ae7c6166d62354f04ca25391d8',
                                                                           'ba37b62f122aca2aaff8ad84244df273', 'bd5b2a1bc31a73a4f37561a50e8e238c', 'c0664e6873e08cb34f7333015aabb07c',
                                                                           'c36dba3bdc497773e638b8c441a9a0eb', 'c752096e70b99e9b3feeb36fb2beb2e1', 'd8a16afe0d36d2502e504377df9e7e44'))
new_OBJ1_ASVsubset <- merge_phyloseq(OBJ1_ASVsubset, tax_table(OBJ1), sample_data(OBJ1), ...)
#print out as csv to test export? Was easy enough with deseq results, but w/ phylo objects may need some extra steps to test between taxa and otu etc..
write.csv(as.data.frame(specialist_res_subset), 
          file = "IndivdualASV_subset_TEST.csv")
#This results in a new phyloseq object that is subset to the list of OTU ids you provided.
## WIP

#### Rarefy / TSS ---------------------------------------------
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
OBJ1_ts_rounded = transform_sample_counts(OBJ1_ts, function(OTU) round(OTU,digits = 0) )
OBJ1_ts_rounded

#Re: TSS proportions and integers
#Microbiome analysis says can multiply by 1,000,000 for easier interpretation...let's test..?
OBJ1_ts_multiplied = transform_sample_counts(OBJ1_ts, function(OTU) OTU*1000000 )
otu_table(OBJ1_ts_multiplied)
#Josh Suggested trying by 100 insetad....may as well try a few variations 1000 etc... but question remains what should really be done
#e.g are shannons and simpsons impacted ONLY by ratios, or do they need whole numbers? Need to test some comparisons and see what other people do

#### Sheet Exports (test) -----------------------------------------------------------

##Small interlude to test exporting the TSS's OTU by itself, to check the transform and agglomerating wksheets
#export TSS sheet
# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(OBJ1_exp_tss), "matrix")
# transpose if necessary
#only if transposing:# if(taxa_are_rows(OBJ1_exp_tss)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
write.csv(OTUdf,"TSS_objexp_transform.csv", row.names = TRUE)

##Same as above but for agglomerated taxa
#family
OBJ1_tss_FAM <- tax_glom(OBJ1_exp_tss,taxrank = "Family")
OTU2 = as(otu_table(OBJ1_tss_FAM), "matrix")
# transpose if necessary
#only if transposing:# if(taxa_are_rows(physeq1)){OTU2 <- t(OTU2)}
# Coerce to data.frame
OTUdf2 = as.data.frame(OTU2)
write.csv(OTUdf,"TSS_Family_transform.csv", row.names = TRUE)
#write corresopnding tax table
TAX2 = as(tax_table(OBJ1_tss_FAM), "matrix")
TAXdf2 = as.data.frame(TAX2)
write.csv(TAXdf2,"TSS_Family_tax.csv", row.names = TRUE)

#genus
OBJ1_tss_GEN <- tax_glom(OBJ1_exp_tss,taxrank = "Genus")
OTU3 = as(otu_table(OBJ1_tss_GEN), "matrix")
# transpose if necessary
#if(taxa_are_rows(physeq1)){OTU3 <- t(OTU3)}
# Coerce to data.frame
OTUdf3 = as.data.frame(OTU3)
write.csv(OTUdf3,"TSS_Genus_transform.csv", row.names = TRUE)
#write corresopnding tax table
TAX3 = as(tax_table(OBJ1_tss_GEN), "matrix")
TAXdf3 = as.data.frame(TAX3)
write.csv(TAXdf3,"TSS_Genus_tax.csv", row.names = TRUE)


##### Agglomerated taxa sheet export for some manual digging 1/12/22 ------------------
#agglomerated, taxa, family , export , csv, downstream
Family_Agglomerated <- tax_glom(OBJ_W14_TRIM,taxrank = "Family")
Family_Agglomerated
FamilyAgglomerated_OTUs = as(otu_table(Family_Agglomerated), "matrix")
FamGlomframe = as.data.frame(FamilyAgglomerated_OTUs)
write.csv(as.data.frame(FamGlomframe), 
          file = "Family_AglloomeratedOTUsfor_VisualisingCalcs.csv")
Fam_AgllomTax = as(tax_table(Family_Agglomerated), "matrix")
FamGlomTaxFrame = as.data.frame(Fam_AgllomTax)
write.csv(FamGlomTaxFrame,"Family_AgloomeratedTaxa_forvisualisingcalcs.csv", row.names = TRUE)


# Alpha Diversity ---------------------------------------------------------
##ALPHA DIVERSITY_____________________________________________________________________________
##BOXPLOTS, ANOVA, TUKEY TEST

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

## Microbiome - Core Taxa and Plots ADAPTED from Jaq , and added Core ASV  Venn plot -------------------------------------------------------------

##So...the workflow suggested for analysis of the core microbiome is to transform it into compositional
#For the above experiment I did my own manual relative abundance 

#In light of that for the purposes of this experiment I am going to take un-TSS data and transform it
#and go through the process and see what happens, what the bar charts look like

library(RColorBrewer)
library(ggsci)
library(viridis)
library(microbiome)
library(knitr) #if using the kable i.e. table function

## first - transform to compositional abundance (relative abundance - yes)

OBJ1_mbiocomp <- microbiome::transform(OBJ_W14, "compositional")

#calculates prevalences for all tax groups - gives you an idea of whats there
head(prevalence(OBJ1_mbiocomp, detection = 0.2/100, sort = TRUE))

#This makes a phyloseq obj with only those taxa > 0.2% relative abundance and in >50% of samples
OBJ1_core <- core(OBJ1_mbiocomp, detection = 0.2/100, prevalence = 50/100)
OBJ1_core

#SKIP THIS IS FILTERING
#Just to TEST Raw core without filtering with the below steps
OBJ1_ord2 <- tax_glom(OBJ1_core,taxrank = "Genus")
OBJ1_ord2

### - NEED TO CUT THIS DOWN BY A BIT FOR VISUALS AT GENUS - see filter_taxa

## making the bar charts 
OBJ1_ord <- tax_glom(OBJ1_core,taxrank = "Genus")
OBJ1_ord

## Filter taxa to make visuals better
#by adjusting the number in here adds or takes away the number of taxa. Is the cut off. 
OBJ1_ord2 <- filter_taxa(OBJ1_ord, function(x) mean(x) > 0.002, TRUE)
OBJ1_ord2 

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col = sample(col_vector, n))

#without black borders....
p <- plot_bar(OBJ1_ord2, "Location", fill = "Genus")
p <- p + scale_fill_manual(values = col_vector)
p <- p + geom_bar(stat = "identity", position = "stack")
p

##PLOT - all samples
p <- plot_bar(OBJ1_ord2, "Location", fill = "Genus")
p <- p + scale_fill_manual(values = col_vector)
p <- p + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + labs(y = "Relative abundance")
p <- p + geom_bar(stat = "identity", position = "stack")
p <- p + guides(fill = guide_legend(title = "Genus")) # only need at genus due to lots of fam. pull across
p

#testing colour palettes because don't know what palette colvec will be big enough and allow distinguishing
p <- p + scale_fill_manual(values = col_vector)
p <- p + scale_fill_brewer("dark2")
p2 + scale_fill_igv()
p <- p + scale_fill_viridis_d(direction = -1)

p <- plot_frequencies(sample_data(OBJ1_ord2), "Anode", "Cathode")
print(p)

#The figures came out well - the genus maybe needs trimming slightly as has lots of OTUs

#be interesting to compare this to what the output from experiment 1
### welp lol - checked both the outcome of my manual relative abundance in excel and this
## they are basically the same. The phyloseq OBJs look the same and the outcome on 
# barcharts look the same at the same detection and prevalence. 

#Start at making venn diagram for core members at w14

#NOTE THAT THIS IS INDENEPENDENT OF ABOVE PERVANLENCE SETTINGS, 
#e.g build from the original subset, NOT the already filtered ones

#install.packages("eulerr") # If not installed
devtools::install_github('microsud/microbiomeutilities')

library(RColorBrewer)
library(eulerr)
library(microbiome)
library(microbiomeutilities)

location_list <- unique(as.character(meta(OBJ1_mbiocomp)$Location))
print(location_list)

#Loop for identifying core taxa for each location
list_core <- c() # an empty object to store information

#To Print Taxa names (instead of asv ID)
# format names
OBJ_mbio_taxa <- format_to_besthit(OBJ1_mbiocomp)
# check names
taxa_names(OBJ_mbio_taxa)[1:5]

for (n in location_list) { # for each variable n in location
    #print(paste0("Identifying Core Taxa for ", n))
    
    Loc.sub <- subset_samples(OBJ_mbio_taxa, location_list == n) # Choose sample from location by n
    
    core_m <- core_members(Loc.sub, # loc.sub is phyloseq selected with only samples from g 
                           detection = 0.002, # 0.002 in atleast 50% samples 
                           prevalence = 0.50)
    print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each location.
    list_core[[n]] <- core_m # add to a list core taxa for each group.
    #print(list_core)
}
#JustID
print(list_core)

#Make Venn
# Specify colors and plot venn
# supplying colors in the order they appear in list_core
mycols <- c(Anode = "#eb5717", Cathode = "#352ea8", Root = "#0fa676") 
plot(venn(list_core),
     fills = mycols)

### Normality Tests ---------------------------------------------------------

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
ANOVA1 <- aov(Div$Shannon ~ Div$Location * Div$Connection)
summary(ANOVA1)

ANOVA1 <- aov(Div$Shannon ~ Div$Location)
summary(ANOVA1)

TUKEY <- TukeyHSD(ANOVA1,'Div$Location', conf.level = 0.95)
TUKEY

### WIP - Boxplot labels ----------------------------------------------------

##WORK IN PROGRESS - AUTOMATING A WAY TO PUT THE LETTERS ON THE BOXPLOTS
# library
library(multcompView)
# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
# Extract labels and factor levels from Tukey post-hoc
     Tukey.levels <- TUKEY[[variable]][,4]
     Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])

     #I need to put the labels in the same order as in the boxplot :
     Tukey.labels$treatment = rownames(Tukey.labels)
     Tukey.labels = Tukey.labels[order(Tukey.labels$treatment) , ]
     return(Tukey.labels)
     }

LABELS <- generate_label_df(TUKEY, "Div$Location")
LABELS

my_colors <- c("green4","red3","darkorange1","green4","red3","darkorange1") 
#make sure you colours appear in the order that corresponds to your factor levels

# Draw the basic boxplot
a <- boxplot(Div$Shannon ~ Div$Location , ylim = c(min(Div$Shannon), 1.1*max(Div$Shannon)) ,
             col = my_colors[as.numeric(LABELS[,1])] , ylab = "Shannon Diversity" , main = "")
# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )
text( c(1:nlevels(Div$Group)) , a$stats[nrow(a$stats),] + over , LABELS[,1]  , col = my_colors[as.numeric(LABELS[,1])] )

##End WORK IN PROGRESS

##example Kruskal wallis with post hoc test
kruskal.test(Shannon ~ Location, Div)

library(FSA)
PT = dunnTest(Shannon ~ Location, data = Div, method = "bh")
PT

#### Boxplots ----------------------------------------------------------------

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
s1 <- s1  + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
s1 <- s1 + stat_summary(geom = 'text', label = c("a","ab","a","b","b","ab"), fun.y = max, vjust = -1)
s1 <- s1 + scale_x_discrete(name = "Title", limits = lim, labels = lim)
s1 <- s1 + scale_y_continuous(name = "Shannon diversty",  breaks = seq(5, 6.5, 0.5), limits = c(5, 6.5))
s1 <- s1 + geom_boxplot(fill = fill, colour = line, alpha = 0.5)
s1

c1 <- ggplot(Div, aes(Group,Chao1)) +  geom_boxplot()
c1 <- c1 + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
#c1 <- c1 + stat_summary(geom = 'text', label = c("a","a","a","b","b","b"), fun.y = max, vjust = -1)
c1 <- c1 + scale_x_discrete(name = "Title", limits = lim, labels = lab)
c1 <- c1 + scale_y_continuous(name = "Chao1",  breaks = seq(0, 2000, 500), limits = c(0, 2000))
c1 <- c1 + geom_boxplot(fill = fill, colour = line, alpha = 0.5)
c1

p1 <- ggplot(Div, aes(Group, simpson)) +  geom_boxplot()
p1 <- p1 + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
#p1 <- p1 + stat_summary(geom = 'text', label = c("a","ab","a","b","b","ab"), fun.y = max, vjust = -1)
p1 <- p1 + scale_x_discrete(name = "Title", limits = lim, labels = lab)
p1 <- p1 + scale_y_continuous(name = "Simpson",  breaks = seq(0,0.5,0.1), limits = c(0,0.5))
p1 <- p1 + geom_boxplot(fill = fill, colour = line, alpha = 0.5)
p1

#make a panel image
library(Rmisc)
multiplot(s1,c1,p1, cols = 3)

# Heatmaps ----------------------------------------------------------------
##PHYLOSEQ HEATMAPS (WIP)_____________________________________________
##INPROGRESS

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())

#Have your phyloseq object already made, TSS'd

#Overall
gpt_all <- subset_taxa(OBJ1_exp_tss, Kingdom == "Bacteria")
gpt_all <- prune_taxa(names(sort(taxa_sums(gpt_all),TRUE)[1:30]), gpt_all)
plot_heatmap(gpt_all, sample.label = "Location")
plot_heatmap(gpt_all, "NMDS", "unifrac", "Connection", "Family", sample.label = "Location", weighted = FALSE)
plot_heatmap(gpt_all, "NMDS", "unifrac", "Connection", "Family", sample.label = "Location", weighted = TRUE)
heatmap(otu_table(gpt))

gpt_prot <- subset_taxa(OBJ1_exp_tss, Phylum == "Proteobacteria")
gpt_prot <- prune_taxa(names(sort(taxa_sums(gpt_prot),TRUE)[1:30]), gpt_prot)
plot_heatmap(gpt_prot)
plot_heatmap(gpt_prot, "NMDS", "unifrac", "Connection", "Family", weighted = FALSE)
plot_heatmap(gpt_prot, "NMDS", "unifrac", "Connection", "Family", weighted = TRUE)
plot_heatmap(gpt_prot)

#Subsetting...test to try to get labels legible and see time series
gpt_sub_an <- subset_taxa(OBJ1_anode_tss, Kingdom == "Bacteria")
gpt_sub_an <- prune_taxa(names(sort(taxa_sums(gpt_sub_an),TRUE)[1:30]), gpt_sub_an)
plot_heatmap(gpt_sub_an, sample.label = "Connection")
plot_heatmap(gpt_sub_an, method = "NMDS", distance = "unifrac", sample.label = "Connection", taxa.label = "Genus")
heatmap(otu_table(gpt_sub_an))

#Subsetting w14...test 
gpt_sub14 <- subset_taxa(OBJ_W14_tss, Kingdom == "Bacteria")
gpt_sub14 <- prune_taxa(names(sort(taxa_sums(gpt_sub14),TRUE)[1:30]), gpt_sub14)
plot_heatmap(gpt_sub14, method = "MDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.label = "Family", weighted = TRUE)
heatmap(otu_table(gpt_sub14))

#Try everything by time for those asvs
gpt_alltime <- subset_taxa(OBJ1_exp_tss, Kingdom == "Bacteria")
gpt_alltime <- prune_taxa(names(sort(taxa_sums(gpt_alltime),TRUE)[1:30]), gpt_alltime)
plot_heatmap(gpt_alltime, method = "MDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted = TRUE)
heatmap(otu_table(gpt_alltime))

#Now agglomerate to diff taxa levels instead of asv
OBJ_all_tss_FAM <- tax_glom(OBJ1_exp_tss,taxrank = "Family")
OBJ_all_tss_FAM <- subset_taxa(OBJ_all_tss_FAM, Kingdom == "Bacteria")
OBJ_all_tss_FAM <- prune_taxa(names(sort(taxa_sums(OBJ_all_tss_FAM),TRUE)[1:30]), OBJ_all_tss_FAM)
plot_heatmap(OBJ_all_tss_FAM, method = "MDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted = TRUE)
heatmap(otu_table(OBJ_all_tss_FAM))

## Heatmap for specialisation index -------------------------------------
#Separating out most variable asvs across trial for paper figure
#Import heatmap subset .csv's (in this case: top30 ASVs sorted by abundance across locations and selected when above 0.8 specialisation index as caluclated)

#Prep libraries and metadata regardless of which subset is being run
library(phyloseq)
library(ape)
library(ggplot2)
theme_set(theme_bw())

#metadata applies for both top and hyper specialists
treat_spec <- as.data.frame(read.csv("mapping_file_spec_index.csv", row.names = 1, header = TRUE))
treat <- as.data.frame(read.csv("mapping_file.csv", row.names = 1, header = TRUE))
TREE <- read.tree("rooted_tree.nwk")

### Top specialists ---------------------------------------------------------

# FINAL dataset - specialist list derived from set without montebello -------

### FINAL dataset - specialist list derived from set with montebello treatment cut out....
#will make a few different sheet options here....
#One will be the same as previously made specilaist sheets, just the top 30 by average proportion from specialisation nidex
#next is a top 100 asvs
#and finallly a sheet with ALL Goebacter ASVs from the final dataset...this is ~204 total and also consists of all Geobacteraceae and Desulforomanas (because Geobacter is the only representative in this set)
#asv_id for top specialists by abundance across all samples (noteusing same object names past here as with the full set, so make sure to nly run one orthe other)

#Top 30
otu_table_spec <- as.data.frame(read.csv("ReadMap_FINAL_SPEC.csv", header = TRUE,row.names = "OTU_ID"))
taxmat_spec <- as.matrix(read.csv("Tax_FINAL_SPEC.csv", row.names = 1, header = TRUE))

#Top100
otu_table_spec <- as.data.frame(read.csv("ReadMap_FinalT100_SPEC.csv", header = TRUE,row.names = "OTU_ID"))
taxmat_spec <- as.matrix(read.csv("Tax_FinalT100_SPEC.csv", row.names = 1, header = TRUE))

#Geobacter
otu_table_spec <- as.data.frame(read.csv("ReadMap_FINALGeoONLY.csv", header = TRUE,row.names = "OTU_ID"))
taxmat_spec <- as.matrix(read.csv("Tax_FINALGeoONLY.csv", row.names = 1, header = TRUE))

#not done for new cut...yet, use this one instead of the line above now as I have appended higher taxonomy to levels that were missing to make figures clearer
#taxmat_spec <- as.matrix(read.csv("tax_spec_index_append_hi_tax.csv", row.names = 1, header = TRUE))


#####Old spec import
#asv_id for top specialists by abundance across all samples
otu_table_spec <- as.data.frame(read.csv("readmap_spec_index.csv", header = TRUE,row.names = "OTU_ID"))
#taxonomy for top specialists by abundance across all samples
taxmat_spec <- as.matrix(read.csv("tax_spec_index.csv", row.names = 1, header = TRUE))
#use this one instead of the line above now as I have appended higher taxonomy to levels that were missing to make figures clearer
taxmat_spec <- as.matrix(read.csv("tax_spec_index_append_hi_tax.csv", row.names = 1, header = TRUE))

OTU = otu_table(otu_table_spec, taxa_are_rows = TRUE)
TAX = tax_table(taxmat_spec)
TREAT = sample_data(treat)
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
High_spec_gp <- subset_taxa(OBJ1_spec_ts, Kingdom == "Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:300]), High_spec_gp)
plot_heatmap(High_spec_gp, sample.label = "Location") # default phyloseq ordination based sorting with ASV ID labels
plot_heatmap(High_spec_gp, sample.label = "Location", taxa.label = "Genus") # replace ASV labels with species
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.label = "Connection", taxa.label = "Species") #these arguments let you reorganise order/labelling
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species")
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
plot_heatmap(High_spec_gp, sample.label = "Location", sample.order = "Location", taxa.order = "Species") # default phyloseq ordination based sorting with ASV ID labels
heatmap(otu_table(High_spec_gp))
#label legibility
ASVheatplot_ORD <- plot_heatmap(High_spec_gp, sample.label = "Location", taxa.label = "Species") # this will do the default phyloseq ordination based sorting
ASVheatplot_ORD <- plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Species")
ASVheatplot_ORD + theme(axis.text.x = element_text(size = 9, angle = 80, hjust = 0.4), axis.text.y = element_text(size = 11))

#Redo for figure (and a second asv one to compare with deseq manually)
plot_ORD <- plot_heatmap(High_spec_gp, sample.label = "Location", taxa.label = "Species", sample.order = "Location", taxa.order = "Species") # replace ASV labels with species
plot_ORD <- plot_heatmap(High_spec_gp, sample.label = "Location", sample.order = "Location", taxa.order = "Species") # replace ASV labels with species
plot_ORD + theme(axis.text.x = element_text(size = 9, angle = 80, hjust = 0.4), axis.text.y = element_text(size = 11))

#Lets try glom to group the asv's together
#will give a  smaller heatmap but might be an interesting/quick way of seeing if trends hold
#e.g will definitely combine all the individual Geobacter asv's but I'm unsure of the shared taxa between other ones

#Species
High_spec_SPE <- tax_glom(OBJ1_spec_ts,taxrank = "Species")
High_spec_SPE  <- subset_taxa(High_spec_SPE, Kingdom == "Bacteria")
High_spec_SPE  <- prune_taxa(names(sort(taxa_sums(High_spec_SPE),TRUE)[1:50]), High_spec_SPE)
plot_heatmap(High_spec_SPE, sample.label = "Location", taxa.label = "Species") # this will do the default phyloseq ordination based sorting
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
heatmap(otu_table(High_spec_SPE))
#lets take out the individual species map ordered by the ordination for more legible font
speciesheatplot_ORD <- plot_heatmap(High_spec_SPE, sample.label = "Location", taxa.label = "Species") # this will do the default phyloseq ordination based sorting
speciesheatplot_ORD + theme(axis.text.x = element_text(size = 9, angle = 80, hjust = 0.4), axis.text.y = element_text(size = 11))

#Genus
High_spec_GEN <- tax_glom(OBJ1_spec_ts,taxrank = "Genus")
High_spec_GEN  <- subset_taxa(High_spec_GEN, Kingdom == "Bacteria")
High_spec_GEN  <- prune_taxa(names(sort(taxa_sums(High_spec_GEN),TRUE)[1:50]), High_spec_GEN)
plot_heatmap(High_spec_GEN, sample.label = "Location", taxa.label = "Genus") # this will do the default phyloseq ordination based sorting
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
heatmap(otu_table(High_spec_GEN))

#Family
High_spec_FAM <- tax_glom(OBJ1_spec_ts,taxrank = "Family")
High_spec_FAM  <- subset_taxa(High_spec_FAM, Kingdom == "Bacteria")
High_spec_FAM  <- prune_taxa(names(sort(taxa_sums(High_spec_FAM),TRUE)[1:50]), High_spec_FAM)
plot_heatmap(High_spec_FAM, sample.label = "Location", taxa.label = "Family")
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Family", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Family", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted = TRUE, cex = 2.5)
heatmap(otu_table(High_spec_FAM))
#test for labelling font size
heatplot <- plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted = TRUE)
heatplot + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 0.4), axis.text.y = element_text(size = 11))

#Order
High_spec_ORD <- tax_glom(OBJ1_spec_ts,taxrank = "Order")
High_spec_ORD  <- subset_taxa(High_spec_ORD, Kingdom == "Bacteria")
High_spec_ORD  <- prune_taxa(names(sort(taxa_sums(High_spec_ORD),TRUE)[1:50]), High_spec_ORD)
plot_heatmap(High_spec_ORD, sample.label = "Location", taxa.label = "Order")
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order", weighted = TRUE)
heatmap(otu_table(High_spec_ORD))

#Class
High_spec_CLA <- tax_glom(OBJ1_spec_ts,taxrank = "Class")
High_spec_CLA  <- subset_taxa(High_spec_CLA, Kingdom == "Bacteria")
High_spec_CLA  <- prune_taxa(names(sort(taxa_sums(High_spec_CLA),TRUE)[1:50]), High_spec_CLA)
plot_heatmap(High_spec_CLA, sample.label = "Location", taxa.label = "Class")
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Class", taxa.label = "Class", weighted =  TRUE)
heatmap(otu_table(High_spec_CLA))

#Phyla
High_spec_PHY <- tax_glom(OBJ1_spec_ts,taxrank = "Phylum")
High_spec_PHY  <- subset_taxa(High_spec_PHY, Kingdom == "Bacteria")
High_spec_PHY  <- prune_taxa(names(sort(taxa_sums(High_spec_PHY),TRUE)[1:50]), High_spec_PHY)
plot_heatmap(High_spec_PHY, sample.label = "Location", taxa.label = "Phylum")
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Phylum", taxa.label = "Phylum", weighted = TRUE)
heatmap(otu_table(High_spec_PHY))

#### Top specialists (pt2 - treatment CUT) ---------------------------------------------------------

#asv_id for top specialists by abundance across all samples
otu_table_spec <- as.data.frame(read.csv("readmap_newcut_spec.csv", header = TRUE,row.names = "OTU_ID"))
#taxonomy for top specialists by abundance across all samples
taxmat_spec <- as.matrix(read.csv("tax_newcut_spec.csv", row.names = 1, header = TRUE))
treat_spec <- as.data.frame(read.csv("mapping_file_spec_index.csv", row.names = 1, header = TRUE))

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
High_spec_gp <- subset_taxa(OBJ1_spec_ts, Kingdom == "Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:50]), High_spec_gp)
plot_heatmap(High_spec_gp, sample.label = "Location") # this will do the default phyloseq ordination based sorting with ASV ID labels
plot_heatmap(High_spec_gp, sample.label = "Location", taxa.label = "Species") # replace ASV labels with species
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.label = "Connection", taxa.label = "Species") #these arguments will let you reorganise order and labelling
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species")
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Species")
heatmap(otu_table(High_spec_gp))
#label legibility
ASVheatplot_ORD <- plot_heatmap(High_spec_gp, sample.label = "Location", taxa.label = "Species") # this will do the default phyloseq ordination based sorting
ASVheatplot_ORD + theme(axis.text.x = element_text(size = 9, angle = 80, hjust = 0.4), axis.text.y = element_text(size = 11))

#Redo for figure (and a second asv one to compare with deseq manually)
plot_ORD <- plot_heatmap(High_spec_gp, sample.label = "Location", taxa.label = "Species", sample.order = "Location", taxa.order = "Species") # replace ASV labels with species
plot_ORD <- plot_heatmap(High_spec_gp, sample.label = "Location", sample.order = "Location", taxa.order = "Species") # replace ASV labels with species
plot_ORD + theme(axis.text.x = element_text(size = 9, angle = 80, hjust = 0.4), axis.text.y = element_text(size = 11))

#Lets try glom to group the asv's together
#will give a  smaller heatmap but might be an interesting/quick way of seeing if trends hold
#e.g will definitely combine all the individual Geobacter asv's but I'm unsure of the shared taxa between other ones

#Species
High_spec_SPE <- tax_glom(OBJ1_spec_ts,taxrank = "Species")
High_spec_SPE  <- subset_taxa(High_spec_SPE, Kingdom == "Bacteria")
High_spec_SPE  <- prune_taxa(names(sort(taxa_sums(High_spec_SPE),TRUE)[1:50]), High_spec_SPE)
plot_heatmap(High_spec_SPE, sample.label = "Location", taxa.label = "Species") # this will do the default phyloseq ordination based sorting
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
heatmap(otu_table(High_spec_SPE))
#lets take out the individual species map ordered by the ordination for more legible font
speciesheatplot_ORD <- plot_heatmap(High_spec_SPE, sample.label = "Location", taxa.label = "Species") # this will do the default phyloseq ordination based sorting
speciesheatplot_ORD + theme(axis.text.x = element_text(size = 9, angle = 80, hjust = 0.4), axis.text.y = element_text(size = 11))

#Genus
High_spec_GEN <- tax_glom(OBJ1_spec_ts,taxrank = "Genus")
High_spec_GEN  <- subset_taxa(High_spec_GEN, Kingdom == "Bacteria")
High_spec_GEN  <- prune_taxa(names(sort(taxa_sums(High_spec_GEN),TRUE)[1:50]), High_spec_GEN)
plot_heatmap(High_spec_GEN, sample.label = "Location", taxa.label = "Genus") # this will do the default phyloseq ordination based sorting
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
heatmap(otu_table(High_spec_GEN))

#Family
High_spec_FAM <- tax_glom(OBJ1_spec_ts,taxrank = "Family")
High_spec_FAM  <- subset_taxa(High_spec_FAM, Kingdom == "Bacteria")
High_spec_FAM  <- prune_taxa(names(sort(taxa_sums(High_spec_FAM),TRUE)[1:50]), High_spec_FAM)
plot_heatmap(High_spec_FAM, sample.label = "Location", taxa.label = "Family")
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Family", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Family", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted = TRUE, cex = 2.5)
heatmap(otu_table(High_spec_FAM))
#test for labelling font size
heatplot <- plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted = TRUE)
heatplot + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 0.4), axis.text.y = element_text(size = 11))

#Order
High_spec_ORD <- tax_glom(OBJ1_spec_ts,taxrank = "Order")
High_spec_ORD  <- subset_taxa(High_spec_ORD, Kingdom == "Bacteria")
High_spec_ORD  <- prune_taxa(names(sort(taxa_sums(High_spec_ORD),TRUE)[1:50]), High_spec_ORD)
plot_heatmap(High_spec_ORD, sample.label = "Location", taxa.label = "Order")
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order", weighted = TRUE)
heatmap(otu_table(High_spec_ORD))

#Class
High_spec_CLA <- tax_glom(OBJ1_spec_ts,taxrank = "Class")
High_spec_CLA  <- subset_taxa(High_spec_CLA, Kingdom == "Bacteria")
High_spec_CLA  <- prune_taxa(names(sort(taxa_sums(High_spec_CLA),TRUE)[1:50]), High_spec_CLA)
plot_heatmap(High_spec_CLA, sample.label = "Location", taxa.label = "Class")
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Class", taxa.label = "Class", weighted = TRUE)
heatmap(otu_table(High_spec_CLA))

#Phyla
High_spec_PHY <- tax_glom(OBJ1_spec_ts,taxrank = "Phylum")
High_spec_PHY  <- subset_taxa(High_spec_PHY, Kingdom == "Bacteria")
High_spec_PHY  <- prune_taxa(names(sort(taxa_sums(High_spec_PHY),TRUE)[1:50]), High_spec_PHY)
plot_heatmap(High_spec_PHY, sample.label = "Location", taxa.label = "Phylum")
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Time", sample.label = "Time", taxa.label = "Family", weighted = TRUE)
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Phylum", taxa.label = "Phylum", weighted = TRUE)
heatmap(otu_table(High_spec_PHY))

###### "Hyper" specialists -----------------------------------------------------
# These are ASVs chose using the specialist index worksheet (as above). This timelooking specifically for ASVs that appear solely in one location, hence "hyper" specialists. 
# These are not sorted or picked by abundance, but instead by getting the # of samples that each ASV appearred in, and sorting for highest in one location, and lowest in the other two
# Excel allows the sorting for those 3 columns to best match the ascending/descending cimbination, and then ASVs with a sample count appearannce of 6 or greater were selected
#(initally was oging to go with # of 8, but the roots had zero at this and only 1 at 7)
library(phyloseq)
library(ape)

#asv_id for "hyper" specialists
otu_table_hyper_spec <- as.data.frame(read.csv("readmap_hyper_spec_index.csv", header = TRUE,row.names = "OTU_ID"))
#taxonomy for "hyper" specialists
taxmat_hyper_spec <- as.matrix(read.csv("tax_hyper_spec_index.csv", row.names = 1, header = TRUE))
#metadata applies for both top and hyper specialists
treat_spec <- as.data.frame(read.csv("mapping_file_spec_index.csv", row.names = 1, header = TRUE))

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
High_spec_gp <- subset_taxa(OBJ_HYPER_SPEC, Kingdom == "Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:140]), High_spec_gp)
plot_heatmap(High_spec_gp, sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
hyperheatplot <- plot_heatmap(High_spec_gp, sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
hyperheatplot + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.4), axis.text.y = element_text(size = 7, angle = 0))
heatmap(otu_table(High_spec_gp))

#Top ASVs
High_spec_gp <- subset_taxa(OBJ_HYPER_SPEC, Kingdom == "Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:50]), High_spec_gp)
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
heatmap(otu_table(High_spec_gp))

#Lets try glom to group the asv's together
#will give a  smaller heatmap but might be an interesting/quick way of seeing if trends hold
#e.g will definitely combine all the individual Geobacter asv's but I'm unsure of the shared taxa between other ones

#Species
High_spec_SPE <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Species")
High_spec_SPE  <- subset_taxa(High_spec_SPE, Kingdom == "Bacteria")
High_spec_SPE  <- prune_taxa(names(sort(taxa_sums(High_spec_SPE),TRUE)[1:140]), High_spec_SPE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Genus", weighted = TRUE)
heatmap(otu_table(High_spec_SPE))

#Genus
High_spec_GEN <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Genus")
High_spec_GEN  <- subset_taxa(High_spec_GEN, Kingdom == "Bacteria")
High_spec_GEN  <- prune_taxa(names(sort(taxa_sums(High_spec_GEN),TRUE)[1:140]), High_spec_GEN)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
heatmap(otu_table(High_spec_GEN))

#Family
High_spec_FAM <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Family")
High_spec_FAM  <- subset_taxa(High_spec_FAM, Kingdom == "Bacteria")
High_spec_FAM  <- prune_taxa(names(sort(taxa_sums(High_spec_FAM),TRUE)[1:140]), High_spec_FAM)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted = TRUE)
heatmap(otu_table(High_spec_FAM))

#Order
High_spec_ORD <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Order")
High_spec_ORD  <- subset_taxa(High_spec_ORD, Kingdom == "Bacteria")
High_spec_ORD  <- prune_taxa(names(sort(taxa_sums(High_spec_ORD),TRUE)[1:140]), High_spec_ORD)
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order", weighted = TRUE)
hyperheatplot_ORD <- plot_heatmap(High_spec_ORD, sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order")
hyperheatplot_ORD + theme(axis.text.x = element_text(size = 8, angle = 80, hjust = 0.4), axis.text.y = element_text(size = 8.5, angle = 0))
heatmap(otu_table(High_spec_gp))

heatmap(otu_table(High_spec_ORD))

#Class
High_spec_CLA <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Class")
High_spec_CLA  <- subset_taxa(High_spec_CLA, Kingdom == "Bacteria")
High_spec_CLA  <- prune_taxa(names(sort(taxa_sums(High_spec_CLA),TRUE)[1:30]), High_spec_CLA)
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Class", taxa.label = "Class", weighted = TRUE)
heatmap(otu_table(High_spec_CLA))

#Phyla
High_spec_PHY <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Phylum")
High_spec_PHY  <- subset_taxa(High_spec_PHY, Kingdom == "Bacteria")
High_spec_PHY  <- prune_taxa(names(sort(taxa_sums(High_spec_PHY),TRUE)[1:30]), High_spec_PHY)
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Phylum", taxa.label = "Phylum", weighted = TRUE)
heatmap(otu_table(High_spec_PHY))

##### ABRIDGED Hyper meatmap --------------------------------------------------
#Paper figure version (abridged)
#Because of sheer number of cathode specific ASVs when looking at unique ASVs by location putting an upper cap of 20 on most abundant cathode asvs
#rest can appear in supplementary figure if necessary

#Trimmed unique ASV specialist sheet import
library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)

#asv_id for "hyper" specialists
otu_table_hyper_spec <- as.data.frame(read.csv("readmap_hyper_spec_index_cut1.csv", header = TRUE,row.names = "OTU_ID"))
#taxonomy for "hyper" specialists
taxmat_hyper_spec <- as.matrix(read.csv("tax_hyper_spec_index cut1.csv", row.names = 1, header = TRUE))
#metadata applies for both top and hyper specialists
treat_spec <- as.data.frame(read.csv("mapping_file_spec_index.csv", row.names = 1, header = TRUE))

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
High_spec_gp <- subset_taxa(OBJ_HYPER_SPEC, Kingdom == "Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:140]), High_spec_gp)
plot_heatmap(High_spec_gp, sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
hyperheatplot <- plot_heatmap(High_spec_gp, sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
hyperheatplot + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.4), axis.text.y = element_text(size = 9, angle = 0))
heatmap(otu_table(High_spec_gp))

#Top ASVs
High_spec_gp <- subset_taxa(OBJ_HYPER_SPEC, Kingdom == "Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:50]), High_spec_gp)
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
heatmap(otu_table(High_spec_gp))

#Lets try glom to group the asv's together
#will give a  smaller heatmap but might be an interesting/quick way of seeing if trends hold
#e.g will definitely combine all the individual Geobacter asv's but I'm unsure of the shared taxa between other ones

#Species
High_spec_SPE <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Species")
High_spec_SPE  <- subset_taxa(High_spec_SPE, Kingdom == "Bacteria")
High_spec_SPE  <- prune_taxa(names(sort(taxa_sums(High_spec_SPE),TRUE)[1:140]), High_spec_SPE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Genus", weighted = TRUE)
heatmap(otu_table(High_spec_SPE))

#Genus
High_spec_GEN <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Genus")
High_spec_GEN  <- subset_taxa(High_spec_GEN, Kingdom == "Bacteria")
High_spec_GEN  <- prune_taxa(names(sort(taxa_sums(High_spec_GEN),TRUE)[1:140]), High_spec_GEN)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
heatmap(otu_table(High_spec_GEN))

#Family
High_spec_FAM <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Family")
High_spec_FAM  <- subset_taxa(High_spec_FAM, Kingdom == "Bacteria")
High_spec_FAM  <- prune_taxa(names(sort(taxa_sums(High_spec_FAM),TRUE)[1:140]), High_spec_FAM)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted = TRUE)
heatmap(otu_table(High_spec_FAM))

#Order
High_spec_ORD <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Order")
High_spec_ORD  <- subset_taxa(High_spec_ORD, Kingdom == "Bacteria")
High_spec_ORD  <- prune_taxa(names(sort(taxa_sums(High_spec_ORD),TRUE)[1:140]), High_spec_ORD)
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order", weighted = TRUE)
hyperheatplot_ORD <- plot_heatmap(High_spec_ORD, sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order")
hyperheatplot_ORD + theme(axis.text.x = element_text(size = 8, angle = 80, hjust = 0.4), axis.text.y = element_text(size = 8.5, angle = 0))
heatmap(otu_table(High_spec_gp))

heatmap(otu_table(High_spec_ORD))

#Class
High_spec_CLA <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Class")
High_spec_CLA  <- subset_taxa(High_spec_CLA, Kingdom == "Bacteria")
High_spec_CLA  <- prune_taxa(names(sort(taxa_sums(High_spec_CLA),TRUE)[1:30]), High_spec_CLA)
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Class", taxa.label = "Class", weighted = TRUE)
heatmap(otu_table(High_spec_CLA))

#Phyla
High_spec_PHY <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Phylum")
High_spec_PHY  <- subset_taxa(High_spec_PHY, Kingdom == "Bacteria")
High_spec_PHY  <- prune_taxa(names(sort(taxa_sums(High_spec_PHY),TRUE)[1:30]), High_spec_PHY)
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Phylum", taxa.label = "Phylum", weighted = TRUE)
heatmap(otu_table(High_spec_PHY))

##### Unique specialists (pt2 - treatment cut data) -----------------------------------------------------
# These are ASVs chose using the specialist index worksheet (as above). This timelooking specifically for ASVs that appear solely in one location, hence "hyper" specialists. 
# These are not sorted or picked by abundance, but instead by getting the # of samples that each ASV appearred in, and sorting for highest in one location, and lowest in the other two
# Excel allows the sorting for those 3 columns to best match the ascending/descending cimbination, and then ASVs with a sample count appearannce of 6 or greater were selected
#(initally was oging to go with # of 8, but the roots had zero at this and only 1 at 7)
library(phyloseq)
library(ape)

#asv_id for "hyper" specialists
otu_table_hyper_spec <- as.data.frame(read.csv("readmap_newcut_uniques.csv", header = TRUE,row.names = "OTU_ID"))
#taxonomy for "hyper" specialists
taxmat_hyper_spec <- as.matrix(read.csv("tax_newcut_uniques.csv", row.names = 1, header = TRUE))
#metadata applies for both top and hyper specialists
treat_spec <- as.data.frame(read.csv("mapping_file_spec_index.csv", row.names = 1, header = TRUE))

#ALTERNATIVELY: Run this one with the slimmed down cathode list ()
#asv_id for "hyper" specialists
otu_table_hyper_spec <- as.data.frame(read.csv("readmap_newcut_uniques_trim.csv", header = TRUE,row.names = "OTU_ID"))
#taxonomy for "hyper" specialists
taxmat_hyper_spec <- as.matrix(read.csv("tax_newcut_uniques_trim.csv", row.names = 1, header = TRUE))
#metadata applies for both top and hyper specialists
treat_spec <- as.data.frame(read.csv("mapping_file_spec_index.csv", row.names = 1, header = TRUE))

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
High_spec_gp <- subset_taxa(OBJ_HYPER_SPEC, Kingdom == "Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:200]), High_spec_gp)
plot_heatmap(High_spec_gp, sample.order = "Location", sample.label = "Location", taxa.order = "Genus", taxa.label = "Genus")
hyperheatplot <- plot_heatmap(High_spec_gp, sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
hyperheatplot + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.4), axis.text.y = element_text(size = 7, angle = 0))
heatmap(otu_table(High_spec_gp))

#Top ASVs
High_spec_gp <- subset_taxa(OBJ_HYPER_SPEC, Kingdom == "Bacteria")
High_spec_gp <- prune_taxa(names(sort(taxa_sums(High_spec_gp),TRUE)[1:50]), High_spec_gp)
plot_heatmap(High_spec_gp, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family")
heatmap(otu_table(High_spec_gp))

#Lets try glom to group the asv's together
#will give a  smaller heatmap but might be an interesting/quick way of seeing if trends hold
#e.g will definitely combine all the individual Geobacter asv's but I'm unsure of the shared taxa between other ones

#Species
High_spec_SPE <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Species")
High_spec_SPE  <- subset_taxa(High_spec_SPE, Kingdom == "Bacteria")
High_spec_SPE  <- prune_taxa(names(sort(taxa_sums(High_spec_SPE),TRUE)[1:140]), High_spec_SPE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Species", taxa.label = "Species", weighted = TRUE)
plot_heatmap(High_spec_SPE, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Species", taxa.label = "Genus", weighted = TRUE)
heatmap(otu_table(High_spec_SPE))

#Genus
High_spec_GEN <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Genus")
High_spec_GEN  <- subset_taxa(High_spec_GEN, Kingdom == "Bacteria")
High_spec_GEN  <- prune_taxa(names(sort(taxa_sums(High_spec_GEN),TRUE)[1:140]), High_spec_GEN)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Connection", sample.label = "Connection", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Treatment", sample.label = "Treatment", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
plot_heatmap(High_spec_GEN, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Genus", taxa.label = "Genus", weighted = TRUE)
heatmap(otu_table(High_spec_GEN))

#Family
High_spec_FAM <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Family")
High_spec_FAM  <- subset_taxa(High_spec_FAM, Kingdom == "Bacteria")
High_spec_FAM  <- prune_taxa(names(sort(taxa_sums(High_spec_FAM),TRUE)[1:140]), High_spec_FAM)
plot_heatmap(High_spec_FAM, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Family", taxa.label = "Family", weighted = TRUE)
heatmap(otu_table(High_spec_FAM))

#Order
High_spec_ORD <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Order")
High_spec_ORD  <- subset_taxa(High_spec_ORD, Kingdom == "Bacteria")
High_spec_ORD  <- prune_taxa(names(sort(taxa_sums(High_spec_ORD),TRUE)[1:140]), High_spec_ORD)
plot_heatmap(High_spec_ORD, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order", weighted = TRUE)
hyperheatplot_ORD <- plot_heatmap(High_spec_ORD, sample.order = "Location", sample.label = "Location", taxa.order = "Order", taxa.label = "Order")
hyperheatplot_ORD + theme(axis.text.x = element_text(size = 8, angle = 80, hjust = 0.4), axis.text.y = element_text(size = 8.5, angle = 0))
heatmap(otu_table(High_spec_gp))

heatmap(otu_table(High_spec_ORD))

#Class
High_spec_CLA <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Class")
High_spec_CLA  <- subset_taxa(High_spec_CLA, Kingdom == "Bacteria")
High_spec_CLA  <- prune_taxa(names(sort(taxa_sums(High_spec_CLA),TRUE)[1:30]), High_spec_CLA)
plot_heatmap(High_spec_CLA, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Class", taxa.label = "Class", weighted = TRUE)
heatmap(otu_table(High_spec_CLA))

#Phyla
High_spec_PHY <- tax_glom(OBJ_HYPER_SPEC,taxrank = "Phylum")
High_spec_PHY  <- subset_taxa(High_spec_PHY, Kingdom == "Bacteria")
High_spec_PHY  <- prune_taxa(names(sort(taxa_sums(High_spec_PHY),TRUE)[1:30]), High_spec_PHY)
plot_heatmap(High_spec_PHY, method = "NMDS", distance = "unifrac", sample.order = "Location", sample.label = "Location", taxa.order = "Phylum", taxa.label = "Phylum", weighted = TRUE)
heatmap(otu_table(High_spec_PHY))

# Beta diversity ----------------------------------------------------------
##BETA DIVERSITY_____________________________________________________________________________
##ORDINATIONS DENDROGRAMS ANOSIM PERMANOVA
## Ordinations -------------------------------------------------------------
##some ordinations first
library("ggplot2")
library("RColorBrewer")

#NOTE: Added to mapping file additional columns with more groupings based on: 
#Group all samples per soil columns: Soil_Column
#Group based on whether soil based subtrate (root/anode) vs Water for cathode samples: Substrate

# Final Ordinations using half cut set ------------------------------------

#Main Functional object  outputs from this:, double check that these are correctly referenced due to copy pasting etc
OBJ_Overall_TRIM_tss
OBJ_W0_TRIM_tss
OBJ_W14_TRIM_tss

#Note re: Rstudio export, using SVG format and dimensions 1150 x 900 seems to result in a "close enough to" square ordination when accounting for legend size 

# Week 14 weighted TSS
NMDS_W14wTRIM <- ordinate(OBJ_W14_TRIM_tss, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W14wTRIM
# Week 14 unweighted TSS
NMDS_W14uTRIM <- ordinate(OBJ_W14_TRIM_tss, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W14uTRIM

#W14 with TRIMMED data (half set, so FINAL) TSS
#Unweighted
p4uTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14uTRIM, color = "Connection", shape = "Location", label = NULL)
p4uTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14uTRIM, color = "Inoculum", shape = "Location", label = NULL)
p4uTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14uTRIM, color = "Treatment", shape = "Location", label = NULL)
p4uTRIM
p4uTRIM + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p4wTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14wTRIM, color = "Connection", shape = "Location", label = NULL)
p4wTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14wTRIM, color = "Inoculum", shape = "Location", label = NULL)
p4wTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14wTRIM, color = "Treatment", shape = "Location", label = NULL)
p4wTRIM
p4wTRIM + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Week0
#weighted TSS
NMDS_W0wTRIM <- ordinate(OBJ_W0_TRIM_tss, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W0wTRIM
#unweighted TSS
NMDS_W0uTRIM <- ordinate(OBJ_W0_TRIM_tss, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W0uTRIM

#WWeek 0
#Unweighted TSS
p1u <- plot_ordination(OBJ_W0_TRIM_tss, NMDS_W0uTRIM, color = "Connection", shape = "Location", label = NULL)
p1u <- plot_ordination(OBJ_W0_TRIM_tss, NMDS_W0uTRIM, color = "Inoculum", shape = "Location", label = NULL)
p1u <- plot_ordination(OBJ_W0_TRIM_tss, NMDS_W0uTRIM, color = "Treatment", shape = "Location", label = NULL)
p1u
p1u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted TSS
p1w <- plot_ordination(OBJ_W0_TRIM_tss, NMDS_W0wTRIM, color = "Connection", shape = "Location", label = NULL)
p1w <- plot_ordination(OBJ_W0_TRIM_tss, NMDS_W0wTRIM, color = "Inoculum", shape = "Location", label = NULL)
p1w <- plot_ordination(OBJ_W0_TRIM_tss, NMDS_W0wTRIM, color = "Treatment", shape = "Location", label = NULL)
p1w
p1w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#New section to reformat the above ordinations for publication
#### Ordination Publication export -------------------------------------------

#with legend
#EXPORT: using SVG format and dimensions 1150 x 900 ....too big lets try smaller
#734x528 seems quite good

#when plotting without legend...need to test a different set of resolutions....prsumably a square 1:1 pixel ratio SHOULD be feasible in this case..

#Line below for trying to force the axis decimals to beahve the same between diff. plots....so far unseuccessful
# scale_x_continuous(breaks = 0.01) + scale_y_continuous(breaks = 0.1) 

#TESTING FOR TWEAKED THEME
#Week 0 unweighted plots (for supps)
p1u + geom_point(size = 3.5) + xlab("Axis 1") + ylab("Axis 2") + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#f0f0f0","#7A7A7A")) + theme(
  #remove legend (for matching sizes of plots and/or exporting with the aim of combining 2 plots with just one legend for both)
    legend.position = "none",
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot axes and panel background
  axis.text.x = element_text(colour = "black", size = 12),
  axis.text.y = element_text(colour = "black", size = 12),
  axis.ticks = element_line(colour = "black"),
  #axis.line = element_line(colour = "gray"),
  panel.border = element_rect(colour = "black", fill = NA),
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = 'white'),
  # Change legend
  legend.background = element_rect(fill = "white", color = NA),
  legend.key = element_rect(color = "white", fill = "white"),
  legend.title = element_text(color = "black"),
  legend.text = element_text(color = "black")
  )

#Week 0 - weighted 
p1w + geom_point(size = 3.5) + xlab("Axis 1") + ylab("Axis 2") + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#f0f0f0","#7A7A7A")) + theme(
  #remove legend (for matching sizes of plots and/or exporting with the aim of combining 2 plots with just one legend for both)
  legend.position = "none",
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot axes and panel background
  axis.text.x = element_text(colour = "black", size = 12),
  axis.text.y = element_text(colour = "black", size = 12),
  axis.ticks = element_line(colour = "black"),
  #axis.line = element_line(colour = "gray"),
  panel.border = element_rect(colour = "black", fill = NA),
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = 'white'),
  # Change legend
  legend.background = element_rect(fill = "white", color = NA),
  legend.key = element_rect(color = "white", fill = "white"),
  legend.title = element_text(color = "black"),
  legend.text = element_text(color = "black")
  )

#Week 14 - unweighted
p4uTRIM + geom_point(size = 3.5) + xlab("Axis 1") + ylab("Axis 2") + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#f0f0f0","#7A7A7A")) + theme(
  #remove legend (for matching sizes of plots and/or exporting with the aim of combining 2 plots with just one legend for both)
  legend.position = "none",
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot axes and panel background
  axis.text.x = element_text(colour = "black", size = 12),
  axis.text.y = element_text(colour = "black", size = 12),
  axis.ticks = element_line(colour = "black"),
  #axis.line = element_line(colour = "gray"),
  panel.border = element_rect(colour = "black", fill = NA),
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = 'white'),
  # Change legend
  legend.background = element_rect(fill = "white", color = NA),
  legend.key = element_rect(color = "white", fill = "white"),
  legend.title = element_text(color = "black"),
  legend.text = element_text(color = "black")
  )

#Week 14 weighted
p4wTRIM + geom_point(size = 3.5) + xlab("Axis 1") + ylab("Axis 2") + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#f0f0f0","#7A7A7A")) + theme(
  #remove legend (for matching sizes of plots and/or exporting with the aim of combining 2 plots with just one legend for both)
  legend.position = "none",
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot axes and panel background
  axis.text.x = element_text(colour = "black", size = 12),
  axis.text.y = element_text(colour = "black", size = 12),
  axis.ticks = element_line(colour = "black"),
  #axis.line = element_line(colour = "gray"),
  panel.border = element_rect(colour = "black", fill = NA),
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = 'white'),
  # Change legend 
  legend.background = element_rect(fill = "white", color = NA),
  legend.key = element_rect(color = "white", fill = "white"),
  legend.title = element_text(color = "black"),
  legend.text = element_text(color = "black")
  )


########EXPERIEMNT OUT OF CURISOITY - ordination for half cut set with raw non-TSS...just to compare
#HALF CUT SET
#half cut (replacing old full cut lines with this as will be moving forward with the half cu set, including Pseudomonas but cutting out Montebello)
#Note re: Rstudio export, using SVG format and dimensions 1150 x 900 results in square ordination when accounting for legend size 
#weighted RAW
NMDS_W14wTRIM <- ordinate(OBJ_W14_TRIMhalf, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W14wTRIM
#unweighted RAW
NMDS_W14uTRIM <- ordinate(OBJ_W14_TRIMhalf, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W14uTRIM

#W14 with TRIMMED data (half set, so FINAL) RAW
p4uTRIM <- plot_ordination(OBJ_W14_TRIMhalf, NMDS_W14uTRIM, color = "Connection", shape = "Location", label = NULL)
p4uTRIM <- plot_ordination(OBJ_W14_TRIMhalf, NMDS_W14uTRIM, color = "Inoculum", shape = "Location", label = NULL)
p4uTRIM <- plot_ordination(OBJ_W14_TRIMhalf, NMDS_W14uTRIM, color = "Treatment", shape = "Location", label = NULL)
p4uTRIM
p4uTRIM + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

p4wTRIM <- plot_ordination(OBJ_W14_TRIMhalf, NMDS_W14wTRIM, color = "Connection", shape = "Location", label = NULL)
p4wTRIM <- plot_ordination(OBJ_W14_TRIMhalf, NMDS_W14wTRIM, color = "Inoculum", shape = "Location", label = NULL)
p4wTRIM <- plot_ordination(OBJ_W14_TRIMhalf, NMDS_W14wTRIM, color = "Treatment", shape = "Location", label = NULL)
p4wTRIM
p4wTRIM + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Week0
#weighted RAW
NMDS_W0wTRIM <- ordinate(OBJ_W0_TRIMhalf, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W0wTRIM
#unweighted RAW
NMDS_W0uTRIM <- ordinate(OBJ_W0_TRIMhalf, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W0uTRIM

#Unweighted RAW
p1u <- plot_ordination(OBJ_W0_TRIMhalf, NMDS_W0uTRIM, color = "Connection", shape = "Location", label = NULL)
p1u <- plot_ordination(OBJ_W0_TRIMhalf, NMDS_W0uTRIM, color = "Inoculum", shape = "Location", label = NULL)
p1u <- plot_ordination(OBJ_W0_TRIMhalf, NMDS_W0uTRIM, color = "Treatment", shape = "Location", label = NULL)
p1u
p1u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted RAW
p1w <- plot_ordination(OBJ_W0_TRIMhalf, NMDS_W0wTRIM, color = "Connection", shape = "Location", label = NULL)
p1w <- plot_ordination(OBJ_W0_TRIMhalf, NMDS_W0wTRIM, color = "Inoculum", shape = "Location", label = NULL)
p1w <- plot_ordination(OBJ_W0_TRIMhalf, NMDS_W0wTRIM, color = "Treatment", shape = "Location", label = NULL)
p1w
p1w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))



##NMDS:
NMDS1 <- ordinate(OBJ1_exp_tss, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS1 #use to check that stress is  < 0.2
p1overall <- plot_ordination(OBJ1_exp_tss, NMDS1, color = "Treatment", shape = "Location", label = NULL)
p1overall <- plot_ordination(OBJ1_exp_tss, NMDS1, color = "Week", shape = "Location", label = NULL)
p1overall
p1overall + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Subsetted timepoints WEIGHTED
NMDS_W0w <- ordinate(OBJ_W0_tss, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W0w
NMDS_W6w <- ordinate(OBJ_W6_tss, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W6w
NMDS_W8w <- ordinate(OBJ_W8_tss, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W8w
NMDS_W14w <- ordinate(OBJ_W14_tss, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W14w

#subsetted timepoints UNWEIGHTED
NMDS_W0u <- ordinate(OBJ_W0_tss, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W0u
NMDS_W6u <- ordinate(OBJ_W6_tss, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W6u
NMDS_W8u <- ordinate(OBJ_W8_tss, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W8u
NMDS_W14u <- ordinate(OBJ_W14_tss, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W14u


#Subsetted as only geobacter relatived taxa just out of curiosoty....
NMDS_Desulf <- ordinate(OBJ1_exp_tss, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W0

#Plotting the subsetted ordinations to look at individual time points, treatments, locations
#Week0
#Unweighted
p1u <- plot_ordination(OBJ_W0_tss, NMDS_W0u, color = "Connection", shape = "Location", label = NULL)
p1u <- plot_ordination(OBJ_W0_tss, NMDS_W0u, color = "Inoculum", shape = "Location", label = NULL)
p1u <- plot_ordination(OBJ_W0_tss, NMDS_W0u, color = "Treatment", shape = "Location", label = NULL)
p1u
p1u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p1w <- plot_ordination(OBJ_W0_tss, NMDS_W0w, color = "Connection", shape = "Location", label = NULL)
p1w <- plot_ordination(OBJ_W0_tss, NMDS_W0w, color = "Inoculum", shape = "Location", label = NULL)
p1w <- plot_ordination(OBJ_W0_tss, NMDS_W0w, color = "Treatment", shape = "Location", label = NULL)
p1w
p1w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Same plot with sample names added if needing to single out problematic sames etc
p1 <- plot_ordination(OBJ_W0_tss, NMDS_W0, color = "Inoculum", shape = "Location", label = "Sample_Name")
p1

#Test dark theme options for ASM Poster
p1u <- plot_ordination(OBJ_W0_tss, NMDS_W0u, color = "Treatment", shape = "Location", label = NULL)
p1u
p1u + theme_dark() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

p1u + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Use geom for shapes or lines betweeen points...e.g geom_polygon(aes(fill=Location)) will be base don locations..or can leave blank brackets for default
#geom_path()
#geom_polygon()

#Week6
#unweighted
p2u <- plot_ordination(OBJ_W6_tss, NMDS_W6u, color = "Connection", shape = "Location", label = NULL)
p2u <- plot_ordination(OBJ_W6_tss, NMDS_W6u, color = "Inoculum", shape = "Location", label = NULL)
p2u <- plot_ordination(OBJ_W6_tss, NMDS_W6u, color = "Treatment", shape = "Location", label = NULL)
p2u
p2u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p2w <- plot_ordination(OBJ_W6_tss, NMDS_W6w, color = "Connection", shape = "Location", label = NULL)
p2w <- plot_ordination(OBJ_W6_tss, NMDS_W6w, color = "Inoculum", shape = "Location", label = NULL)
p2w <- plot_ordination(OBJ_W6_tss, NMDS_W6w, color = "Treatment", shape = "Location", label = NULL)
p2w
p2w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

p2 <- plot_ordination(OBJ_W6_tss, NMDS_W6, color = "Inoculum", shape = "Location", label = "Sample_Name")
p2

#Week8
#Unweighted
p3u <- plot_ordination(OBJ_W8_tss, NMDS_W8u, color = "Connection", shape = "Location", label = NULL)
p3u <- plot_ordination(OBJ_W8_tss, NMDS_W8u, color = "Inoculum", shape = "Location", label = NULL)
p3u <- plot_ordination(OBJ_W8_tss, NMDS_W8u, color = "Treatment", shape = "Location", label = NULL)
p3u
p3u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p3w <- plot_ordination(OBJ_W8_tss, NMDS_W8w, color = "Connection", shape = "Location", label = NULL)
p3w <- plot_ordination(OBJ_W8_tss, NMDS_W8w, color = "Inoculum", shape = "Location", label = NULL)
p3w <- plot_ordination(OBJ_W8_tss, NMDS_W8w, color = "Treatment", shape = "Location", label = NULL)
p3w
p3w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Week14
#Unweighted
p4u <- plot_ordination(OBJ_W14_tss, NMDS_W14u, color = "Connection", shape = "Location", label = NULL)
p4u <- plot_ordination(OBJ_W14_tss, NMDS_W14u, color = "Inoculum", shape = "Location", label = NULL)
p4u <- plot_ordination(OBJ_W14_tss, NMDS_W14u, color = "Treatment", shape = "Location", label = NULL)
p4u
p4u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#weighted
p4w <- plot_ordination(OBJ_W14_tss, NMDS_W14w, color = "Connection", shape = "Location", label = NULL)
p4w <- plot_ordination(OBJ_W14_tss, NMDS_W14w, color = "Inoculum", shape = "Location", label = NULL)
p4w <- plot_ordination(OBJ_W14_tss, NMDS_W14w, color = "Treatment", shape = "Location", label = NULL)
p4w
p4w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Same plot with sample names added if needing to single out problematic sames etc
p4 <- plot_ordination(OBJ_W14_tss, NMDS_W14, color = "Connection", shape = "Location", label = "Sample_Name")
p4

###
#Try ordinations of JUST connected, and JUST unconected at W14 to further tease apart interactions between connection and community
#Weighted
NMDS_W14C_w <- ordinate(OBJ1_W14_connected_tss , "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W14U_w <- ordinate(OBJ1_W14_unconnected_tss , "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
#Unweighted
NMDS_W14C_u <- ordinate(OBJ1_W14_connected_tss , "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W14U_u <- ordinate(OBJ1_W14_unconnected_tss , "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)

#Plot weighted
p14C_w <- plot_ordination(OBJ1_W14_connected_tss, NMDS_W14C_w, color = "Inoculum", shape = "Location", label = NULL)
p14U_w <- plot_ordination(OBJ1_W14_unconnected_tss, NMDS_W14U_w, color = "Inoculum", shape = "Location", label = NULL)
p14C_w
p14U_w

p14C_w
p14C_w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
p14U_w
p14U_w + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#Plot unweighted
p14C_u <- plot_ordination(OBJ1_W14_connected_tss, NMDS_W14C_u, color = "Inoculum", shape = "Location", label = NULL)
p14U_u <- plot_ordination(OBJ1_W14_unconnected_tss, NMDS_W14U_u, color = "Inoculum", shape = "Location", label = NULL)
p14C_u
p14U_u

p14C_u
p14C_u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
p14U_u
p14U_u + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

OBJ1_W14_connected_tss 
OBJ1_W14_unconnected_tss 

###### Cleaned up plot element controls (for poster format) ----------------------------------------

#TESTING FOR CUSTOM THEME - DARK
p1w + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#f0f0f0","#7A7A7A")) + theme(
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot axes and panel background
  axis.text.x = element_text(colour = "gray", size = 12),
  axis.text.y = element_text(colour = "gray", size = 12),
  axis.ticks = element_line(colour = "gray"),
  #axis.line = element_line(colour = "gray"),
  panel.border = element_rect(colour = "gray", fill = NA),
  plot.background = element_rect(fill = "black"),
  panel.background = element_rect(fill = 'black'),
  # Change legend (only drawback with black is that it doesn't seem to be an easy thing to make the black legend match the while verions...)
  legend.background = element_rect(fill = "black", color = NA),
  legend.key = element_rect(color = "gray", fill = "black"),
  legend.title = element_text(color = "white"),
  legend.text = element_text(color = "white")
  )

#Using the above dark custom dark theme format:
#Week 0 - weighted (run the relevant W0 section above first)
p1w + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#f0f0f0","#7A7A7A")) + theme(
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot axes and panel background
  axis.text.x = element_text(colour = "gray", size = 12),
  axis.text.y = element_text(colour = "gray", size = 12),
  axis.ticks = element_line(colour = "gray"),
  #axis.line = element_line(colour = "gray"),
  panel.border = element_rect(colour = "gray", fill = NA),
  plot.background = element_rect(fill = "black"),
  panel.background = element_rect(fill = 'black'),
  # Change legend (only drawback with black is that it doesn't seem to be an easy thing to make the black legend match the while verions...)
  legend.background = element_rect(fill = "black", color = NA),
  legend.key = element_rect(color = "gray", fill = "black"),
  legend.title = element_text(color = "white"),
  legend.text = element_text(color = "white")
  )
#Week 14 - weighted (run the relevant W14 section above first)
p4w + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900",
                                                             "#008F47","#00B85C","#f0f0f0","#7A7A7A")) + theme(
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot axes and panel background
  axis.text.x = element_text(colour = "gray", size = 12),
  axis.text.y = element_text(colour = "gray", size = 12),
  axis.ticks = element_line(colour = "gray"),
  #axis.line = element_line(colour = "gray"),
  panel.border = element_rect(colour = "gray", fill = NA),
  plot.background = element_rect(fill = "black"),
  panel.background = element_rect(fill = 'black'),
  # Change legend (only drawback with black is that it doesn't seem to be an easy thing to make the black legend match the while verions...)
  legend.background = element_rect(fill = "black", color = NA),
  legend.key = element_rect(color = "gray", fill = "black"),
  legend.title = element_text(color = "white"),
  legend.text = element_text(color = "white")
  )

#PERMANOVA FOR THESE
#By location at end
LocationC <- get_variable(OBJ1_W14_connected_tss , "Location")
LocationU <- get_variable(OBJ1_W14_unconnected_tss , "Location")
#By Inoculum at end
InoculumC <- get_variable(OBJ1_W14_connected_tss, "Inoculum")
InoculumU <- get_variable(OBJ1_W14_unconnected_tss, "Inoculum")
#Distance
W14permC_w <- distance(OBJ1_W14_connected_tss, "wunifrac")
W14permU_w <- distance(OBJ1_W14_unconnected_tss, "wunifrac")
W14permC_u <- distance(OBJ1_W14_connected_tss, "unifrac")
W14permU_u <- distance(OBJ1_W14_unconnected_tss, "unifrac")

#Two ways, location by connection
W14_Cw_ado = adonis(W14permC_w ~ LocationC * InoculumC, permutations = 9999)
W14_Cw_ado
##
W14_Uw_ado = adonis(W14permU_w ~ LocationU * InoculumU, permutations = 9999)
W14_Uw_ado

W14_Cu_ado = adonis(W14permC_u ~ LocationC * InoculumC, permutations = 9999)
W14_Cu_ado
W14_Uu_ado = adonis(W14permU_u ~ LocationU * InoculumU, permutations = 9999)
W14_Uu_ado

##PCoA:
PCoA1 <- ordinate(OBJ1_ts, "PCoA", distance = "unifrac", weighted = TRUE, parallel = TRUE)
PCoA1$values #use to check eigen values if you wish

p1 <- plot_ordination(OBJ1_ts, NMDS1, color = "Treatment", shape = "Location", label = NULL)
p1

p2 <- plot_ordination(OBJ1_ts, NMDS1, color = "Week", shape = "Location", label = NULL)
p2

#work out the order of plot labels wiht in the facor we are using (i.e 'Group'):
levels(p1$data$Week)

colvec <- c("green4","red3","darkorange1","green4","red3","darkorange1") #make sure you colours appear in the order that corresponds to your factor levels
black <- rep("black",6) #colour vector for shape outlines

##now start again:
help(ordinate)

NMDS1 <- ordinate(OBJ1_ts, "PCoA", distance = "unifrac", weighted = TRUE, parallel = TRUE)
p1 <- plot_ordination(OBJ1_ts, NMDS1 , color = "Week", shape = "Location", label = NULL) + theme_bw()
print(p1)
p1 <- p1 +  geom_point(aes(colour = factor(Group), fill = factor(Group), shape = factor(Location)), size = 3.5) + ggtitle(NULL) + theme(legend.position = "none")
p1 <- p1 + scale_shape_manual(values = c(21,24))
p1 <- p1 + scale_colour_manual(values = black)
p1 <- p1 + scale_fill_manual(values = colvec, breaks = c("Bulk_0", "Bulk_20", "Bulk_100","Rhizo_0", "Rhizo_20", "Rhizo_100"))
p1

NMDS2 <- ordinate(OBJ1_r, "PCoA", distance = "unifrac", weighted = FALSE, parallel = TRUE)
p2 <- plot_ordination(OBJ1_r, NMDS2 , color = "Group", shape = "Location", label = NULL) + theme_bw()
p2 <- p2 +  geom_point(aes(colour = factor(Group), fill = factor(Group), shape = factor(Location)), size = 3.5) + ggtitle(NULL) + theme(legend.position = "none")
p2 <- p2 + scale_shape_manual(values = c(21,24))
p2 <- p2 + scale_colour_manual(values = black)
p2 <- p2 + scale_fill_manual(values = colvec, breaks = c("Bulk_0", "Bulk_20", "Bulk_100","Rhizo_0", "Rhizo_20", "Rhizo_100"))
p2

#make a panel image
library(Rmisc)
multiplot(p1,p2, cols = 2)

# ANOSIM AND PERMANOVA ----------------------------

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

#Using TRIMMED W14 Data
#By location at end
Location <- get_variable(OBJ_W14_TRIM, "Location")
#By connection at end
Connection <- get_variable(OBJ_W14_TRIM, "Connection")
#By Inoculum at end
Inoculum <- get_variable(OBJ_W14_TRIM, "Inoculum")

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

W14_Conn_ano <- anosim(distance(OBJ_W14_tss, "unifrac"), Connection)
W14_Conn_ano

W14_Inoc_ano <- anosim(distance(OBJ_W14_tss, "unifrac"), Inoculum)
W14_Inoc_ano

#Using TRIMMED data
#Weighted
W14_Loc_ano_w <- anosim(distance(OBJ_W14_TRIM_tss, "wunifrac"), Location)
W14_Loc_ano_w

W14_Conn_ano_w <- anosim(distance(OBJ_W14_TRIM_tss, "wunifrac"), Connection)
W14_Conn_ano_w

W14_Inoc_ano_w <- anosim(distance(OBJ_W14_TRIM_tss, "wunifrac"), Inoculum)
W14_Inoc_ano_w

#Unweighted
W14_Loc_ano_u <- anosim(distance(OBJ_W14_TRIM_tss, "unifrac"), Location)
W14_Loc_ano_u

W14_Conn_ano_u <- anosim(distance(OBJ_W14_TRIM_tss, "unifrac"), Connection)
W14_Conn_ano_u

W14_Inoc_ano_u <- anosim(distance(OBJ_W14_TRIM_tss, "unifrac"), Inoculum)
W14_Inoc_ano_u

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

# FINAL DATA - TRIM PERMANOVA ----------------------------------------------------------
library(vegan)
#Cleaned up with FINAL cut dataset
##make distance matricies (these are weighted and unweighted unifrac. can also use 'bray')
OBJ1_W0perm_W <- distance(OBJ_W0_TRIM_tss, "wunifrac")
OBJ1_W0perm_U <- distance(OBJ_W0_TRIM_tss, "unifrac")
OBJ1_W14perm_W <- distance(OBJ_W14_TRIM_tss, "wunifrac")
OBJ1_W14perm_U <- distance(OBJ_W14_TRIM_tss, "unifrac")

#Set up variable objects
Location0 <- get_variable(OBJ_W0_TRIM_tss, "Location")
Location14 <- get_variable(OBJ_W14_TRIM_tss, "Location")
#By connection
Connection0 <- get_variable(OBJ_W0_TRIM_tss, "Connection")
Connection14 <- get_variable(OBJ_W14_TRIM_tss, "Connection")
#By Inoculum 
Inoculum0 <- get_variable(OBJ_W0_TRIM_tss, "Inoculum")
Inoculum14 <- get_variable(OBJ_W14_TRIM_tss, "Inoculum")

#Week 0
#W0_3Group_ado_U = adonis(OBJ1_W0perm_U ~ Location0 * Connection0 * Inoculum0, permutations = 9999)
#W0_3Group_ado_U
#compare adonis2
W0_3Group_ado2_U = adonis2(OBJ1_W0perm_U ~ Location0 * Connection0 * Inoculum0, permutations = 9999)
W0_3Group_ado2_U

#W0_3Group_ado_W = adonis(OBJ1_W0perm_W ~ Location0 * Connection0 * Inoculum0, permutations = 9999)
#W0_3Group_ado_W
#compare adonis2
W0_3Group_ado2_W = adonis2(OBJ1_W0perm_W ~ Location0 * Connection0 * Inoculum0, permutations = 9999)
W0_3Group_ado2_W

#W0 weighted adonis2 test without ordering (check if model as a whole is significant: model)
W0_3Group_ado2_W_model = adonis2(OBJ1_W0perm_W ~ Location0 * Connection0 * Inoculum0, permutations = 9999, by = NULL)
W0_3Group_ado2_W_model

#Week 14
#W14_3Group_ado_U = adonis(OBJ1_W14perm_U ~ Location14 * Connection14 * Inoculum14, permutations = 9999)
#W14_3Group_ado_U
#compare adonis2
W14_3Group_ado2_U = adonis2(OBJ1_W14perm_U ~ Location14 * Connection14 * Inoculum14, permutations = 9999)
W14_3Group_ado2_U

#W14_3Group_ado_W = adonis(OBJ1_W14perm_W ~ Location14 * Connection14 * Inoculum14, permutations = 9999)
#W14_3Group_ado_W
#compare adonis2
W14_3Group_ado2_W = adonis2(OBJ1_W14perm_W ~ Location14 * Connection14 * Inoculum14, permutations = 9999)
W14_3Group_ado2_W

#adonis2 test without ordering (check if model as a whole is significant: model)
W14_3Group_ado2_W_TEST = adonis2(OBJ1_W14perm_W ~ Location14 * Connection14 * Inoculum14, permutations = 9999, by = NULL)
W14_3Group_ado2_W_TEST

#adonis2 test without ordering (when don't want order to matter: margin)
W14_3Group_ado2_W_TEST = adonis2(OBJ1_W14perm_W ~ Location14 * Connection14 * Inoculum14, permutations = 9999, by = "margin")
W14_3Group_ado2_W_TEST


#### Older run of PERMANOVA for TRIMMMED DATA for reference.... may have inconsistencyes re: raw vs tss 
##make distance matricies (these are weighted and unweighted unifrac. can also use 'bray')

OBJ1_W14perm_wu <- distance(OBJ_W14_TRIM, "wunifrac")
OBJ1_W14perm_u <- distance(OBJ_W14_TRIM, "unifrac")

#Two ways, location by connection
W14_Group_adoW = adonis(OBJ1_W14perm_wu ~ Location * Connection, permutations = 9999)
W14_Group_adoW

W14_Group_adoU = adonis(OBJ1_W14perm_u ~ Location * Connection, permutations = 9999)
W14_Group_adoU

#Three way comparison, should be able to just add in inoculum (order of factors shouldn't matter)
#Also added in permutations
#adonis Analysis of Variance, Perumatatioal
#NOTE must make factor for inoculum as above for location etc (three way version), *now fixed*
W14_3Group_ado_w = adonis(OBJ1_W14perm_wu ~ Location * Connection * Inoculum, permutations = 9999)
W14_3Group_ado_w

W14_3Group_ado_u = adonis(OBJ1_W14perm_u ~ Location * Connection * Inoculum, permutations = 9999)
W14_3Group_ado_u

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

#html and bioclite may just be leftovers from old install process, but keeping here in case they end up being needed
# BiocManager should handle it all though
install.packages("htmltools")
library(htmltools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

library(ggplot2)
library(DESeq2)

#Reminder if not already done: remove unwanted samples (the sampling control set, not relevant for logchange etc)
#i.e don't use the TSS transformed data
#Option 1
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
#Option 2
#Alternative subset for TRIM/TREATMENT CUT DATASSET
OBJ1_exp <- subset_samples(OBJ1, Treatment_Trim == "Retain")
#Option 3
#Alternative dataset with half cut (retaining Pseudo, but still cutting Montebello)
OBJ1_exp <- subset_samples(OBJ1, Treatment_Half_Trim == "Retain")


#Phyloseq to DESEQ for testing location differences (accounting for connection)
diagdds = phyloseq_to_deseq2(OBJ1_exp, ~ Connection + Location) #Re: order of factors here, this would be testing for the effect of location, controlling for connection
#e.g for the above ~ Connection + Location , the order of these matters, first one is what is controlled for and second one after the + is one is tested need to double check which is which
#Controlling for the presence of absence of connection, what was the effect of location
#run this on everything instead of just the specialists subset, because context of the full dataset is important for this differential abundance 

#Use Option 1, or 2, not both because they override,, just testing which comparisonsmight make a difference, this sets the baseline/control
diagdds$Location <- relevel(diagdds$Location, ref = "Root") # sets the reference point, baseline or control to be compared against

#option two...to test, because it seems some combinations are not tested as valid compairions when useing "name" to compare
diagdds$Location <- relevel(diagdds$Location, ref = "Anode") # sets the reference point, baseline or control to be compared against

#look at the p adjusted value from the output NOT the regular p 

#Subset Week 14
diagdds <- diagdds[ , diagdds$Week == "Fourteen" ]
#check that subset was done as expected
as.data.frame( colData(diagdds) )

#Run model and factors
diagdds = DESeq(diagdds, test = "Wald", fitType = "parametric")
res = results(diagdds, cooksCutoff = FALSE)
res # print out results

#Different Comparison Direction Sheets
sigtabR_A = results(diagdds, contrast = c("Location","Root","Anode")) #Direction of logfold change is shown as change from anode TO root (so for e.g Geobacters are all negative/decreases)
sigtabA_R = results(diagdds, contrast = c("Location","Anode","Root")) #Increase in Anode FROM Root. Test to see if this does actually switch direction of logfold change.
#Should make logfold of Geobacteris positive as they are increasing from root to anode
sigtabR_C = results(diagdds, contrast = c("Location","Root","Cathode")) #Should be to Root from Cathode (Positive # = Increase in Root FROM Cathode)
sigtabA_C = results(diagdds, contrast = c("Location","Anode","Cathode")) #Should be to Anode from Cathode (Postive # = Increase in Anode FROM Cathode)
sigtabR_A
sigtabA_R
sigtabR_C
sigtabA_C

#Bind results sheets with OTU Taxa
sigtab_taxR_A = cbind(as(sigtabR_A, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtabR_A), ], "matrix"))
sigtab_taxR_A
#for printing out results, below X significance  can skip this if you want to get everything and sort them out yourself
sigtab_taxR_A <- subset(sigtab_taxR_A, sigtab_taxR_A$padj < 0.2)
write.csv(as.data.frame(sigtab_taxR_A), 
          file = "DESeq2_Roots_relativeto_Anode.csv")

sigtab_taxA_R = cbind(as(sigtabA_R, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtabA_R), ], "matrix"))
sigtab_taxA_R
#for printing out results, below X significance  can skip this if you want to get everything and sort them out yourself
sigtab_taxA_R <- subset(sigtab_taxA_R, sigtab_taxA_R$padj < 0.2)
write.csv(as.data.frame(sigtab_taxA_R), 
          file = "DESeq2_Anode_relativeto_Roots.csv")

sigtab_taxR_C = cbind(as(sigtabR_C, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtabR_C), ], "matrix"))
sigtab_taxR_C
#for printing out results, below X significance  can skip this if you want to get everything and sort them out yourself
sigtab_taxR_C <- subset(sigtab_taxR_C, sigtab_taxR_C$padj < 0.2)
write.csv(as.data.frame(sigtab_taxR_C), 
          file = "DESeq2_Roots_relativeto_Cathode.csv")

sigtab_taxA_C = cbind(as(sigtabA_C, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtabA_C), ], "matrix"))
sigtab_taxA_C
#for printing out results, below X significance  can skip this if you want to get everything and sort them out yourself
sigtab_taxA_C <- subset(sigtab_taxA_C, sigtab_taxA_C$padj < 0.2)
write.csv(as.data.frame(sigtab_taxA_C), 
          file = "DESeq2_Anode_relativeto_Cathode.csv")

#Test - Cut out only specialsits (currently only ASV specialists from the FULL dataset) from each sheet (when jsut looking at the top 30 specialists) Otherwise export whole list (and/or with cutoff)
specialist_res_subsetR_A <- subset((sigtab_taxR_A), rownames((sigtab_taxR_A)) %in% c('0238e0e03ffd3faa629954545d336e61', '052f174cab37be599600e58c78283fa2', '0592d48e1457e0fafb90b4158fa522c2',
                                                                                     '084e19e29dc9ebe03f4401003d10bff6', '26be318a519d7fe51cc6cf5d2378d1c8', '275c04dcabb809d184f0b6838763e20e',
                                                                                     '2e7210652ae3b77f31c42743c148406b', '321da2e457e5a9b769814b25476c4b10', '33d9e0d38932e8ce14c359de566b7d08',
                                                                                     '34c559c02664a1ac5ece941ca9000309', '4b8d75e30b64a18cf561c75cdc17043f', '4bb91812872f443514a3977be8c58774',
                                                                                     '4ce53584fbaa2aa3650f10bbe615c714', '57e66fdc78f87cd026455c6394730932', '5fca9caffa56a57fcc31d7ccba92d008',
                                                                                     '770af6feae23fae0ab3f1282a0ccbf18', '79eb38e43351aa3b12eae197935b81fb', '893c52ddb9dd678876c58f35b7ecbad6',
                                                                                     '89514854a16cbb3269c2e9e94a05e9d7', '908664fbed1b4350a0be7c1dc38094a8', '9690acde73fcd49a81835468c4b92397',
                                                                                     '9f9b3c564446dde56bbc0ef69261369f', 'a79e5838a5e0b50040e7f9aecf923028', 'b7f611ae7c6166d62354f04ca25391d8',
                                                                                     'ba37b62f122aca2aaff8ad84244df273', 'bd5b2a1bc31a73a4f37561a50e8e238c', 'c0664e6873e08cb34f7333015aabb07c',
                                                                                     'c36dba3bdc497773e638b8c441a9a0eb', 'c752096e70b99e9b3feeb36fb2beb2e1', 'd8a16afe0d36d2502e504377df9e7e44'))
specialist_res_subsetA_R <- subset((sigtab_taxA_R), rownames((sigtab_taxA_R)) %in% c('0238e0e03ffd3faa629954545d336e61', '052f174cab37be599600e58c78283fa2', '0592d48e1457e0fafb90b4158fa522c2',
                                                                                     '084e19e29dc9ebe03f4401003d10bff6', '26be318a519d7fe51cc6cf5d2378d1c8', '275c04dcabb809d184f0b6838763e20e',
                                                                                     '2e7210652ae3b77f31c42743c148406b', '321da2e457e5a9b769814b25476c4b10', '33d9e0d38932e8ce14c359de566b7d08',
                                                                                     '34c559c02664a1ac5ece941ca9000309', '4b8d75e30b64a18cf561c75cdc17043f', '4bb91812872f443514a3977be8c58774',
                                                                                     '4ce53584fbaa2aa3650f10bbe615c714', '57e66fdc78f87cd026455c6394730932', '5fca9caffa56a57fcc31d7ccba92d008',
                                                                                     '770af6feae23fae0ab3f1282a0ccbf18', '79eb38e43351aa3b12eae197935b81fb', '893c52ddb9dd678876c58f35b7ecbad6',
                                                                                     '89514854a16cbb3269c2e9e94a05e9d7', '908664fbed1b4350a0be7c1dc38094a8', '9690acde73fcd49a81835468c4b92397',
                                                                                     '9f9b3c564446dde56bbc0ef69261369f', 'a79e5838a5e0b50040e7f9aecf923028', 'b7f611ae7c6166d62354f04ca25391d8',
                                                                                     'ba37b62f122aca2aaff8ad84244df273', 'bd5b2a1bc31a73a4f37561a50e8e238c', 'c0664e6873e08cb34f7333015aabb07c',
                                                                                     'c36dba3bdc497773e638b8c441a9a0eb', 'c752096e70b99e9b3feeb36fb2beb2e1', 'd8a16afe0d36d2502e504377df9e7e44'))
specialist_res_subsetR_C <- subset((sigtab_taxR_C), rownames((sigtab_taxR_C)) %in% c('0238e0e03ffd3faa629954545d336e61', '052f174cab37be599600e58c78283fa2', '0592d48e1457e0fafb90b4158fa522c2',
                                                                                     '084e19e29dc9ebe03f4401003d10bff6', '26be318a519d7fe51cc6cf5d2378d1c8', '275c04dcabb809d184f0b6838763e20e',
                                                                                     '2e7210652ae3b77f31c42743c148406b', '321da2e457e5a9b769814b25476c4b10', '33d9e0d38932e8ce14c359de566b7d08',
                                                                                     '34c559c02664a1ac5ece941ca9000309', '4b8d75e30b64a18cf561c75cdc17043f', '4bb91812872f443514a3977be8c58774',
                                                                                     '4ce53584fbaa2aa3650f10bbe615c714', '57e66fdc78f87cd026455c6394730932', '5fca9caffa56a57fcc31d7ccba92d008',
                                                                                     '770af6feae23fae0ab3f1282a0ccbf18', '79eb38e43351aa3b12eae197935b81fb', '893c52ddb9dd678876c58f35b7ecbad6',
                                                                                     '89514854a16cbb3269c2e9e94a05e9d7', '908664fbed1b4350a0be7c1dc38094a8', '9690acde73fcd49a81835468c4b92397',
                                                                                     '9f9b3c564446dde56bbc0ef69261369f', 'a79e5838a5e0b50040e7f9aecf923028', 'b7f611ae7c6166d62354f04ca25391d8',
                                                                                     'ba37b62f122aca2aaff8ad84244df273', 'bd5b2a1bc31a73a4f37561a50e8e238c', 'c0664e6873e08cb34f7333015aabb07c',
                                                                                     'c36dba3bdc497773e638b8c441a9a0eb', 'c752096e70b99e9b3feeb36fb2beb2e1', 'd8a16afe0d36d2502e504377df9e7e44'))
specialist_res_subsetA_C <- subset((sigtab_taxA_C), rownames((sigtab_taxA_C)) %in% c('0238e0e03ffd3faa629954545d336e61', '052f174cab37be599600e58c78283fa2', '0592d48e1457e0fafb90b4158fa522c2',
                                                                                     '084e19e29dc9ebe03f4401003d10bff6', '26be318a519d7fe51cc6cf5d2378d1c8', '275c04dcabb809d184f0b6838763e20e',
                                                                                     '2e7210652ae3b77f31c42743c148406b', '321da2e457e5a9b769814b25476c4b10', '33d9e0d38932e8ce14c359de566b7d08',
                                                                                     '34c559c02664a1ac5ece941ca9000309', '4b8d75e30b64a18cf561c75cdc17043f', '4bb91812872f443514a3977be8c58774',
                                                                                     '4ce53584fbaa2aa3650f10bbe615c714', '57e66fdc78f87cd026455c6394730932', '5fca9caffa56a57fcc31d7ccba92d008',
                                                                                     '770af6feae23fae0ab3f1282a0ccbf18', '79eb38e43351aa3b12eae197935b81fb', '893c52ddb9dd678876c58f35b7ecbad6',
                                                                                     '89514854a16cbb3269c2e9e94a05e9d7', '908664fbed1b4350a0be7c1dc38094a8', '9690acde73fcd49a81835468c4b92397',
                                                                                     '9f9b3c564446dde56bbc0ef69261369f', 'a79e5838a5e0b50040e7f9aecf923028', 'b7f611ae7c6166d62354f04ca25391d8',
                                                                                     'ba37b62f122aca2aaff8ad84244df273', 'bd5b2a1bc31a73a4f37561a50e8e238c', 'c0664e6873e08cb34f7333015aabb07c',
                                                                                     'c36dba3bdc497773e638b8c441a9a0eb', 'c752096e70b99e9b3feeb36fb2beb2e1', 'd8a16afe0d36d2502e504377df9e7e44'))
#Export .csv of specialists only (fingers crossed)
write.csv(as.data.frame(specialist_res_subsetR_A), 
          file = "DESeq2_specialists_R_A.csv")
write.csv(as.data.frame(specialist_res_subsetA_R), 
          file = "DESeq2_specialists_A_R.csv")
write.csv(as.data.frame(specialist_res_subsetR_C), 
          file = "DESeq2_specialists_R_C.csv")
write.csv(as.data.frame(specialist_res_subsetA_C), 
          file = "DESeq2_specialists_A_C.csv")
#okay this all worked, is good. One point that has become apparent from this is that I need to double check which direction the logfold change is based on the contrasts chosen..
#I *think* they might actually be the revserve of what previously wrote
#

####Leftover DESEQ testing and removed commands

#for printing out results, alpha is a filter, dont run this if you want to get everything and sort them out yourself
alpha = 0.06
sigtab = res[which(res$padj < alpha), ]

# WIP FOR SUBSETTING: Subset deseq results BY ASV to bypass need for manual creation of specialist OTU and TAX tables..
# Should be possible using the base R subset command (as opposed to subset_taxa) and specifiying your list of individual ASVs
# Functions fine for results themselves, but still working on fixing agglomerated taxa for the subset...may now be fixed!
#BPrepare W14 subsets
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")

#Bind taxonomy to results
res = cbind(as(res, "data.frame"), as(tax_table(OBJ_W14)[rownames(res), ], "matrix"))
res

#Aglommerate to certain taxa intead of ASV #I THINK THIS IS THE STUMBLING BLOCK, is overwriting deseq with the entire physeq objsect
#when we just want it to append taxa at the family level to match
res <- tax_glom(OBJ_W14,taxrank = "Family")
res
write.csv(as.data.frame(res), 
          file = "DESeq2_specialist_subset_TEST.csv")

#WIP TEST equivalence vs the contrast command, name
#their example was these two commands being essentially equivalent except for some stuff about how contrast sets LFC to zero...
#not sure what this means in practice, hence why I want to test what differs
res <- results(dds, name = "condition_treated_vs_untreated")
res <- results(dds, contrast = c("condition","treated","untreated"))
#This should print out valid comparison names
resultsNames(diagdds)

sigtabR_A <- results(diagdds, name = "Location_Anode_vs_Root")
sigtabR_C <- results(diagdds, name = "Location_Cathode_vs_Root")

#These ones work only when using Anode as the baseline
sigtabA_R <- results(diagdds, name = "Location_Root_vs_Anode")
sigtabA_C <- results(diagdds, name = "Location_Cathode_vs_Anode")


#TO DO (fix) Export .csv (now that I've figured otu speicialist subset need to reformat this to correctly point to the full datasheet...TO DO)
write.csv(as.data.frame(sigtab_otuA), 
          file = "DESeq2_resultsA.csv")
write.csv(as.data.frame(sigtab_otuB), 
          file = "DESeq2_resultsB.csv")
write.csv(as.data.frame(sigtab_otuC), 
          file = "DESeq2_resultsC.csv")

#convert to matrix
sigtab_otu = cbind(as(sigtab, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtab), ], "matrix"))

write.csv(as.data.frame(res), 
          file = "raw_results.csv")

#only print out results with adjusted P higher than:
resSig <- subset(res, padj < 0.1)
resSig
write.csv(as.data.frame(resSig), 
          file = "sig_results.csv")

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
sigtab_otu$Phylum = factor(as.character(sigtab_otu$Phylum), levels = names(x))
# Class order
x = tapply(sigtab_otu$log2FoldChange, sigtab_otu$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_otu$Phylum = factor(as.character(sigtab_otu$Phylum), levels = names(x))
 
#Visulise
ggplot(sigtab_otu, aes(x = Phylum, y = log2FoldChange, color = Phylum)) + geom_point(size = 6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))


### DESeq for connection alone ----------------------------------------------
library(ggplot2)
library("DESeq2")
#Start here if haven't already (use raw vaules, but still subset out bad samples)
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")

#Agglomerate at desired taxa level (if you want to, otherwise proceed to and will get individual ASV changes)
OBJ1_exp_Genus <- tax_glom(OBJ1_exp,taxrank = "Genus")
OBJ1_exp_Family <- tax_glom(OBJ1_exp,taxrank = "Family")
OBJ1_exp_Order <- tax_glom(OBJ1_exp,taxrank = "Order")
OBJ1_exp_Class <- tax_glom(OBJ1_exp,taxrank = "Class")
#And create subset for use later when assigning taxa (NOT for analysis, just to bind taxa to results)
OBJ_Genus_W14 <- subset_samples(OBJ1_exp_Genus, Week == "Fourteen")

##Import to deseq and order factors to be compared (overall, not agglomerated)
diagdds = phyloseq_to_deseq2(OBJ1_exp, ~ Location + Connection) #Re: order of factors here, this would be testing for the effect of location, controlling for connection
#This time control for location differences, looking for connection response.

##Import agglomerated genus
diagdds = phyloseq_to_deseq2(OBJ1_exp_Genus, ~ Location + Connection) #Re: order of factors here, this would be testing for the effect of location, controlling for connection
#import agglomerated family
diagdds = phyloseq_to_deseq2(OBJ1_exp_Family , ~ Location + Connection) #Re: order of factors here, this would be testing for the effect of location, controlling for connection

diagdds$Location <- relevel(diagdds$Connection, ref = "Unconnected") # sets the reference point, baseline or control to be compared against

#Subset Week 14
diagdds <- diagdds[ , diagdds$Week == "Fourteen" ]
#check that subset was done as expected
as.data.frame( colData(diagdds) )

#Because of the way the data is nested between groups need to do some finangling to allow the model to distinguish between them correctly
#Have added a new column that groups samples by the soil column they were in (without location data, so e.g W14UP4 for all three cathode/anode/root sampeles)

#Run model and factors
diagdds = DESeq(diagdds, test = "Wald", fitType = "parametric")
res = results(diagdds, cooksCutoff = FALSE)
res # print out results

#Shouldnt need these anymore for agglomerated, but retain for overall
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")

#Bind taxonomy to results
res = cbind(as(res, "data.frame"), as(tax_table(OBJ_Genus_W14)[rownames(res), ], "matrix"))
res

#Export .csv
write.csv(as.data.frame(res), 
          file = "DESeq2_Genus_14.csv")

#For Agglomerated
#Different Comparison Direction Sheets
sigtabC_U = results(diagdds, contrast = c("Connection","Connected","Unconnected")) #a positive number here should represent an Increase in Connected FROM Unconnected
sigtabU_C = results(diagdds, contrast = c("Connection","Unconnected","Connected")) #postive here should be switched, so a postive = Increase in Unconnected FROM Connected
# So i think A should be the one to use...feels more logical to look at change towards connection (and makes discussion of electroactive enrichemnet easier)

#Bind results sheets with OTU Taxa
sigtab_taxC_U = cbind(as(sigtabC_U, "data.frame"), as(tax_table(OBJ_Genus_W14)[rownames(sigtabC_U), ], "matrix"))
sigtab_taxC_U

sigtab_taxU_C = cbind(as(sigtabU_C, "data.frame"), as(tax_table(OBJ_Genus_W14)[rownames(sigtabU_C), ], "matrix"))
sigtab_taxU_C

#Export .csv
write.csv(as.data.frame(sigtab_taxC_U), 
          file = "DESeq2_resultsC_U.csv")
write.csv(as.data.frame(sigtab_taxU_C), 
          file = "DESeq2_resultsU_C.csv")


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
diagdds$Cd_dose <- relevel(diagdds$Cd_dose, ref = "0_mg") ##reset reference level
diagdds = DESeq(diagdds, test = "Wald", fitType = "parametric")
resultsNames(diagdds)

sigtab = results(diagdds, contrast = c("Cd_dose","100_mg","0_mg"))
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

melt <- psmelt(OBJ1_DE2)
head(melt)
melt <- melt[sort.list(melt[,11]), ]

##ordering the x axis:

melt$Fig1 <- as.character(melt$Fig1)
#Then turn it back into an ordered factor
melt$Fig1 <- factor(melt$Fig1, levels = unique(melt$Fig1))
melt$Fig1 <- factor(melt$Fig1, levels = c("0_mg1","0_mg2","0_mg3","0_mg4","0_mg5","20_mg1","20_mg2","20_mg3","20_mg4",
                                          "20_mg5","100_mg1","100_mg2","100_mg3","100_mg4","100_mg5"))

t1 <- ggplot(melt, aes(Fig1, Abundance/5, fill = Order)) + geom_bar(stat = "identity",)
t1 <- t1 + facet_grid(Location~.) + theme(axis.text.x = element_text(angle = 45, size = 9,vjust = 0.7)) + scale_fill_manual(values = colvec)
t1 <- t1 + labs(title = NULL, x = expression(Cadmium~dose~(mg~kg^{-1}~soil)), y = "Relative abundance")
t1 <- t1 + scale_x_discrete(labels = c(rep("0 mg",5),rep("20 mg",5), rep("100 mg",5)))
t1

# ANCOM-BC ----------------------------------------------------------------
#As per: https://github.com/FrederickHuangLin/ANCOMBC/issues/19
#Using transofrmed/fractional counts is NOT recommended for ANCOM, it is 
#The ANCOM-BC methodology is developed for differential abundance analysis with regards to absolute abundances,
#i.e. reading in raw counts data is required. Using fractional relative abundances is highly not recommended.
#Re: formula structure when you are putting in your data, see: https://github.com/FrederickHuangLin/ANCOMBC/discussions/24
#yuor FIRST variable in the fomrula is what you are looking for a difference in, while adjusting for the effects of further variables
#Therefore, in your case, if the variable of interest is location, and you specify formula = "location + connection", it means you are trying to detect differentially abundant taxa with regards to location while adjusting connection effect

#Note on "setting" your reference point if you wish to do so.....ANCOMBC always does it's compairons in alphabetical order
#So if for example you want logfold change AWAY FROM uninoculated to Geo/Pseudo...it will not do this by default, and will instead use Geobacter inoculum as its;'s baseline/reference
#to fix this simply rename the uninoculated variable to A_uninoculated so that it is compareed correctly

#IMPORTANT NOTE: because of the above note about alphabetical ordering of groups will no be using a customised treatment table file
#3from now on, when the matter of concern is to change the ordering of ANCOM results in any way you must modifiy the column names of the .csv
#The relevant .csv file for re-ordering groups alphabetically is called: ANCOM_mapping_file.csv
#As of now, the renamed columns will be:
#Location column: A_Root , B_Anode (this means that the plant root is used as the the baseline for location (this will allow for all comparison directions)
#Inoculum column: A_Uninoculated (this will mean that the uninolculated treatment will now correctly be used as the baseline for inoculum)

setwd("~/Documents/University/Analysis/PMFC_18/2020 rerun outputs/Format for phyloseq")
library(phyloseq)
library(ape)
otu_table <- as.data.frame(read.csv("raw_readmap.csv", header = TRUE,row.names = "OTU_ID"))
taxmat <- as.matrix(read.csv("tax_table.csv", row.names = 1, header = TRUE))
treat <- as.data.frame(read.csv("ANCOM_mapping_file.csv", row.names = 1, header = TRUE))
TREE <- read.tree("rooted_tree.nwk")
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
TREAT = sample_data(treat)
OBJ1 = phyloseq(OTU,TAX,TREAT,TREE)
library(magrittr)
OBJ1 <- OBJ1 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
    Family  != "mitochondria" &
    Class   != "Chloroplast"
  )
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
#Half cut that retains Pseudomonas (use these ones moving forward from now!)
OBJ_Overall_TRIM <- subset_samples(OBJ1_exp, Treatment_Half_Trim == "Retain")
OBJ_W0 <- subset_samples(OBJ1_exp, Week == "Zero")
OBJ_W0_TRIM <- subset_samples(OBJ_W0, Treatment_Half_Trim == "Retain")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")
OBJ_W14_TRIM <- subset_samples(OBJ_W14, Treatment_Half_Trim == "Retain")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ANCOMBC")

library(ANCOMBC)
library(microbiome)

#### Final Cut Set ANCOM-BC LOCATION --------------------------------------------------
#FINAL half-cut DATASET
#HALFCUT DATASET, looking at location, controlling for connection/inoculum
out4 = ancombc(phyloseq = OBJ_W14_TRIM, formula = "Location + Connection + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 10000,group = "Location", 
              struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
res_global_tax4

#binding the tax table from phyloseq is suddenly no longer working so maybe will have to export seperately?
Tax_trim_forstitch = as(tax_table(OBJ_W14_TRIM), "matrix")
# transpose if necessary
#only if transposing:# if(taxa_are_rows(OBJ1_exp_tss)){OTU1 <- t(OTU1)}
# Coerce to data.frame
Tax_trim_forstitch = as.data.frame(Tax_trim_forstitch)
write.csv(Tax_trim_forstitch,"Tax_ANCOM_stitch.csv", row.names = TRUE)



#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#Export as single excel file in separate worksheets
library(openxlsx)

res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_coef4), ], "matrix"))
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_q4), ], "matrix"))
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_se4), ], "matrix"))
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_w4), ], "matrix"))
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_diff4), ], "matrix"))
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_p4), ], "matrix"))

#Name blank ASV dataframe column as "ASV" so that it exports along with rest of dataframe
colnames(res_global_tax4)
res_global_tax4 <- data.frame("ASV" = rownames(res_global_tax4), res_global_tax4)
tab_coef4 <- data.frame("ASV" = rownames(tab_coef4), tab_coef4)
tab_q4 <- data.frame("ASV" = rownames(tab_q4), tab_q4)
tab_se4 <- data.frame("ASV" = rownames(tab_se4), tab_se4)
tab_w4 <- data.frame("ASV" = rownames(tab_w4), tab_w4)
tab_diff4 <- data.frame("ASV" = rownames(tab_diff4), tab_diff4)
tab_p4 <- data.frame("ASV" = rownames(tab_p4), tab_p4)
colnames(res_global_tax4)

dataset_names <- list('Global' = res_global_tax4, 'Coefficients' = tab_coef4,'P_value' = tab_p4, 'Q_AdjustedP' = tab_q4,'SE_Residuals' = tab_se4,'W_stat' = tab_w4,'Differential_Abundance' = tab_diff4)
write.xlsx(dataset_names, file = 'ANCOMBC_W14_Location_FINALcut_ASV.xlsx', asTable = TRUE, firstRow = TRUE)

#agglomerate half cut dataset to Family for higher level look at differential abundance
#glom first
OBJ_W14_TRIM_fam <- tax_glom(OBJ_W14_TRIM,taxrank = "Family")
#ensure referring to agglomerated object
out4 = ancombc(phyloseq = OBJ_W14_TRIM_fam, formula = "Location + Connection + Inoculum", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Location", 
              struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#Export as single excel file in separate worksheets
library(openxlsx)

res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_coef4), ], "matrix"))
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_q4), ], "matrix"))
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_se4), ], "matrix"))
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_w4), ], "matrix"))
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_diff4), ], "matrix"))
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_p4), ], "matrix"))

colnames(res_global_tax4)
res_global_tax4 <- data.frame("ASV" = rownames(res_global_tax4), res_global_tax4)
tab_coef4 <- data.frame("ASV" = rownames(tab_coef4), tab_coef4)
tab_q4 <- data.frame("ASV" = rownames(tab_q4), tab_q4)
tab_se4 <- data.frame("ASV" = rownames(tab_se4), tab_se4)
tab_w4 <- data.frame("ASV" = rownames(tab_w4), tab_w4)
tab_diff4 <- data.frame("ASV" = rownames(tab_diff4), tab_diff4)
tab_p4 <- data.frame("ASV" = rownames(tab_p4), tab_p4)
colnames(res_global_tax4)

dataset_names <- list('Global' = res_global_tax4, 'Coefficients' = tab_coef4,'P_value' = tab_p4, 'Q_AdjustedP' = tab_q4,'SE_Residuals' = tab_se4,'W_stat' = tab_w4,'Differential_Abundance' = tab_diff4)
write.xlsx(dataset_names, file = 'ANCOMBC_W14_Location_FINALcut_Family.xlsx', asTable = TRUE, firstRow = TRUE)

###### WIP Connection ancom family level -------------------------------------------
OBJ_W14_TRIM_fam <- tax_glom(OBJ_W14_TRIM,taxrank = "Family")
#ensure referring to agglomerated object
out4 = ancombc(phyloseq = OBJ_W14_TRIM_fam, formula = "Connection + Location + Inoculum", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Connection", 
              struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global
out4

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#Export as single excel file in separate worksheets
library(openxlsx)

res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_coef4), ], "matrix"))
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_q4), ], "matrix"))
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_se4), ], "matrix"))
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_w4), ], "matrix"))
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_diff4), ], "matrix"))
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_p4), ], "matrix"))

colnames(res_global_tax4)
res_global_tax4 <- data.frame("ASV" = rownames(res_global_tax4), res_global_tax4)
tab_coef4 <- data.frame("ASV" = rownames(tab_coef4), tab_coef4)
tab_q4 <- data.frame("ASV" = rownames(tab_q4), tab_q4)
tab_se4 <- data.frame("ASV" = rownames(tab_se4), tab_se4)
tab_w4 <- data.frame("ASV" = rownames(tab_w4), tab_w4)
tab_diff4 <- data.frame("ASV" = rownames(tab_diff4), tab_diff4)
tab_p4 <- data.frame("ASV" = rownames(tab_p4), tab_p4)
colnames(res_global_tax4)

dataset_names <- list('Coefficients' = tab_coef4,'P_value' = tab_p4, 'Q_AdjustedP' = tab_q4,'SE_Residuals' = tab_se4,'W_stat' = tab_w4,'Differential_Abundance' = tab_diff4)
write.xlsx(dataset_names, file = 'ANCOMBC_W14_Connection_FINALcut_Family.xlsx', asTable = TRUE, firstRow = TRUE)

##### WIP inocula ANCOM family level -----------------------------------------------------------------
OBJ_W14_TRIM_fam <- tax_glom(OBJ_W14_TRIM,taxrank = "Family")
#ensure referring to agglomerated object
out4 = ancombc(phyloseq = OBJ_W14_TRIM_fam, formula = "Inoculum + Location + Connection", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Inoculum", 
              struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global
out4

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#Export as single excel file in separate worksheets
library(openxlsx)

res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_coef4), ], "matrix"))
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_q4), ], "matrix"))
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_se4), ], "matrix"))
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_w4), ], "matrix"))
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_diff4), ], "matrix"))
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_p4), ], "matrix"))

colnames(res_global_tax4)
res_global_tax4 <- data.frame("ASV" = rownames(res_global_tax4), res_global_tax4)
tab_coef4 <- data.frame("ASV" = rownames(tab_coef4), tab_coef4)
tab_q4 <- data.frame("ASV" = rownames(tab_q4), tab_q4)
tab_se4 <- data.frame("ASV" = rownames(tab_se4), tab_se4)
tab_w4 <- data.frame("ASV" = rownames(tab_w4), tab_w4)
tab_diff4 <- data.frame("ASV" = rownames(tab_diff4), tab_diff4)
tab_p4 <- data.frame("ASV" = rownames(tab_p4), tab_p4)
colnames(res_global_tax4)

dataset_names <- list('Global' = res_global_tax4, 'Coefficients' = tab_coef4,'P_value' = tab_p4, 'Q_AdjustedP' = tab_q4,'SE_Residuals' = tab_se4,'W_stat' = tab_w4,'Differential_Abundance' = tab_diff4)
write.xlsx(dataset_names, file = 'ANCOMBC_W14_Inoculum_FINALcut_Family.xlsx', asTable = TRUE, firstRow = TRUE)

#### WIP inocula ANCOM GENUS level -----------------------------------------------------------------
OBJ_W14_TRIM_Genus <- tax_glom(OBJ_W14_TRIM,taxrank = "Genus")
#ensure referring to agglomerated object
out4 = ancombc(phyloseq = OBJ_W14_TRIM_Genus, formula = "Inoculum + Location + Connection", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Inoculum", 
              struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global
out4

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_Genus)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#Export as single excel file in separate worksheets
library(openxlsx)

res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_Genus)[rownames(res_global4), ], "matrix"))
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM_Genus)[rownames(tab_coef4), ], "matrix"))
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM_Genus)[rownames(tab_q4), ], "matrix"))
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM_Genus)[rownames(tab_se4), ], "matrix"))
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM_Genus)[rownames(tab_w4), ], "matrix"))
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM_Genus)[rownames(tab_diff4), ], "matrix"))
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM_Genus)[rownames(tab_p4), ], "matrix"))

colnames(res_global_tax4)
res_global_tax4 <- data.frame("ASV" = rownames(res_global_tax4), res_global_tax4)
tab_coef4 <- data.frame("ASV" = rownames(tab_coef4), tab_coef4)
tab_q4 <- data.frame("ASV" = rownames(tab_q4), tab_q4)
tab_se4 <- data.frame("ASV" = rownames(tab_se4), tab_se4)
tab_w4 <- data.frame("ASV" = rownames(tab_w4), tab_w4)
tab_diff4 <- data.frame("ASV" = rownames(tab_diff4), tab_diff4)
tab_p4 <- data.frame("ASV" = rownames(tab_p4), tab_p4)
colnames(res_global_tax4)

dataset_names <- list('Global' = res_global_tax4, 'Coefficients' = tab_coef4,'P_value' = tab_p4, 'Q_AdjustedP' = tab_q4,'SE_Residuals' = tab_se4,'W_stat' = tab_w4,'Differential_Abundance' = tab_diff4)
write.xlsx(dataset_names, file = 'ANCOMBC_W14_Inoculum_FINALcut_Genus.xlsx', asTable = TRUE, firstRow = TRUE)

# Heirarchical Clustering -------------------------------------------------

##HEIRACHICAL CLUSTERING____________________________________________________________________
library(vegan)

##make disimialrity matrix
OBJ1_r_b <- distance(OBJ1_r, "bray")

# This is the actual hierarchical clustering call, specifying average-link clustering

OBJ1.hclust <- hclust(OBJ1_r_b, method = "average")
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
colorscale <- c("green4","red3", "darkorange1","green4", "red3","darkorange1")
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
par(mar = c(1,1,1,7))
dend %>%
  set("labels_col", "black")  %>% #label colours
  set("labels_cex", 0.6)  %>% #label text size
  set("branches_k_color", value = rep(c("#138d75","goldenrod1", "darkorange" ),2), k = 6) %>% #colur k major splits
  set("branches_lwd", 1.8) %>% #branch line weight
  set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 1) %>%  # node point size
  set("leaves_col", labels_colors(dend)) %>% # node point color
  plot(horiz = TRUE, axes = TRUE)
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

#TEMP
OBJ1_Gen <- transform(OBJ1_exp_Genus, "compositional")#transform to relative abundance -microbiom pkg
OBJ1_Gen <- transform_sample_counts(OBJ1_r, function(x) x / sum(x) )#alternative way transform to relative abundance - phyloseq pkg
OBJ1_Gen 
OBJ1_Gen <- filter_taxa(OBJ1_Gen, function(x) mean(x) > 0.001, TRUE)
OBJ1_Fam <- tax_glom(OBJ1_Gen,taxrank = "Family")#concatenate OTUs at a higher phylogenetic level
OBJ1_Ord <- tax_glom(OBJ1_Gen,taxrank = "Order")#concatenate OTUs at a higher phylogenetic level

p1 <- plot_bar(OBJ1_Ord, "Location", fill = "Order")
p1 <- plot_bar(OBJ1_Ord, "Connection", fill = "Order", facet_grid = ~Location)
p1

p1 <- plot_bar(OBJ1_Fam, "Location", fill = "Family")
p1 <- plot_bar(OBJ1_Fam, "Connection", fill = "Family", facet_grid = ~Location)
p1

p1 <- plot_bar(OBJ1_Gen, "Location", fill = "Genus")
p1 <- plot_bar(OBJ1_Gen, "Connection", fill = "Genus", facet_grid = ~Location)
p1

OBJ1_ord <- tax_glom(OBJ1_rr,taxrank = "Order")#concatenate OTUs at a higher phylogenetic level
OBJ1_ord #check how make OTUs you have - too amny makes an untidy fig legend
#you can either concatenate at a higher tax rank or remove really rare otus thus:
OBJ1_ord2 <- filter_taxa(OBJ1_ord, function(x) mean(x) > 0.0005, TRUE)
OBJ1_ord2 #much more manageable number of OTUs

##PLOT1 - all samples
p1 <- plot_bar(OBJ1_ord2, "Sample", fill = "Class")
p1 <- p1 + scale_fill_manual(values = colvec)
p1 <- p1 + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1

##sometimes changes are subtle and hard to see
##one optinon is to plot only OTUs that are differetnially abundant (a lesson for later)

##PLOT2 - all samples faceted by phylum
OBJ1_ord2 <- filter_taxa(OBJ1_ord, function(x) mean(x) > 0.005, TRUE)
p2 <- plot_bar(OBJ1_ord2, "Sample", fill = "Class", facet_grid = ~Phylum)
p2 <- p2 + scale_fill_manual(values = colvec)
p2 <- p2 + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2

lim <- c("Bulk_0","Bulk_20","Bulk_100","Rhizo_0","Rhizo_20","Rhizo_100")#order the axis how you want them
lim_Cd_dose <- c("0_mg","20_mg","100_mg")#order the axis how you want them

##PLOT3 - Cd-dose groups faceted by Location
OBJ1_ord2 <- filter_taxa(OBJ1_ord, function(x) mean(x) > 0.005, TRUE)
p3 <- plot_bar(OBJ1_ord2, "Cd_dose", fill = "Order", facet_grid = ~Location)
p3 <- p3 + scale_x_discrete(name = "Title", limits = lim_Cd_dose)
p3 <- p3 + scale_fill_manual(values = colvec)
p3 <- p3 + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
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

melt <- psmelt(OBJ1_ord2)#extract the information from the phyloseq object
head(melt) #check column headings
melt <- melt[sort.list(melt[,10]), ] #reorder 'melt' by the column we will be plotting by (in this cae 'order si in col 10)

##ordering the x axis labels:
melt$Cd_dose <- as.character(melt$Cd_dose)
#Then turn it back into an ordered factor
melt$Cd_dose <- factor(melt$Cd_dose, levels = unique(melt$Cd_dose))
melt$Cd_dose <- factor(melt$Cd_dose, levels = c("0_mg","20_mg","100_mg"))

t1 <- ggplot(melt, aes(Cd_dose, Abundance/5, fill = Order)) + geom_bar(stat = "identity",)
t1 <- t1 + facet_grid(Location~.) + theme(axis.text.x = element_text(angle = 45, size = 9,vjust = 0.7)) + scale_fill_manual(values = colvec)
t1 <- t1 + labs(title = NULL, x = expression(Cadmium~dose~(mg~kg^{-1}~soil)), y = "Relative abundance")
t1 <- t1 + scale_x_discrete(labels = c(rep("0 mg",5),rep("20 mg",5), rep("100 mg",5)))
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
diagdds$Cd_dose <- relevel(diagdds$Cd_dose, ref = "0_mg") ##reset reference level
diagdds = DESeq(diagdds, test = "Wald", fitType = "parametric")
resultsNames(diagdds)

sigtab = results(diagdds, contrast = c("Cd_dose","100_mg","0_mg"))
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

melt <- psmelt(OBJ1_DE2)
head(melt)
melt <- melt[sort.list(melt[,11]), ]

##ordering the x axis:

melt$Fig1 <- as.character(melt$Fig1)
#Then turn it back into an ordered factor
melt$Fig1 <- factor(melt$Fig1, levels = unique(melt$Fig1))
melt$Fig1 <- factor(melt$Fig1, levels = c("0_mg1","0_mg2","0_mg3","0_mg4","0_mg5","20_mg1",
                                          "20_mg2","20_mg3","20_mg4","20_mg5","100_mg1","100_mg2","100_mg3","100_mg4","100_mg5"))

t1 <- ggplot(melt, aes(Fig1, Abundance/5, fill = Order)) + geom_bar(stat = "identity",)
t1 <- t1 + facet_grid(Location~.) + theme(axis.text.x = element_text(angle = 45, size = 9,vjust = 0.7)) + scale_fill_manual(values = colvec)
t1 <- t1 + labs(title = NULL, x = expression(Cadmium~dose~(mg~kg^{-1}~soil)), y = "Relative abundance")
t1 <- t1 + scale_x_discrete(labels = c(rep("0 mg",5),rep("20 mg",5), rep("100 mg",5)))
t1

# PICRUSt2 ----------------------------------------------------------------

#Starting note on conversion and import before starting.
#Ouputs of PICRUSt2 pipeline with add.description. step should result in your tables for EC's, KO's, and Pathways.
#Each of these can be manually copied and pasted into pre-existing template .csv's (e.g any that you haave used with OTU's ASV's above)
#When doing this, main points to note are:
#Make sure to retain original readmap and tax sheet headers (e.g OTU_ID, and Taxa levels), as these are what phyloseq expects for import.
#Copy your function labels as OTUs for both sheets.
#Copy past function descriptions into EACH taxa level on your taxa sheet
#And finally, your abundances and sample labels for readmap
#Unless you've manually selected subsets, or done som eautomatic subsetting that differs from your above analysis you should be able to substitute the same metadata/mapping file as above

#### Import ------------------------------------------------------------------

setwd("~/Documents/University/Analysis/PMFC_18/2020 rerun outputs/Format for phyloseq")

# Readmap
EC_table <- as.data.frame(read.csv("readmap_EC.csv", header = TRUE,row.names = "OTU_ID"))
KO_table <- as.data.frame(read.csv("readmap_KO.csv", header = TRUE,row.names = "OTU_ID"))
Path_table <- as.data.frame(read.csv("readmap_Path.csv", header = TRUE,row.names = "OTU_ID"))

# "Tax" function table
EC_mat <- as.matrix(read.csv("tax_EC.csv", row.names = 1, header = TRUE))
KO_mat <- as.matrix(read.csv("tax_KO.csv", row.names = 1, header = TRUE))
Path_mat <- as.matrix(read.csv("tax_Path.csv", row.names = 1, header = TRUE))

# Metadata/Mapping/Treatment file
treat <- as.data.frame(read.csv("mapping_file.csv", row.names = 1, header = TRUE))

# Make Phylo object
library(phyloseq)

#EC Phylo object
EC_OTU = otu_table(EC_table, taxa_are_rows = TRUE)
EC_TAX = tax_table(EC_mat)
TREAT = sample_data(treat)
EC_PHYLO = phyloseq(EC_OTU,EC_TAX,TREAT)

#KO Phylo object
KO_OTU = otu_table(KO_table, taxa_are_rows = TRUE)
KO_TAX = tax_table(KO_mat)
TREAT = sample_data(treat)
KO_PHYLO = phyloseq(KO_OTU,KO_TAX,TREAT)

#Path Phylo object
Path_OTU = otu_table(Path_table, taxa_are_rows = TRUE)
Path_TAX = tax_table(Path_mat)
TREAT = sample_data(treat)
Path_PHYLO = phyloseq(Path_OTU,Path_TAX,TREAT)

#BEFORE PROCEEDING need to transform, first: TSS AND then LOG+1/anderson log on the phyloseq object 
#(this just here to test that the standardise function was doing the exact same thing)
Path_PHYLO_tss_manual = transform_sample_counts(Path_PHYLO, function(OTU) OTU/sum(OTU) )
Path_PHYLO_tss_manual
#ALTERNATIVELY can use metagMisc  package to do both TSS and anderson log
#https://github.com/vmikk/metagMisc/blob/master/man/phyloseq_standardize_otu_abundance.Rd
#https://github.com/vmikk/metagMisc
#https://rdrr.io/github/vmikk/metagMisc/man/phyloseq_standardize_otu_abundance.html
#anderson log
#Old version of the below commands
#Path_PHYLO_tss <- phyloseq_standardize_otu_abundance(Path_PHYLO, method = "total")
#Path_PHYLO_tss 
#Path_PHYLO_log <- phyloseq_standardize_otu_abundance(Path_PHYLO_tss, method = "log")
#Path_PHYLO_log
#devtools::install_github("vmikk/metagMisc")

install.packages("remotes")
remotes::install_github("vmikk/metagMisc")
library(metagMisc)

#and cut out unwanted samples
Path_PHYLO_log <- subset_samples(Path_PHYLO_log, Experiment == "Y")

#### Subset ------------------------------------------------------------------

#Combined treatments
PATH_Unin_Conn <- subset_samples(Path_PHYLO_log, Treatment == "Uninoculated Connected")
PATH_Unin_Unconn <- subset_samples(Path_PHYLO_log, Treatment == "Uninoculated Unconnected")
PATH_Geo_Conn <- subset_samples(Path_PHYLO_log, Treatment == "Geobacter Connected")
PATH_Geo_Unconn <- subset_samples(Path_PHYLO_log, Treatment == "Geobacter Unconnected")
PATH_Mont_Conn <- subset_samples(Path_PHYLO_log, Treatment == "Montebello Connected")
PATH_Mont_Unconn <- subset_samples(Path_PHYLO_log, Treatment == "Montebello Unconnected")
PATH_Pseudo_Conn <- subset_samples(Path_PHYLO_log, Treatment == "Pseudomonas Connected")
PATH_Pseudo_Unconn <- subset_samples(Path_PHYLO_log, Treatment == "Pseudomonas Unconnected")

#Locations
PATH_anode_tss <- subset_samples(Path_PHYLO_log, Location == "Anode")
PATH_cathode_tss <- subset_samples(Path_PHYLO_log, Location == "Cathode")

#Connection
PATH_connected_tss <- subset_samples(Path_PHYLO_log, Connection == "Connected")
PATH_unconnected_tss <- subset_samples(Path_PHYLO_log, Connection == "Unconnected")

#Week/Time
PATH_W0 <- subset_samples(Path_PHYLO_log, Week == "Zero")
PATH_W6 <- subset_samples(Path_PHYLO_log, Week == "Six")
PATH_W8 <- subset_samples(Path_PHYLO_log, Week == "Eight")
PATH_W14 <- subset_samples(Path_PHYLO_log, Week == "Fourteen")

#Wk14 Connected and Unconnected Subsets
PATH_W14_connected <- subset_samples(PATH_W14, Connection == "Connected")
PATH_W14_unconnected <- subset_samples(PATH_W14, Connection == "Unconnected")

#Treatment
PATH_Unin <- subset_samples(Path_PHYLO_log, Inoculum == "Uninoculated")
PATH_Geo <- subset_samples(Path_PHYLO_log, Inoculum == "Geobacter")
PATH_Mont <- subset_samples(Path_PHYLO_log, Inoculum == "Montebello")
PATH_Pseudo <- subset_samples(Path_PHYLO_log, Inoculum == "Pseudomonas")

### NMDS + Ordination Plots -------------------------------------------------------------

library("ggplot2")
library("RColorBrewer")

#subsetted timepoints WEIGHTED
PATH_NMDS_W14w <- ordinate(PATH_W14, "NMDS", distance = "bray", binary = FALSE)
PATH_NMDS_W14w

#Week 0
PATH_NMDS_W0w <- ordinate(PATH_W0, "NMDS", distance = "bray", binary = FALSE)
PATH_NMDS_W0w

#Subsetted timepoints UNWEIGHTED
PATH_NMDS_W14u <- ordinate(PATH_W14, "NMDS", distance = "bray", binary = TRUE)
PATH_NMDS_W14u

#Week 0
PATH_NMDS_W0u <- ordinate(PATH_W0, "NMDS", distance = "bray", binary = TRUE)
PATH_NMDS_W0u

#Plot Ordination weighted
path_W14w_ord <- plot_ordination(PATH_W14, PATH_NMDS_W14w, color = "Treatment", shape = "Location", label = NULL)
path_W14w_ord
path_W14w_ord + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47",
                                                                                                                              "#00B85C","#141414","#7A7A7A"))

#Plot Ordination unweighted
path_W14u_ord <- plot_ordination(PATH_W14, PATH_NMDS_W14u, color = "Treatment", shape = "Location", label = NULL)
path_W14u_ord
path_W14u_ord + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47",
                                                                                                                              "#00B85C","#141414","#7A7A7A"))

#Time Zero weighted
path_W0w_ord <- plot_ordination(PATH_W0, PATH_NMDS_W0w, color = "Treatment", shape = "Location", label = NULL)
path_W0w_ord
path_W0w_ord + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47",
                                                                                                                             "#00B85C","#141414","#7A7A7A"))

#Time Zero unweighted
path_W0u_ord <- plot_ordination(PATH_W0, PATH_NMDS_W0u, color = "Treatment", shape = "Location", label = NULL)
path_W0u_ord
path_W0u_ord + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47",
                                                                                                                             "#00B85C","#141414","#7A7A7A"))

#Connection alone at W14
#Unconnected
PATH_NMDS_W14_un_w <- ordinate(PATH_W14_unconnected, "NMDS", distance = "bray", binary = FALSE)
PATH_NMDS_W14_un_u <- ordinate(PATH_W14_unconnected, "NMDS", distance = "bray", binary = TRUE)

PATH_NMDS_W14_un_w
PATH_NMDS_W14_un_u

#Connected
PATH_NMDS_W14_con_w <- ordinate(PATH_W14_connected, "NMDS", distance = "bray", binary = FALSE)
PATH_NMDS_W14_con_u <- ordinate(PATH_W14_connected, "NMDS", distance = "bray", binary = TRUE)

PATH_NMDS_W14_con_w
PATH_NMDS_W14_con_u

## PERMANOVA ---------------------------------------------------------------

library(vegan)
#Setup objects for testing significance, effect size, correlation.
#First pull the variables you want to test into objects:
#By location at end and start
Location <- get_variable(PATH_W14, "Location")
Location0 <- get_variable(PATH_W0, "Location")

#By connection at end
Connection <- get_variable(PATH_W14, "Connection")
Connection0 <- get_variable(PATH_W0, "Connection")

#By Inoculum at end
Inoculum <- get_variable(PATH_W14, "Inoculum")
Inoculum0 <- get_variable(PATH_W0, "Inoculum")

#W14
PATH_W14perm_w <- distance(PATH_W14, "bray" , binary = FALSE)
PATH_W14perm_u <- distance(PATH_W14, "bray" , binary = TRUE)

PATH_W14_ado_w = adonis(PATH_W14perm_w ~ Location * Connection * Inoculum, permutations = 9999)
PATH_W14_ado_w

PATH_W14_ado_u = adonis(PATH_W14perm_u ~ Location * Connection * Inoculum, permutations = 9999)
PATH_W14_ado_u

#W0
PATH_W0perm_w <- distance(PATH_W0, "bray" , binary = FALSE)
PATH_W0perm_u <- distance(PATH_W0, "bray" , binary = TRUE)

PATH_W0_ado_w = adonis(PATH_W0perm_w ~ Location0 * Connection0 * Inoculum0, permutations = 9999)
PATH_W0_ado_w

PATH_W0_ado_u = adonis(PATH_W0perm_u ~ Location0 * Connection0 * Inoculum0, permutations = 9999)
PATH_W0_ado_u

# CLAM Clamtest for Specialist and Generalist categorisation - WIP --------
#https://rdrr.io/rforge/vegan/man/clamtest.html
#https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/clamtest
#commands
#documentation and example data:

library(vegan)

clamtest(comm, groups, coverage.limit = 10, specialization = 2/3,  npoints = 20, alpha = 0.05/20)
"summary"(object, ...)
"plot"(x, xlab, ylab, main,  pch = 21:24, col.points = 1:4,  col.lines = 2:4, lty = 1:3, position = "bottomright", ...)
#example

# Example using two mite groups. The mite data are available in vegan

data(mite)
data(mite.env)
sol <- with(mite.env, clamtest(mite, Shrub == "None", alpha = 0.005)) # with Shrub==none , Spec_TRUE species fall on the lower right,
#and Spec_FLASE species fall on upper left of plot, Shrub==NONE in this instance results in NONE specialists being called Spec_TRUE,
#all others lumped as Spec_FALS
summary(sol)
head(sol)
plot(sol, xlab = "None", ylab = "Other", position = NULL)

#Try to adapt to my data....
#first otu table must be transposed so that species are columns, samples = rows
#doesnt seem possible with excel...exceeds maximum column #
#transpose in R

t_otu_table <- as.data.frame(t(otu_table))

sol <- with(treat, clamtest(t_otu_table, Substrate == "Water", alpha = 0.005))

sol <- with(treat, clamtest(t_otu_table, Location == "Root", alpha = 0.005))

summary(sol)
head(sol)
plot(sol, position = "topright", bty = "n", xpd = TRUE)
plot(sol, xlab = "Anode", ylab = "Other", position = NULL, bty = "n", xpd = TRUE)
plot(sol, xlab = "Water", ylab = "Soil")

write.csv(sol,"File_.csv", row.names = TRUE)

# Other WIP / General Useful script section -------------------------------

## A Colour Palette Interlude ------------------------------------------------
## This section used to decide on what colour palettes ill apply to each condition/treaatment, and based on ordination sections primarily
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
#Uninoculated Connected: Black  = "#141414"
#Uninoculated Unconnected: Grey  = "#7A7A7A"
my_colvec <- c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A")

#### To use for fills, add
bp + scale_fill_manual(values = cbp2)
### To use for line and point colors, add
sp + scale_colour_manual(values = cbp2)
#Themes
# https://ggplot2.tidyverse.org/reference/ggtheme.html
#Colours
#### https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/ , https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html

## Image export resolution testing  -------------------------------------
##is based on species heeatmap plot section for specialists, but can replot and test with others
### WIP TEST for exporting control, DPI, etc,.... trying to find how to get 300 dpi exports to be correctly recognised as such by photoshop etc

library(Cairo)
# using the package Cairo
Cairo::Cairo(
  30, #length
  30, #width
  file = paste("Resolution_Test_Export", ".jpeg", sep = ""), #also can use .png
  type = "jpeg", #png is another option
  bg = "transparent", #white or transparent depending on your requirements 
  dpi = 300,
  units = "cm" #you can change this to pixels using px etc 
)
plot(speciesheatplot_ORD) #this is your graph object plotted beforehand, so should already appear in your data pane
dev.off()

###using tiff function
tiff("PLOTExample.tif", width = 4200, height = 2200, units = "px", res = 300)
plot(speciesheatplot_ORD)
dev.off()
#Neither if these options seem to genuinely export an image that reads as 300 dpi when opened in other apps though....

###One more test... The bitmap one seems to be the only one that genuinely registers as being 300dpi BUT it kinda looks like crap by compairson....
pdf(file = "FileName.pdf", width = 12, height = 17, family = "Helvetica") # defaults to 7 x 7 inches
plot(speciesheatplot_ORD)
dev.off()
 
bitmap("FileName.tiff", height = 2200, width = 4200, units = 'px', type = "tiff24nc", res = 300)
plot(speciesheatplot_ORD)
dev.off()
 
tiff("FileName.tiff", height = 2200, width = 4200, units = 'px', 
     compression = "lzw", res = 300)
plot(speciesheatplot_ORD)
dev.off()
####

# New notes and/or functions to fix up and integrate ----------------------

# Pairise perma, log transofrm, ancom-bc are main things to integrate as of now....

#Note to self: re: remember to add in pairwise permanova....speak to Josh if can't find and work out..

###################ANCOMBC################### DO NOT DO THIS YET (better for over time analysis)
(!requireNamespace("BiocManager", quietly = TRUE))
+     install.packages("BiocManager")
BiocManager::install("ANCOMBC")
 
library("ANCOMBC")
library(microbiome)
library(xlsx)
 
#####Example ANCOM-BC from Sarah....FAMILY/PICRUST
#Before ANCOM-BC do TSS and then Anderson log

#do TSS

#do log transform
OBJ2 <- phyloseq_standardize_otu_abundance(OBJ1, method = "log")

#Where OBJ1_X_PICRUST is your functional Phyloseq object
OBJ1_PICRUST = aggregate_taxa(OBJ1_X_PICRUST, "Family")
outPICRUST = ancombc(phyloseq = OBJ1_PICRUST, formula = "Treatment + CollectionPoint", p_adj_method = "BH", zero_cut = 0.90, lib_cut = 1000,group = "Treatment",
                     struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5,max_iter = 100, conserve = TRUE, alpha = 0.0504, global = TRUE)
resPICRUST = outPICRUST$res
resPICRUST_global = outPICRUST$res_global
outPICRUST2 = ancombc(phyloseq = OBJ1_PICRUST, formula = "CollectionPoint + Treatment", p_adj_method = "BH", zero_cut = 0.90, lib_cut = 1000, group = "CollectionPoint",
                      struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5,max_iter = 100, conserve = TRUE, alpha = 0.0504, global = TRUE)
resPICRUST2 = outPICRUST2$res
resPICRUST_global2 = outPICRUST2$res_global
write.xlsx(resPICRUST_global, "ANCOM_BC-PICRUST-BH.xlsx", sheetName = "Treatment",append = TRUE)
write.xlsx(resPICRUST_global2, "ANCOM_BC-PICRUST-BH.xlsx", sheetName = "CollectionPoint",append = TRUE)
 
#collating the pairwise results:
tab_coef = resPICRUST$beta
write.xlsx(tab_coef, "ANCOM_BC-PICRUST-BH.xlsx", sheetName = "coef",append = TRUE)
tab_se = resPICRUST$se
write.xlsx(tab_se, "ANCOM_BC-PICRUST-BH.xlsx", sheetName = "residuals",append = TRUE)
tab_w = resPICRUST$W
write.xlsx(tab_w, "ANCOM_BC-PICRUST-BH.xlsx", sheetName = "W_stat",append = TRUE)
tab_p = resPICRUST$p_val
write.xlsx(tab_p, "ANCOM_BC-PICRUST-BH.xlsx", sheetName = "p",append = TRUE)
tab_q = resPICRUST$q
write.xlsx(tab_q, "ANCOM_BC-PICRUST-BH.xlsx", sheetName = "adjusted_p",append = TRUE)
tab_diff = resPICRUST$diff_abn
write.xlsx(tab_diff, "ANCOM_BC-PICRUST-BH.xlsx", sheetName = "DA",append = TRUE)

# Logarithmic transformation as in Anderson et al., 2006
OBJ2 <- phyloseq_standardize_otu_abundance(OBJ1, method = "log")
#google this to check that it's log X+1 google
 
#PAIRWISE PERMANOVA
####Example documentation
install.packages(devtools)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
data(iris)
pairwise.adonis(iris[,1:4],iris$Species)

# For strata (blocks), following example of Jari Oksanen in adonis2. 
dat <- expand.grid(rep = gl(2,1), NO3 = factor(c(0,10,30)),field = gl(3,1) )
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3) + 2) + rnorm(18)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3) + 2) + rnorm(18)/2
Y <- data.frame(Agropyron, Schizachyrium)

pairwise.adonis2(Y ~ NO3, data = dat, strata = 'field')
#or more examples see also 
?pairwise.adonis() 
?pairwise.adonis2()
####Example documentation

#put data into vegan object, and put metadata into vegan object, $ID for column of metadata, eg Location
pairwise.adonis(OBJ1BacVegan, VeganSam$ID, perm = 9999, p.adjust.m = "BH")

#FOR ANCOM BC pre-treatment/transform
#TSS then Log+1
#Despite its name simply being "log", the Anderson "log" function from Anderson 2006 does seem to be a modified log transofrm
#it should accounts for zeros, it is not technically just a log+1 either accordingto documentation,
#but modified based on the number encountered and adjusted to repvent unusable transofmation (presumably)

# Broken or Deprecated script ---------------------------------------------
#Double check nothing useful within tnaythinghere before deleting

#DESEQ Subsetting for top specialists (NOT WORKING, MAY BE EASIER TO SUBSET THE FULL SET ABOVE AFTER BINDING TAX etc)
#SUBSETTING FOR SPECIALISTS ALONE
#This one will only ever work at ASV level, but may come in handy for any other work we do on the specilaist or "hyper" specilast subsets, lets us pull them out immediately
specialist_res_subset <- subset((res), rownames((res)) %in% c('0238e0e03ffd3faa629954545d336e61', '052f174cab37be599600e58c78283fa2', '0592d48e1457e0fafb90b4158fa522c2',
                                                              '084e19e29dc9ebe03f4401003d10bff6', '26be318a519d7fe51cc6cf5d2378d1c8', '275c04dcabb809d184f0b6838763e20e',
                                                              '2e7210652ae3b77f31c42743c148406b', '321da2e457e5a9b769814b25476c4b10', '33d9e0d38932e8ce14c359de566b7d08',
                                                              '34c559c02664a1ac5ece941ca9000309', '4b8d75e30b64a18cf561c75cdc17043f', '4bb91812872f443514a3977be8c58774',
                                                              '4ce53584fbaa2aa3650f10bbe615c714', '57e66fdc78f87cd026455c6394730932', '5fca9caffa56a57fcc31d7ccba92d008',
                                                              '770af6feae23fae0ab3f1282a0ccbf18', '79eb38e43351aa3b12eae197935b81fb', '893c52ddb9dd678876c58f35b7ecbad6',
                                                              '89514854a16cbb3269c2e9e94a05e9d7', '908664fbed1b4350a0be7c1dc38094a8', '9690acde73fcd49a81835468c4b92397',
                                                              '9f9b3c564446dde56bbc0ef69261369f', 'a79e5838a5e0b50040e7f9aecf923028', 'b7f611ae7c6166d62354f04ca25391d8',
                                                              'ba37b62f122aca2aaff8ad84244df273', 'bd5b2a1bc31a73a4f37561a50e8e238c', 'c0664e6873e08cb34f7333015aabb07c',
                                                              'c36dba3bdc497773e638b8c441a9a0eb', 'c752096e70b99e9b3feeb36fb2beb2e1', 'd8a16afe0d36d2502e504377df9e7e44'))
specialistsR_A = results(specialist_res_subset, contrast = c("Location","Root","Anode"))
specialistsR_C = results(specialist_res_subset, contrast = c("Location","Root","Cathode"))
specialistsA_C = results(specialist_res_subset, contrast = c("Location","Anode","Cathode"))

#Bind results sheets with OTU Taxa
specialistsR_A_tax = cbind(as(specialistsR_A, "data.frame"), as(tax_table(OBJ1_exp)[rownames(specialistsR_A), ], "matrix"))
specialistsR_A_tax

sigtab_otuB = cbind(as(sigtabB, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtabB), ], "matrix"))
sigtab_otuB

sigtab_otuC = cbind(as(sigtabC, "data.frame"), as(tax_table(OBJ1_exp)[rownames(sigtabC), ], "matrix"))
sigtab_otuC

#Export .csv
write.csv(as.data.frame(specialistsR_A_tax), 
          file = "DESeq2_results_R_A_spec.csv")
write.csv(as.data.frame(sigtab_otuB), 
          file = "DESeq2_resultsB.csv")
write.csv(as.data.frame(sigtab_otuC), 
          file = "DESeq2_resultsC.csv")

write.csv(as.data.frame(specialist_res_subset), 
          file = "DESeq2_specialist_subset_TEST.csv")
## END SUBSETTING WIP


###### Deprecated: ANCOM of Full Dataset OLD -----------------------------------------------


#FULL DATASET
#Unsure of formula structure re: direction change, etc...test if you can just test for location alone
out = ancombc(phyloseq = OBJ_W14, formula = "Location", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Location",
              struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5,max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res1 = out$res
res_global1 = out$res_global

#This one should be essentially the same but adjust for connection diffferences? Overall seems VERY similar comparong the two, some very slight shifts
out2 = ancombc(phyloseq = OBJ_W14, formula = "Location + Connection", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Location",
               struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5,max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res2 = out2$res
res_global2 = out2$res_global

#Append coefficients, adjusted-P values etc
tab_coef = res1$beta
tab_q = res1$q
tab_se = res1$se
tab_w = res1$W
tab_diff = res1$diff_abn
tab_p = res1$p_val

#ignore the binds for now, just export as seperate .csv's, columns get too confusing with the same names
res_global2 <- cbind(res_global2, tab_coef, tab_q)
res_global2 <- cbind(res_global2, tab_q)

res_global_tax1 = cbind(as(res_global1, "data.frame"), as(tax_table(OBJ_W14)[rownames(res_global1), ], "matrix"))
write.csv(as.data.frame(res_global_tax1), file = "ANCOM-BC_ALL_ASV_W14_Loc_GLOBAL.csv")
tab_coef = cbind(as(tab_coef, "data.frame"), as(tax_table(OBJ_W14)[rownames(tab_coef), ], "matrix"))
write.csv(as.data.frame(tab_coef), file = "ANCOM-BC_ALL_ASV_W14_Loc_Coef.csv")
tab_q = cbind(as(tab_q, "data.frame"), as(tax_table(OBJ_W14)[rownames(tab_q), ], "matrix"))
write.csv(as.data.frame(tab_q), file = "ANCOM-BC_ALL_ASV_W14_Loc_Q_adjustedP.csv")
tab_se = cbind(as(tab_se, "data.frame"), as(tax_table(OBJ_W14)[rownames(tab_se), ], "matrix"))
write.csv(as.data.frame(tab_se), file = "ANCOM-BC_ALL_ASV_W14_Loc_resid.csv")
tab_w = cbind(as(tab_w, "data.frame"), as(tax_table(OBJ_W14)[rownames(tab_w), ], "matrix"))
write.csv(as.data.frame(tab_w), file = "ANCOM-BC_ALL_ASV_W14_Loc_Wstat.csv")
tab_diff = cbind(as(tab_diff, "data.frame"), as(tax_table(OBJ_W14)[rownames(tab_diff), ], "matrix"))
write.csv(as.data.frame(tab_diff), file = "ANCOM-BC_ALL_ASV_W14_Loc_DA.csv")
tab_p = cbind(as(tab_p, "data.frame"), as(tax_table(OBJ_W14)[rownames(tab_p), ], "matrix"))
write.csv(as.data.frame(tab_p), file = "ANCOM-BC_ALL_ASV_W14_Loc_Pval.csv")

#for printing out results of the second comparison.,ignore for now
res_global_tax2 = cbind(as(res_global2, "data.frame"), as(tax_table(OBJ_W14)[rownames(res_global2), ], "matrix"))
res_global_tax2
write.csv(as.data.frame(res_global_tax2), file = "ANCOM-BC_ALL_ASV_W14_Loc-Conn.csv")

##Full CUT DATASET
out3 = ancombc(phyloseq = OBJ_W14_TRIMfull, formula = "Location", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Location",
              struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5,max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res3 = out3$res
res_global3 = out3$res_global

#Append coefficients, adjusted-P values etc
tab_coef3 = res3$beta
tab_q3 = res3$q
tab_se3 = res3$se
tab_w3 = res3$W
tab_diff3 = res3$diff_abn
tab_p3 = res3$p_val

#bind tax and print results
res_global_tax3 = cbind(as(res_global3, "data.frame"), as(tax_table(OBJ_W14_TRIMfull)[rownames(res_global3), ], "matrix"))
write.csv(as.data.frame(res_global_tax3), file = "ANCOM-BC_ALL_ASV_W14_Loc_GLOBAL_CUT.csv")
tab_coef3 = cbind(as(tab_coef3, "data.frame"), as(tax_table(OBJ_W14_TRIMfull)[rownames(tab_coef3), ], "matrix"))
write.csv(as.data.frame(tab_coef3), file = "ANCOM-BC_ALL_ASV_W14_Loc_Coef_CUT.csv")
tab_q3 = cbind(as(tab_q3, "data.frame"), as(tax_table(OBJ_W14_TRIMfull)[rownames(tab_q3), ], "matrix"))
write.csv(as.data.frame(tab_q3), file = "ANCOM-BC_ALL_ASV_W14_Loc_Q_adjustedP_CUT.csv")
tab_se3 = cbind(as(tab_se3, "data.frame"), as(tax_table(OBJ_W14_TRIMfull)[rownames(tab_se3), ], "matrix"))
write.csv(as.data.frame(tab_se3), file = "ANCOM-BC_ALL_ASV_W14_Loc_resid_CUT.csv")
tab_w3 = cbind(as(tab_w3, "data.frame"), as(tax_table(OBJ_W14_TRIMfull)[rownames(tab_w3), ], "matrix"))
write.csv(as.data.frame(tab_w3), file = "ANCOM-BC_ALL_ASV_W14_Loc_Wstat_CUT.csv")
tab_diff3 = cbind(as(tab_diff3, "data.frame"), as(tax_table(OBJ_W14_TRIMfull)[rownames(tab_diff3), ], "matrix"))
write.csv(as.data.frame(tab_diff3), file = "ANCOM-BC_ALL_ASV_W14_Loc_DA_CUT.csv")
tab_p3 = cbind(as(tab_p3, "data.frame"), as(tax_table(OBJ_W14_TRIMfull)[rownames(tab_p3), ], "matrix"))
write.csv(as.data.frame(tab_p3), file = "ANCOM-BC_ALL_ASV_W14_Loc_Pval_CUT.csv")

#Bind results sheets with OTU Taxa
res_global_tax3 = cbind(as(res_global3, "data.frame"), as(tax_table(OBJ_W14_TRIMfull)[rownames(res_global3), ], "matrix"))
res_global_tax3
#for printing out results
write.csv(as.data.frame(res_global_tax), file = "ANCOM-BC_ALL_ASV_W14_Loc_Cut.csv")
write.csv(as.data.frame(tab_coef), file = "ANCOM-BC_ALL_ASV_W14_Loc_Cut_coefficients.csv")
write.csv(as.data.frame(tab_q), file = "ANCOM-BC_ALL_ASV_W14_Loc_Cut_adjusted_P.csv")

#HALFCUT DATASET
out4 = ancombc(phyloseq = OBJ_W14_TRIMhalf, formula = "Location", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Location", 
              struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5,max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIMhalf)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#bind tax and print results
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIMhalf)[rownames(res_global4), ], "matrix"))
write.csv(as.data.frame(res_global_tax4), file = "ANCOM-BC_ALL_ASV_W14_Loc_GLOBAL_cuthalf.csv")
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIMhalf)[rownames(tab_coef4), ], "matrix"))
write.csv(as.data.frame(tab_coef4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Coef_cuthalf.csv")
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIMhalf)[rownames(tab_q4), ], "matrix"))
write.csv(as.data.frame(tab_q4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Q_adjustedP_cuthalf.csv")
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIMhalf)[rownames(tab_se4), ], "matrix"))
write.csv(as.data.frame(tab_se4), file = "ANCOM-BC_ALL_ASV_W14_Loc_resid_cuthalf.csv")
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIMhalf)[rownames(tab_w4), ], "matrix"))
write.csv(as.data.frame(tab_w4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Wstat_cuthalf.csv")
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIMhalf)[rownames(tab_diff4), ], "matrix"))
write.csv(as.data.frame(tab_diff4), file = "ANCOM-BC_ALL_ASV_W14_Loc_DA_cuthalf.csv")
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIMhalf)[rownames(tab_p4), ], "matrix"))
write.csv(as.data.frame(tab_p4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Pval_cuthalf.csv")



# MaAsLin2 ----------------------------------------------------------------
#as a comparison point with deseq and ancombc
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Maaslin2")
library(Maaslin2)
mas_1 <- Maaslin2(
  input_data = data.frame(otu_table(ps)),
  input_metadata = data.frame(sample_data(ps)),
  output = "/Users/olljt2/desktop/Maaslin2_default_output",
  min_abundance = 0.0,
  min_prevalence = 0.0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = "location",
  correction = "BH",
  standardize = FALSE,
  cores = 1)

#Indicator species - Indicspecies-------------------------------------------------------------
# first seems to be recommended that you have otus as columns instead of rows..so we need to transpose..
#starting with otu_table, lets make a new transposed one
library(indicspecies)
#Overall
otu_transposed <- data.frame(t(otu_table))
abund = otu_transposed[,2:ncol(otu_transposed)]
location = treat$Location
inv = multipatt(abund, location, func = "r.g", control = how(nperm = 9999))

#Week 14 full

#Week 14 cut

#Week 14 half-cut

####### FINAL Cut set CONNECTION ------------------------------------------------
##   ANCOM OF CONNECTION / LACK ODF CONNECTION

#### Final Cut Set ANCOM-BC LOCATION --------------------------------------------------
#FINAL half-cut DATASET
#HALFCUT DATASET
out4 = ancombc(phyloseq = OBJ_W14_TRIM, formula = "Location + Connection + Inoculum", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Location", 
              struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#bind tax and print results
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
write.csv(as.data.frame(res_global_tax4), file = "ANCOM-BC_ASV_W14_Loc_GLOBAL_FINALcut.csv")
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_coef4), ], "matrix"))
write.csv(as.data.frame(tab_coef4), file = "ANCOM-BC_ASV_W14_Loc_Coef_FINALcut.csv")
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_q4), ], "matrix"))
write.csv(as.data.frame(tab_q4), file = "ANCOM-BC_ASV_W14_Loc_Q_adjustedP_FINALcut.csv")
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_se4), ], "matrix"))
write.csv(as.data.frame(tab_se4), file = "ANCOM-BC_ASV_W14_Loc_resid_FINALcut.csv")
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_w4), ], "matrix"))
write.csv(as.data.frame(tab_w4), file = "ANCOM-BC_ASV_W14_Loc_Wstat_FINALcut.csv")
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_diff4), ], "matrix"))
write.csv(as.data.frame(tab_diff4), file = "ANCOM-BC_ASV_W14_Loc_DA_FINALcut.csv")
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_p4), ], "matrix"))
write.csv(as.data.frame(tab_p4), file = "ANCOM-BC_ASV_W14_Loc_Pval_FINALcut.csv")

## final cut Version Agglomerated to family ---------------------------------------------------
#agglomerate half cut dataset to Family for higher level look at differential abundance
#glom first
OBJ_W14_TRIM_fam <- tax_glom(OBJ_W14_TRIM,taxrank = "Family")
#ensure referring to agglomerated object
out4 = ancombc(phyloseq = OBJ_W14_TRIM_fam, formula = "Location + Connection + Inoculum", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Location", 
              struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#Export as single excel file in separate worksheets
library(openxlsx)

res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_coef4), ], "matrix"))
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_q4), ], "matrix"))
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_se4), ], "matrix"))
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_w4), ], "matrix"))
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_diff4), ], "matrix"))
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_p4), ], "matrix"))

colnames(res_global_tax4)
res_global_tax4 <- data.frame("ASV" = rownames(res_global_tax4), res_global_tax4)
tab_coef4 <- data.frame("ASV" = rownames(tab_coef4), tab_coef4)
tab_q4 <- data.frame("ASV" = rownames(tab_q4), tab_q4)
tab_se4 <- data.frame("ASV" = rownames(tab_se4), tab_se4)
tab_w4 <- data.frame("ASV" = rownames(tab_w4), tab_w4)
tab_diff4 <- data.frame("ASV" = rownames(tab_diff4), tab_diff4)
tab_p4 <- data.frame("ASV" = rownames(tab_p4), tab_p4)
colnames(res_global_tax4)

dataset_names <- list('Global' = res_global_tax4, 'Coefficients' = tab_coef4,'P_value' = tab_p4, 'Q_AdjustedP' = tab_q4,'SE_Residuals' = tab_se4,'W_stat' = tab_w4,'Differential_Abundance' = tab_diff4)
write.xlsx(dataset_names, file = 'ANCOMBC_W14_Loc_FINALcut_Family.xlsx', asTable = TRUE, firstRow = TRUE)

#bind tax and print results
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
write.csv(as.data.frame(res_global_tax4), file = "ANCOM-BC_W14_Loc_GLOBAL_FINALcut_fam.csv")
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_coef4), ], "matrix"))
write.csv(as.data.frame(tab_coef4), file = "ANCOM-BC_W14_Loc_Coef_FINALcut_fam.csv")
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_q4), ], "matrix"))
write.csv(as.data.frame(tab_q4), file = "ANCOM-BC_W14_Loc_Q_adjustedP_FINALcut_fam.csv")
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_se4), ], "matrix"))
write.csv(as.data.frame(tab_se4), file = "ANCOM-BC_W14_Loc_resid_FINALcut_fam.csv")
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_w4), ], "matrix"))
write.csv(as.data.frame(tab_w4), file = "ANCOM-BC_W14_Loc_Wstat_FINALcut_fam.csv")
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_diff4), ], "matrix"))
write.csv(as.data.frame(tab_diff4), file = "ANCOM-BC_W14_Loc_DA_FINALcut_fam.csv")
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_p4), ], "matrix"))
write.csv(as.data.frame(tab_p4), file = "ANCOM-BC_W14_Loc_Pval_FINALcut_fam.csv")




out4 = ancombc(phyloseq = OBJ_W14_TRIM, formula = "Connection", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Connection", 
              struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#bind tax and print results
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
write.csv(as.data.frame(res_global_tax4), file = "ANCOM-BC_ALL_ASV_W14_Con_GLOBAL_FINALcuthalf.csv")
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_coef4), ], "matrix"))
write.csv(as.data.frame(tab_coef4), file = "ANCOM-BC_ALL_ASV_W14_Con_Coef_FINALcuthalf.csv")
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_q4), ], "matrix"))
write.csv(as.data.frame(tab_q4), file = "ANCOM-BC_ALL_ASV_W14_Con_Q_adjustedP_FINALcuthalf.csv")
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_se4), ], "matrix"))
write.csv(as.data.frame(tab_se4), file = "ANCOM-BC_ALL_ASV_W14_Con_resid_FINALcuthalf.csv")
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_w4), ], "matrix"))
write.csv(as.data.frame(tab_w4), file = "ANCOM-BC_ALL_ASV_W14_Con_Wstat_FINALcuthalf.csv")
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_diff4), ], "matrix"))
write.csv(as.data.frame(tab_diff4), file = "ANCOM-BC_ALL_ASV_W14_Con_DA_FINALcuthalf.csv")
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_p4), ], "matrix"))
write.csv(as.data.frame(tab_p4), file = "ANCOM-BC_ALL_ASV_W14_Con_Pval_FINALcuthalf.csv")

#WIP agglomerate half cut dataset to Family for higher level look at differential abundance
#glom first
OBJ_W14_TRIM_fam <- tax_glom(OBJ_W14_TRIM,taxrank = "Family")
#ensure referring to agglomerated object
out4 = ancombc(phyloseq = OBJ_W14_TRIM_fam, formula = "Connection", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Connection", 
              struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#bind tax and print results
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(res_global4), ], "matrix"))
write.csv(as.data.frame(res_global_tax4), file = "ANCOM-BC_ALL_ASV_W14_Con_GLOBAL_FINALcuthalf_fam.csv")
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_coef4), ], "matrix"))
write.csv(as.data.frame(tab_coef4), file = "ANCOM-BC_ALL_ASV_W14_Con_Coef_FINALcuthalf_fam.csv")
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_q4), ], "matrix"))
write.csv(as.data.frame(tab_q4), file = "ANCOM-BC_ALL_ASV_W14_Con_Q_adjustedP_FINALcuthalf_fam.csv")
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_se4), ], "matrix"))
write.csv(as.data.frame(tab_se4), file = "ANCOM-BC_ALL_ASV_W14_Con_resid_FINALcuthalf_fam.csv")
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_w4), ], "matrix"))
write.csv(as.data.frame(tab_w4), file = "ANCOM-BC_ALL_ASV_W14_Con_Wstat_FINALcuthalf_fam.csv")
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_diff4), ], "matrix"))
write.csv(as.data.frame(tab_diff4), file = "ANCOM-BC_ALL_ASV_W14_Con_DA_FINALcuthalf_fam.csv")
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM_fam)[rownames(tab_p4), ], "matrix"))
write.csv(as.data.frame(tab_p4), file = "ANCOM-BC_ALL_ASV_W14_Con_Pval_FINALcuthalf_fam.csv")





# Test of final cut dataset using additional facotrs in the ancom  --------

#HALFCUT DATASET
out4 = ancombc(phyloseq = OBJ_W14_TRIM, formula = "Location + Connection + Inoculum", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Location", 
              struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#bind tax and print results
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
write.csv(as.data.frame(res_global_tax4), file = "ANCOM-BC_ALL_ASV_W14_Loc_GLOBAL_FINALcuthalf.csv")
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_coef4), ], "matrix"))
write.csv(as.data.frame(tab_coef4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Coef_FINALcuthalf.csv")
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_q4), ], "matrix"))
write.csv(as.data.frame(tab_q4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Q_adjustedP_FINALcuthalf.csv")
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_se4), ], "matrix"))
write.csv(as.data.frame(tab_se4), file = "ANCOM-BC_ALL_ASV_W14_Loc_resid_FINALcuthalf.csv")
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_w4), ], "matrix"))
write.csv(as.data.frame(tab_w4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Wstat_FINALcuthalf.csv")
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_diff4), ], "matrix"))
write.csv(as.data.frame(tab_diff4), file = "ANCOM-BC_ALL_ASV_W14_Loc_DA_FINALcuthalf.csv")
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_p4), ], "matrix"))
write.csv(as.data.frame(tab_p4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Pval_FINALcuthalf.csv")


#HALFCUT DATASET - account for all facgors and global test for conection changes
out4 = ancombc(phyloseq = OBJ_W14_TRIM, formula = "Location + Connection + Inoculum", p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,group = "Connection", 
              struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.0549, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global

#Bind results sheets with OTU Taxa
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
res_global_tax4

#Append coefficients, adjusted-P values etc
tab_coef4 = res4$beta
tab_q4 = res4$q
tab_se4 = res4$se
tab_w4 = res4$W
tab_diff4 = res4$diff_abn
tab_p4 = res4$p_val

#bind tax and print results
res_global_tax4 = cbind(as(res_global4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(res_global4), ], "matrix"))
write.csv(as.data.frame(res_global_tax4), file = "ANCOM-BC_ALL_ASV_W14_Loc_GLOBAL_FINALcuthalf.csv")
tab_coef4 = cbind(as(tab_coef4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_coef4), ], "matrix"))
write.csv(as.data.frame(tab_coef4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Coef_FINALcuthalf.csv")
tab_q4 = cbind(as(tab_q4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_q4), ], "matrix"))
write.csv(as.data.frame(tab_q4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Q_adjustedP_FINALcuthalf.csv")
tab_se4 = cbind(as(tab_se4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_se4), ], "matrix"))
write.csv(as.data.frame(tab_se4), file = "ANCOM-BC_ALL_ASV_W14_Loc_resid_FINALcuthalf.csv")
tab_w4 = cbind(as(tab_w4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_w4), ], "matrix"))
write.csv(as.data.frame(tab_w4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Wstat_FINALcuthalf.csv")
tab_diff4 = cbind(as(tab_diff4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_diff4), ], "matrix"))
write.csv(as.data.frame(tab_diff4), file = "ANCOM-BC_ALL_ASV_W14_Loc_DA_FINALcuthalf.csv")
tab_p4 = cbind(as(tab_p4, "data.frame"), as(tax_table(OBJ_W14_TRIM)[rownames(tab_p4), ], "matrix"))
write.csv(as.data.frame(tab_p4), file = "ANCOM-BC_ALL_ASV_W14_Loc_Pval_FINALcuthalf.csv")






# TEST section for ordination and permanova of css vs tss vs ander --------

#some ordinations first
library("ggplot2")
library("RColorBrewer")

#NOTE: Added to mapping file additional columns with more groupings based on: 
#Group all samples per soil columns: Soil_Column
#Group based on whether soil based subtrate (root/anode) vs Water for cathode samples: Substrate


#Note re: Rstudio export, using SVG format and dimensions 1150 x 900 results in square ordination when accounting for legend size 

#weighted TSS
NMDS_W14wTRIM <- ordinate(OBJ_W14_TRIM_tss, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W14wTRIM
#unweighted TSS
NMDS_W14uTRIM <- ordinate(OBJ_W14_TRIM_tss, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W14uTRIM

#W14 with TRIMMED data (half set, so FINAL) TSS
#Unweighted
p4uTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14uTRIM, color = "Connection", shape = "Location", label = NULL)
p4uTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14uTRIM, color = "Inoculum", shape = "Location", label = NULL)
p4uTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14uTRIM, color = "Treatment", shape = "Location", label = NULL)
p4uTRIM
p4uTRIM + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p4wTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14wTRIM, color = "Connection", shape = "Location", label = NULL)
p4wTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14wTRIM, color = "Inoculum", shape = "Location", label = NULL)
p4wTRIM <- plot_ordination(OBJ_W14_TRIM_tss, NMDS_W14wTRIM, color = "Treatment", shape = "Location", label = NULL)
p4wTRIM
p4wTRIM + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))



#CSS

#weighted CSS
NMDS_W14wTRIMcss <- ordinate(OBJ_W14_TRIM_css, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W14wTRIMcss
#unweighted CSS
NMDS_W14uTRIMcss <- ordinate(OBJ_W14_TRIM_css, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W14uTRIMcss

#W14 with TRIMMED data (half set, so FINAL) CSS
#Unweighted
p4uTRIMcss <- plot_ordination(OBJ_W14_TRIM_css, NMDS_W14uTRIMcss, color = "Connection", shape = "Location", label = NULL)
p4uTRIMcss <- plot_ordination(OBJ_W14_TRIM_css, NMDS_W14uTRIMcss, color = "Inoculum", shape = "Location", label = NULL)
p4uTRIMcss <- plot_ordination(OBJ_W14_TRIM_css, NMDS_W14uTRIMcss, color = "Treatment", shape = "Location", label = NULL)
p4uTRIMcss
p4uTRIMcss + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p4wTRIMcss <- plot_ordination(OBJ_W14_TRIM_css, NMDS_W14wTRIMcss, color = "Connection", shape = "Location", label = NULL)
p4wTRIMcss <- plot_ordination(OBJ_W14_TRIM_css, NMDS_W14wTRIMcss, color = "Inoculum", shape = "Location", label = NULL)
p4wTRIMcss <- plot_ordination(OBJ_W14_TRIM_css, NMDS_W14wTRIMcss, color = "Treatment", shape = "Location", label = NULL)
p4wTRIMcss
p4wTRIMcss + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))


#Anderson log

#weighted Anderson log
NMDS_W14wTRIManderson <- ordinate(OBJ_W14_TRIM_anderson, "NMDS", distance = "unifrac", weighted = TRUE, parallel = TRUE)
NMDS_W14wTRIManderson
#unweighted Anderson log
NMDS_W14uTRIManderson <- ordinate(OBJ_W14_TRIM_anderson, "NMDS", distance = "unifrac", weighted = FALSE, parallel = TRUE)
NMDS_W14uTRIManderson

#W14 with TRIMMED data (half set, so FINAL) Anderson log
#Unweighted
p4uTRIManderson <- plot_ordination(OBJ_W14_TRIM_anderson, NMDS_W14uTRIManderson, color = "Connection", shape = "Location", label = NULL)
p4uTRIManderson <- plot_ordination(OBJ_W14_TRIM_anderson, NMDS_W14uTRIManderson, color = "Inoculum", shape = "Location", label = NULL)
p4uTRIManderson <- plot_ordination(OBJ_W14_TRIM_anderson, NMDS_W14uTRIManderson, color = "Treatment", shape = "Location", label = NULL)
p4uTRIManderson
p4uTRIManderson + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))
#Weighted
p4wTRIManderson <- plot_ordination(OBJ_W14_TRIM_anderson, NMDS_W14wTRIManderson, color = "Connection", shape = "Location", label = NULL)
p4wTRIManderson <- plot_ordination(OBJ_W14_TRIM_anderson, NMDS_W14wTRIManderson, color = "Inoculum", shape = "Location", label = NULL)
p4wTRIManderson <- plot_ordination(OBJ_W14_TRIM_anderson, NMDS_W14wTRIManderson, color = "Treatment", shape = "Location", label = NULL)
p4wTRIManderson
p4wTRIManderson + theme_grey() + theme(text = element_text(size = 14)) + geom_point(size = 3.5) + scale_color_manual(values = c("#177BB5","#56B4E9","#BF8300","#E09900","#008F47","#00B85C","#141414","#7A7A7A"))

#PERMANOVA TSS
library(vegan)
#Cleaned up with FINAL cut dataset
##make distance matricies (these are weighted and unweighted unifrac. can also use 'bray')
OBJ1_W14perm_W <- distance(OBJ_W14_TRIM_tss, "wunifrac")
OBJ1_W14perm_U <- distance(OBJ_W14_TRIM_tss, "unifrac")

#Set up variable objects
Location14 <- get_variable(OBJ_W14_TRIM_tss, "Location")
#By connection
Connection14 <- get_variable(OBJ_W14_TRIM_tss, "Connection")
#By Inoculum 
Inoculum14 <- get_variable(OBJ_W14_TRIM_tss, "Inoculum")

#Week 14
#compare adonis2
W14_3Group_ado2_U = adonis2(OBJ1_W14perm_U ~ Location14 * Connection14 * Inoculum14, permutations = 9999)
W14_3Group_ado2_U

#compare adonis2
W14_3Group_ado2_W = adonis2(OBJ1_W14perm_W ~ Location14 * Connection14 * Inoculum14, permutations = 9999)
W14_3Group_ado2_W

#adonis2 test without ordering
W14_3Group_ado2_W_TEST = adonis2(OBJ1_W14perm_W ~ Location14 * Connection14 * Inoculum14, permutations = 9999, by = NULL)
W14_3Group_ado2_W_TEST

#PERMANOVA CSS
library(vegan)
#Cleaned up with FINAL cut dataset
##make distance matricies (these are weighted and unweighted unifrac. can also use 'bray')
OBJ1_W14perm_Wcss <- distance(OBJ_W14_TRIM_css, "wunifrac")
OBJ1_W14perm_Ucss <- distance(OBJ_W14_TRIM_css, "unifrac")

#Set up variable objects
Location14css <- get_variable(OBJ_W14_TRIM_css, "Location")
#By connection
Connection14css <- get_variable(OBJ_W14_TRIM_css, "Connection")
#By Inoculum 
Inoculum14css <- get_variable(OBJ_W14_TRIM_tss, "Inoculum")

#Week 14
#compare adonis2
W14_3Group_ado2_Ucss = adonis2(OBJ1_W14perm_Ucss ~ Location14css * Connection14css * Inoculum14css, permutations = 9999)
W14_3Group_ado2_Ucss

#compare adonis2
W14_3Group_ado2_Wcss = adonis2(OBJ1_W14perm_Wcss ~ Location14css * Connection14css * Inoculum14css, permutations = 9999)
W14_3Group_ado2_Wcss

#adonis2 test without ordering
W14_3Group_ado2_W_TESTcss = adonis2(OBJ1_W14perm_Wcss ~ Location14css * Connection14css * Inoculum14css, permutations = 9999, by = NULL)
W14_3Group_ado2_W_TESTcss

#PERMANOVA anderson log
library(vegan)
#Cleaned up with FINAL cut dataset
##make distance matricies (these are weighted and unweighted unifrac. can also use 'bray')
OBJ1_W14perm_Wanderson <- distance(OBJ_W14_TRIM_anderson, "wunifrac")
OBJ1_W14perm_Uanderson <- distance(OBJ_W14_TRIM_anderson, "unifrac")

#Set up variable objects
Location14anderson <- get_variable(OBJ_W14_TRIM_anderson, "Location")
#By connection
Connection14anderson <- get_variable(OBJ_W14_TRIM_anderson, "Connection")
#By Inoculum 
Inoculum14anderson <- get_variable(OBJ_W14_TRIM_anderson, "Inoculum")

#compare adonis2
W14_3Group_ado2_Wanderson = adonis2(OBJ1_W14perm_Wanderson ~ Location14anderson * Connection14anderson * Inoculum14anderson, permutations = 9999)
W14_3Group_ado2_Wanderson

#adonis2 test without ordering
W14_3Group_ado2_W_TESTanderson = adonis2(OBJ1_W14perm_Wanderson ~ Location14anderson * Connection14anderson * Inoculum14anderson, permutations = 9999, by = NULL)
W14_3Group_ado2_W_TESTanderson



### quick dirty test for t-SNE plotting
library(vegan)
library(microbiome)
library(Rtsne) # Load package

method <- "tsne"
trans <- "hellinger"
distance <- "euclidean"

# Distance matrix for samples
ps <- microbiome::transform(OBJ1_exp, trans)

# Calculate sample similarities
dm <- vegdist(otu_table(OBJ1_exp), distance)

# Run TSNE
tsne_out <- Rtsne(dm, dims = 2) 
proj <- tsne_out$Y
rownames(proj) <- rownames(otu_table(OBJ1_exp))

library(ggplot2)
p <- plot_landscape(proj, legend = T, size = 1) 
print(p)


# anocom2 adaptation ------------------------------------------------------

#Quick test for adapting old ANCOMBC script to ANCOMBC2 
#at minimum need to double check that prv_cut in v2 is is equal to zero_cut in v1, in anycase zero_cut no longer exists so that feature must now either be an unchanging default or have it's name changed

setwd("~/Documents/University/Analysis/PMFC_18/2020 rerun outputs/Format for phyloseq")
library(phyloseq)
library(ape)
otu_table <- as.data.frame(read.csv("raw_readmap.csv", header = TRUE,row.names = "OTU_ID"))
taxmat <- as.matrix(read.csv("tax_table.csv", row.names = 1, header = TRUE))
treat <- as.data.frame(read.csv("ANCOM_mapping_file.csv", row.names = 1, header = TRUE))
TREE <- read.tree("rooted_tree.nwk")
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
TREAT = sample_data(treat)
OBJ1 = phyloseq(OTU,TAX,TREAT,TREE)
library(magrittr)
OBJ1 <- OBJ1 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
    Family  != "mitochondria" &
    Class   != "Chloroplast"
  )
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
#Half cut that retains Pseudomonas (final paper#1 data set)
OBJ_Overall_TRIM <- subset_samples(OBJ1_exp, Treatment_Half_Trim == "Retain")
OBJ_W0 <- subset_samples(OBJ1_exp, Week == "Zero")
OBJ_W0_TRIM <- subset_samples(OBJ_W0, Treatment_Half_Trim == "Retain")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")
OBJ_W14_TRIM <- subset_samples(OBJ_W14, Treatment_Half_Trim == "Retain")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ANCOMBC")

library(ANCOMBC)
library(microbiome)
library(dplyr)
library(ggplot2)
library(data.table)

#### Final Cut Dataset import for ANCOM-BC2 test --------------------------------------------------
#FINAL half-cut DATASET - DONT RUN THIS ONE WAS JUST FOR INITIAL TESTING AND KEEPING FOR REFERENCE< 
#START AT THE BELOW VERIOSNS FOR EAACH TAXAT LEVEL AND FACTOR
#HALFCUT DATASET, looking at location, controlling for connection/inoculum
out4 = ancombc2(data = OBJ_W14_TRIM, fix_formula = "Location + Connection + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Location", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)
res4 = out4$res
res_global4 = out4$res_global



# Start HERE for fixing updated anocom2bc with FINAL cut and trim data sets for version ----------------------------------------------------


#CONNECTION, using unconnnected as base

#So i'll need 3 (4 if doing species too) sets of results for each group and taxa level, and each of those (except connection), will result in a primary and a global

#Connection at FAMILY
output_CON_FAM = ancombc2(data = OBJ_W14_TRIM, tax_level = "Family" , fix_formula = "Connection + Location + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Connection", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_CON_FAM = output_CON_FAM$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_CON_FAM %>%
    data.table(caption = "ANCOM-BC2 Primary Results Connection Family")
write.csv(res_prim_CON_FAM,"ANCOM-BC2 Primary Results_Connection_Family.csv", row.names = TRUE)


#Connection at GENUS
output_CON_GEN = ancombc2(data = OBJ_W14_TRIM, tax_level = "Genus" , fix_formula = "Connection + Location + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Connection", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_CON_GEN = output_CON_GEN$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_CON_GEN %>%
    data.table(caption = "ANCOM-BC2 Primary Results_Connection_Genus")
write.csv(res_prim_CON_GEN,"ANCOM-BC2 Primary Results_Connection_Genus.csv", row.names = TRUE)


#Connection at SPECIES
output_CON_Spe = ancombc2(data = OBJ_W14_TRIM, tax_level = "Species" , fix_formula = "Connection + Location + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Connection", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_CON_Spe = output_CON_Spe$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_CON_Spe %>%
    data.table(caption = "ANCOM-BC2 Primary Results_Conn_Species")
write.csv(res_prim_CON_Spe,"ANCOM-BC2 Primary Results_Conn_Species.csv", row.names = TRUE)


#Connection  at ASV
output_CON_ASV = ancombc2(data = OBJ_W14_TRIM, tax_level = NULL , fix_formula = "Connection + Location + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Connection", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_CON_ASV = output_CON_ASV$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_CON_ASV %>%
    data.table(caption = "ANCOM-BC2 Primary Results_Conn_ASV")
write.csv(res_prim_CON_ASV,"ANCOM-BC2 Primary Results_Conn_ASV.csv", row.names = TRUE)




#LOCATION
#Location, controlling for inoculum and connection,

#Location: controlling for conn/in at FAMILY level and using roots as base
output_LOC_FAM = ancombc2(data = OBJ_W14_TRIM, tax_level = "Family" , fix_formula = "Location + Connection + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Location", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_LOC_FAM = output_LOC_FAM$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_LOC_FAM %>%
    data.table(caption = "ANCOM-BC2 Primary Results_LOC_FAM")
write.csv(res_prim_LOC_FAM,"ANCOM-BC2 Primary Results_LOC_FAM.csv", row.names = TRUE)

res_global_LOC_FAM = output_LOC_FAM$res_global %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_global_LOC_FAM %>%
    data.table(caption = "ANCOM-BC2 Global Results_LOC_FAM")
write.csv(res_global_LOC_FAM,"ANCOM-BC2 Global Results_LOC_FAM.csv", row.names = TRUE)

#Location: controlling for conn/in at GENUS level and using uninoculated as base
output_LOC_GEN = ancombc2(data = OBJ_W14_TRIM, tax_level = "Genus" , fix_formula = "Location + Connection + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Location", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_LOC_GEN = output_LOC_GEN$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_LOC_GEN %>%
    data.table(caption = "ANCOM-BC2 Primary Results_LOC_GEN")
write.csv(res_prim_LOC_GEN,"ANCOM-BC2 Primary Results_LOC_GEN.csv", row.names = TRUE)

res_global_LOC_GEN = output_LOC_GEN$res_global %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_global_LOC_GEN %>%
    data.table(caption = "ANCOM-BC2 Global Results_LOC_GEN")
write.csv(res_global_LOC_GEN,"ANCOM-BC2 Global Results_LOC_GEN.csv", row.names = TRUE)

#Location: controlling for conn/in at SPECIES level and using uninoculated as base
output_LOC_SPE = ancombc2(data = OBJ_W14_TRIM, tax_level = "Species" , fix_formula = "Location + Connection + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Location", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_LOC_SPE = output_LOC_SPE$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_LOC_SPE %>%
    data.table(caption = "ANCOM-BC2 Primary Results_LOC_SPE")
write.csv(res_prim_LOC_SPE,"ANCOM-BC2 Primary Results_LOC_SPE.csv", row.names = TRUE)

res_global_LOC_SPE = output_LOC_SPE$res_global %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_global_LOC_SPE %>%
    data.table(caption = "ANCOM-BC2 Global Results_LOC_SPE")
write.csv(res_global_LOC_SPE,"ANCOM-BC2 Global Results_LOC_SPE.csv", row.names = TRUE)

#Location: controlling for conn/in at ASV level and using uninoculated as base
output_LOC_ASV = ancombc2(data = OBJ_W14_TRIM, tax_level = NULL , fix_formula = "Location + Connection + Inoculum", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Location", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_LOC_ASV = output_LOC_ASV$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_LOC_ASV %>%
    data.table(caption = "ANCOM-BC2 Primary Results_LOC_ASV")
write.csv(res_prim_LOC_ASV,"ANCOM-BC2 Primary Results_LOC_ASV.csv", row.names = TRUE)

res_global_LOC_ASV = output_LOC_ASV$res_global %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_global_LOC_ASV %>%
    data.table(caption = "ANCOM-BC2 Global Results_LOC_ASV")
write.csv(res_global_LOC_ASV,"ANCOM-BC2 Global Results_LOC_ASV.csv", row.names = TRUE)





#INOCULUM

#Inoculum at FAMILY controlling for location and connection, uninoculated as baseline

output_INO_FAM = ancombc2(data = OBJ_W14_TRIM, tax_level = "Family" , fix_formula = "Inoculum + Location + Connection", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Inoculum", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_INO_FAM = output_INO_FAM$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_INO_FAM %>%
    data.table(caption = "ANCOM-BC2 Primary Results_INO_FAM")
write.csv(res_prim_INO_FAM,"ANCOM-BC2 Primary Results_INO_FAM.csv", row.names = TRUE)

res_global_INO_FAM = output_INO_FAM$res_global %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_global_INO_FAM %>%
    data.table(caption = "ANCOM-BC2 Global Results_INO_FAM")
write.csv(res_global_INO_FAM,"ANCOM-BC2 Global Results_INO_FAM.csv", row.names = TRUE)

#Inoculum at GENUS controlling for location and connection, uninoculated as baseline

output_INO_GEN = ancombc2(data = OBJ_W14_TRIM, tax_level = "Genus" , fix_formula = "Inoculum + Location + Connection", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Inoculum", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_INO_GEN = output_INO_GEN$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_INO_GEN %>%
    data.table(caption = "ANCOM-BC2 Primary Results_INO_GEN")
write.csv(res_prim_INO_GEN,"ANCOM-BC2 Primary Results_INO_GEN.csv", row.names = TRUE)

res_global_INO_GEN = output_INO_GEN$res_global %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_global_INO_GEN %>%
    data.table(caption = "ANCOM-BC2 Global Results_INO_GEN")
write.csv(res_global_INO_GEN,"ANCOM-BC2 Global Results_INO_GEN.csv", row.names = TRUE)

#Inoculum at SPECIES controlling for location and connection, uninoculated as baseline

output_INO_SPE = ancombc2(data = OBJ_W14_TRIM, tax_level = "Species" , fix_formula = "Inoculum + Location + Connection", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Inoculum", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_INO_SPE = output_INO_SPE$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_INO_SPE %>%
    data.table(caption = "ANCOM-BC2 Primary Results_INO_SPE")
write.csv(res_prim_INO_SPE,"ANCOM-BC2 Primary Results_INO_SPE.csv", row.names = TRUE)

res_global_INO_SPE = output_INO_SPE$res_global %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_global_INO_SPE %>%
    data.table(caption = "ANCOM-BC2 Global Results_INO_SPE")
write.csv(res_global_INO_SPE,"ANCOM-BC2 Global Results_INO_SPE.csv", row.names = TRUE)

#Inoculum at ASV controlling for location and connection, uninoculated as baseline

output_INO_ASV = ancombc2(data = OBJ_W14_TRIM, tax_level = NULL , fix_formula = "Inoculum + Location + Connection", p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,group = "Inoculum", 
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = TRUE)

res_prim_INO_ASV = output_INO_ASV$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim_INO_ASV %>%
    data.table(caption = "ANCOM-BC2 Primary Results_INO_ASV")
write.csv(res_prim_INO_ASV,"ANCOM-BC2 Primary Results_INO_ASV.csv", row.names = TRUE)

res_global_INO_ASV = output_INO_ASV$res_global %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_global_INO_ASV %>%
    data.table(caption = "ANCOM-BC2 Global Results_INO_ASV")
write.csv(res_global_INO_ASV,"ANCOM-BC2 Global Results_INO_ASV.csv", row.names = TRUE)



# Use deseq to get basemean for each taxa level compared in results -------
library(ggplot2)
library("DESeq2")
#Start above with import for ancombx2 and experimental/cut subesetting etc first

#Agglomerate at desired taxa level (if you want to, otherwise proceed to and will get individual ASV changes)
OBJ_W14_TRIM # use this one for getting individual ASV values
OBJ1_exp_Species <- tax_glom(OBJ_W14_TRIM,taxrank = "Species")
OBJ1_exp_Genus <- tax_glom(OBJ_W14_TRIM,taxrank = "Genus")
OBJ1_exp_Family <- tax_glom(OBJ_W14_TRIM,taxrank = "Family")

#And create subset for use later when assigning taxa (NOT for analysis, just to bind taxa to results)
OBJ_Genus_W14 <- subset_samples(OBJ1_exp_Genus, Week == "Fourteen")


#Week 0 one, for getting baseMean as a comparison point
OBJ1_expW1_Family <- tax_glom(OBJ_W0_TRIM,taxrank = "Family")
diagdds = phyloseq_to_deseq2(OBJ_W0_TRIM, ~ Location + Connection) #Re: order of factors here, this would be testing for the effect of location, controlling for connection



##Import to deseq and order factors to be compared (overall, not agglomerated)
diagdds = phyloseq_to_deseq2(OBJ_W14_TRIM, ~ Location + Connection) #Re: order of factors here, this would be testing for the effect of location, controlling for connection
##Import agglomerated species
diagdds = phyloseq_to_deseq2(OBJ1_exp_Species, ~ Location + Connection) #Re: order of factors here, this would be testing for the effect of location, controlling for connection
##Import agglomerated genus
diagdds = phyloseq_to_deseq2(OBJ1_exp_Genus, ~ Location + Connection) #Re: order of factors here, this would be testing for the effect of location, controlling for connection
#import agglomerated family
diagdds = phyloseq_to_deseq2(OBJ1_exp_Family , ~ Location + Connection) #Re: order of factors here, this would be testing for the effect of location, controlling for connection


## Compare connection with lack of connection (only at wk 14) with trimmed data fir final paper tables looking at
diagdds = phyloseq_to_deseq2(OBJ1_exp_Family , ~ Connection + Location) #Re: order of factors here, this would be testing for the effect of location, controlling for connection


diagdds$Location <- relevel(diagdds$Connection, ref = "Unconnected") # sets the reference point, baseline or control to be compared against

#Subset Week 14
diagdds <- diagdds[ , diagdds$Week == "Fourteen" ]
#check that subset was done as expected
as.data.frame( colData(diagdds) )

#Because of the way the data is nested between groups need to do some finangling to allow the model to distinguish between them correctly
#Have added a new column that groups samples by the soil column they were in (without location data, so e.g W14UP4 for all three cathode/anode/root sampeles)

#Run model and factors
diagdds = DESeq(diagdds, test = "Wald", fitType = "parametric")
res = results(diagdds, cooksCutoff = FALSE)
res # print out results

#Shouldnt need these anymore for agglomerated, but retain for overall
OBJ1_exp <- subset_samples(OBJ1, Experiment == "Y")
OBJ1_exp <- subset_samples(OBJ1, Treatment_Half_Trim == "Retain")
OBJ_W14 <- subset_samples(OBJ1_exp, Week == "Fourteen")

#Bind taxonomy to results
res = cbind(as(res, "data.frame"), as(tax_table(OBJ1_exp_Family)[rownames(res), ], "matrix"))
res

#Export .csv
write.csv(as.data.frame(res), 
          file = "DESeq2_W14_ConnectionBaseMeans.csv")

#For Agglomerated
#Different Comparison Direction Sheets
sigtabC_U = results(diagdds, contrast = c("Connection","Connected","Unconnected")) #a positive number here should represent an Increase in Connected FROM Unconnected
sigtabU_C = results(diagdds, contrast = c("Connection","Unconnected","Connected")) #postive here should be switched, so a postive = Increase in Unconnected FROM Connected
# So i think A should be the one to use...feels more logical to look at change towards connection (and makes discussion of electroactive enrichemnet easier)

#Bind results sheets with OTU Taxa
sigtab_taxC_U = cbind(as(sigtabC_U, "data.frame"), as(tax_table(OBJ_Genus_W14)[rownames(sigtabC_U), ], "matrix"))
sigtab_taxC_U

sigtab_taxU_C = cbind(as(sigtabU_C, "data.frame"), as(tax_table(OBJ_Genus_W14)[rownames(sigtabU_C), ], "matrix"))
sigtab_taxU_C

#Export .csv
write.csv(as.data.frame(sigtab_taxC_U), 
          file = "DESeq2_resultsC_U.csv")
write.csv(as.data.frame(sigtab_taxU_C), 
          file = "DESeq2_resultsU_C.csv")

