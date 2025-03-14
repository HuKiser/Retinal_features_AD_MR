### Script for analysing the causal risk factors of Alzheimer's disease identified from fundal imaging of the eye using bidirectional Mendelian Randomization (MR)
### TwoSampleMR package was applied for the estiamting the bidirectional effect between Alzheimer's disease and ocular features,
### For more information on commands/additional options see: https://mrcieu.github.io/TwoSampleMR/
## For simplicity, the code of bidirectional MR analysis between Alzheimer's disease and RNFL are provided here
######################################################################################################################################
# Load packages
######################################################################################################################################
rm(list=ls(all=TRUE)) #empties your R environment

#Load packages 
install.packages("devtools")
library(devtools)
install_github("mrcieu/ieugwasr")
library(ieugwasr)
install_github("MRCIEU/TwoSampleMR") 
library(TwoSampleMR)
install.packages("ggplot2")
library(ggplot2)


setwd("PATH_TO_FILE")

########################################################################################
########### Estimtation of the causal effect of AD on ocular features ####################
################## AD on RNFL thickness ################################################

#### For the direction of MR casual effect AD to Ocular features, we extracted GWAS summary data of Alzheimer's disease (Kunkle et al, Nature Genetics, 2019) 
## from IEU OpenGWAS project (ID: "ieu-b-2") for the list of variants that showed best level of association after meta analysis of Stage 1 and Stage 2 
## (from the published table 1: https://www.nature.com/articles/s41588-019-0358-2/tables/1). I also copied the column of effect allele from that published table 
## as this column was missing in the summary data of AD in IEU OpenGWAS project.

############# Read exposure data (AD)
ADx<-read_exposure_data(filename = "Kunkle_publised_gwas_AD_summary_data_with_maf.txt",
                              snp_col = "MarkerName",
                              beta_col = "Beta",
                              se_col = "SE",
                              effect_allele_col = "Effect_allele",
                              other_allele_col = "Non_Effect_allele",
                              eaf_col = "MAF",
                              pval_col = "Pvalue",
                              chr_col = "Chromosome",
                              pos_col = "Position"
)
dim(ADx) # 21 SNPs



#### these columns added to perform Steiger filtering for binary exposure (AD) 
ADx$exposure<-"AD"
ADx$prevalence.exposure<-0.1
ADx$ncase.exposure<-21982
ADx$ncontrol.exposure<-41944
ADx$units.exposure<-"log odds"


### no clumping performed as these variants were reported genome wide significant for AD (P<=5E-08)

## Source of RNFL summary data: https://www.ebi.ac.uk/gwas/efotraits/OBA_2050111 
### read outcome data (RNFL)  
RNFL<-read_outcome_data(snps = ADx$SNP, filename = "GCST90014266_buildGRCh37.tsv", sep = "\t",
                        snp_col = "variant_id",
                        beta_col = "beta",
                        se_col = "standard_error",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele",
                        eaf_col = "effect_allele_frequency",
                        pval_col = "p_value",
                        chr_col = "chromosome",
                        pos_col = "base_pair_location"
)


#### these columns added to perform Steiger filtering for the outcome (RNFL)
RNFL$outcome<-"RNFL"
RNFL$units.outcome<-rep("SD", nrow(RNFL))
RNFL$samplesize.outcome<-rep(31434, nrow(RNFL))

## checking if any of the IV are significantly associated to outcome of interst
length(which(RNFL$pvalue.outcome<=5E-8)) 

#### Harmonization ####################################
dat1 <- harmonise_data(
  exposure_dat = ADx, 
  outcome_dat = RNFL
)

#### Main MR analysis ###############################
mr_results<-mr(dat1, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

mr_results

###### adjusting the beta by multiplying 0.693 for MR with binary exposure ################
mr_results$cbeta<-mr_results$b*0.693
mr_results$cse<-mr_results$se*0.693
mr_results$lcl<- mr_results$cbeta-1.96*mr_results$cse
mr_results$ucl<- mr_results$cbeta+1.96*mr_results$cse
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$cbeta,mr_results$lcl, mr_results$ucl,mr_results$pval)
results
names(results)<-c("Outcome", "Exposure", "N_SNPs", "MR_Methods", "Beta", "CI_Lower_Limit", "CI_Upper_Limit", "P_value")
results


#### Measuring the evidence of heterogeneity (significant heterogeneity in the SNP-exposure effects might be suggestive of pleiotropy)
het <- mr_heterogeneity(dat1)
het

########## measuring the evidence of horizontal pleiotropy (MR egger intercept), A significant intercept suggests signifciant pleiotropy.
pleio <- mr_pleiotropy_test(dat1)
pleio

######## confidence interval of MR egger inercept
pleio$egger_intercept-1.96*pleio$se
pleio$egger_intercept+1.96*pleio$se

########### performing single SNP analysis:

##### default single SNP analyses gives the wald ratio
res_single <- mr_singlesnp(dat1)

####### leave one out analyses - by defalut uses IVW
res_loo <- mr_leaveoneout(dat1)

####### performing  MR Steiger filtering for additional sensitivity analysis

##### steiger filtering
dat1_steiger<-steiger_filtering(dat1)
dat1_AD_RNFL_steiger<-subset(dat1_steiger, dat1_steiger$steiger_dir==TRUE)

#### MR analysis after Steiger filtering
mr_results_steiger<-mr(dat1_AD_RNFL_steiger, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
mr_results_steiger

####################################################### Generating MR plots ###################################################
### Generating a Scatter plot of main MR analysis
p1<-mr_scatter_plot(mr_results,dat1)
p1[[1]]
ggsave(p1[[1]], file="AD_RNFL_scatter_plot.pdf", width=7, height=7)

### Generating a forest plot of each of the SNP effects
p2<-mr_forest_plot(res_single)
p2[[1]]
ggsave(p2[[1]], file="AD_RNFL_forest_plot.pdf", width=7, height=7)

# Generating a funnel plot to check asymmetry
p3<-mr_funnel_plot(res_single)
p3[[1]]
ggsave(p3[[1]], file="AD_RNFL_main_funnel_plot.pdf", width=7, height=7)

# generating a leave one out plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
p4<-mr_leaveoneout_plot(res_loo)
p4[[1]]
ggsave(p4[[1]], file="AD_RNFL_main_leaveoneout_plot.pdf", width=7, height=7)


########################################################################################
########### Estimtation of the causal effect of ocular features on AD ####################
################## RNFL thickness on AD ################################################


## Source of RNFL summary data: https://www.ebi.ac.uk/gwas/efotraits/OBA_2050111 
### read exposure data (RNFL)
RNFL<-read_exposure_data(filename = "GCST90014266_buildGRCh37.tsv", sep = "\t",
                         snp_col = "variant_id",
                         beta_col = "beta",
                         se_col = "standard_error",
                         effect_allele_col = "effect_allele",
                         other_allele_col = "other_allele",
                         eaf_col = "effect_allele_frequency",
                         pval_col = "p_value",
                         chr_col = "chromosome",
                         pos_col = "base_pair_location"
)

#### these columns added to perform Steiger filtering for the exposure (RNFL)
RNFL$exposure<-"RNFL"
RNFL$units.exposure<-rep("SD", nrow(RNFL))
RNFL$samplesize.exposure<-rep(31434, nrow(RNFL))
head(RNFL)

##### selecting the SNPs associated with RNFL with p<5E-08
RNFL1<-subset(RNFL, RNFL$pval.exposure<=5E-8)

#### clumping to get independantly associated SNPs with RNFL thickness
RNFL2<- clump_data(RNFL1)

#### Read outcome data using the IEU OpenGWAS project code: "ieu-b-2"
AD1<-extract_outcome_data(snps = RNFL2$SNP, outcomes = "ieu-b-2")

#### checking if any of the IV are significantly associated to outcome of interst
length(which(AD1$pvalue.outcome<=5E-8)) 

#### these columns are newly added to perform Steiger filtering for binary outcome (AD) 
AD1$prevalence.outcome<-0.1
AD1$ncase.outcome<-21982
AD1$ncontrol.outcome<-41944
AD1$units.outcome<-"log odds"

##################### extracting the EAF column from the european reference pannel (1000Genome project) as this column is missing in IEU OpendGWAS dataset
######## this step is only required if existing summary data has missing effect allele (EAF) column 
### install biomaRt package for extracting missing column from 1000Genome project
BiocManager::install("biomaRt")
library(biomaRt)

#### Connect to Ensembl database using biomaRt
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "grch37.ensembl.org")  # Using GRCh37 version

AD1<-AD1[, -which(names(AD1)== "eaf.outcome")]

### List of SNPs of AD 
snp_list <- AD1$SNP

### Fetch SNP information including effect allele frequency (EUR population)
snp_data <- getBM(
  attributes = c("refsnp_id", "allele", "minor_allele", "minor_allele_freq"),
  filters = "snp_filter",
  values = snp_list,
  mart = ensembl
)

snp_data

##### merging a subset of outcome data and reference data
data1<-AD1[, c("SNP","effect_allele.outcome", "other_allele.outcome", "beta.outcome")]
data2<-snp_data[, c("refsnp_id", "minor_allele", "minor_allele_freq")]
data12<-merge(data1, data2, by.x = "SNP", by.y = "refsnp_id")
data12$eaf.outcome<-"NA"

#### converting minor allele frequency (MAF) to effect allele frequnency (EAF) for the subset of outcome data
data123 <- within(data12, {
  # Check if Effect_Allele matches Minor_Allele
  match_effect_minor <- effect_allele.outcome == minor_allele
  
  # If they do not match, switch the  minor allele to effect allele
  minor_allele <- ifelse(match_effect_minor, minor_allele, effect_allele.outcome)
  
  # Update Effect_Allele_Frequency from 1-Minor_Allele_Frequency if they didn't match
  eaf.outcome <- ifelse(match_effect_minor, eaf.outcome, 1-minor_allele_freq)
  # Update Effect_Allele_Frequency from Minor_Allele_Frequency if they match
  eaf.outcome <- ifelse(!match_effect_minor, eaf.outcome, minor_allele_freq)
})

### extracting the SNPs and EAF from the converted dataset and merging the eaf column with main outcome data
eafdata<-data123[, c("SNP", "eaf.outcome")]
AD2<-merge(AD1, eafdata, by = "SNP")

#################################################################################################################

#### Harmonization ####################################
dat2 <- harmonise_data(
  exposure_dat = RNFL2, 
  outcome_dat = AD2
)

#### main MR analysis
mr_results<-mr(dat2, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
mr_results

#Converting the estimatd beta to odds ratio and estimating 95% confidence interval
or<-generate_odds_ratios(mr_results)
or
orresults<-cbind.data.frame(or$outcome,or$exposure, or$nsnp, or$method, or$or, or$or_lci95, or$or_uci95, or$pval)
orresults
names(orresults)<-c("Outcome", "Exposure", "N_SNPs", "MR_Methods", "OR", "CI_Lower_Limit", "CI_Upper_Limit", "P_value")
orresults

# Runing some sensitivity analyses
# measuring the evidence of heterogeneity in the genetic effects
het <- mr_heterogeneity(dat2)


# Measuring the evidence of horizontal pleiotropy
pleio <- mr_pleiotropy_test(dat2)

########### performing single SNP analysis:

##### default single SNP analyses gives the wald ratio
res_single <- mr_singlesnp(dat2)

####### leave one out analyses - by defalut uses IVW
res_loo <- mr_leaveoneout(dat2)

####### performing  MR Steiger filtering for additional sensitivity analysis

##### steiger filtering
dat2_steiger<-steiger_filtering(dat2)
dat1_RNFL_AD_steiger<-subset(dat2_steiger, dat2_steiger$steiger_dir==TRUE)

#### MR analysis after Steiger filtering
mr_results_steiger<-mr(dat2_RNFL_AD_steiger, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
mr_results_steiger

####################################################### Generating MR plots ###################################################
### Generating a Scatter plot of main MR analysis
p1<-mr_scatter_plot(mr_results,dat2)
p1[[1]]
ggsave(p1[[1]], file="RNFL_AD_scatter_plot.pdf", width=7, height=7)

### Generating a forest plot of each of the SNP effects
p2<-mr_forest_plot(res_single)
p2[[1]]
ggsave(p2[[1]], file="RNFL_AD_forest_plot.pdf", width=7, height=7)

# Generating a funnel plot to check asymmetry
p3<-mr_funnel_plot(res_single)
p3[[1]]
ggsave(p3[[1]], file="RNFL_AD_funnel_plot.pdf", width=7, height=7)

# generating a leave one out plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
p4<-mr_leaveoneout_plot(res_loo)
p4[[1]]
ggsave(p4[[1]], file="RNFL_AD_leaveoneout_plot.pdf", width=7, height=7)


################ F statistic #########################################

### selecting the SNPs applied in the main MR analysis
dat_mr<-subset(dat2, dat2$mr_keep==TRUE)

#### defining the vector values from the exposure data

eaf <- dat_mr[,"eaf.exposure"]
b <- dat_mr[,"beta.exposure"]
se <- dat_mr[,"se.exposure"]
p <- dat_mr[,"pval.exposure"]
n <- dat_mr[,"samplesize.exposure"]
snp <- dat_mr[,"SNP"]

# Converting EAF to MAF where necessary
maf <- ifelse(eaf > 0.5, 1 - eaf, eaf)

# Calculating per SNP R2
r2 <- (2 * b^2 * maf * (1 - maf)) / ((2 * b^2 * maf * (1 - maf)) + (se^2 * (2 * n) * maf * (1 - maf)))

# Individual F-stats
k <- 1
F <- r2 * (n - 1 - k) / ((1 - r2) * k)

# Overall R2 and F-stats
k <- length(snp)
all_r2 <- sum(r2)
all_F <- all_r2 * (mean(n) - 1 - k) / ((1 - all_r2) * k)
########################################################################################################
