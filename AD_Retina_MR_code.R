#Script for analysing the causal risk factors of Alzheimer's disease identified from fundal imaging of the eye using bidirectional Mendelian Randomization (MR)
#For more information on commands/additional options see: https://mrcieu.github.io/TwoSampleMR/

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


setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Retina AD MR")


########### Estiamtes of the causal effect of AD on ocular features ####################
##########################################################################
##########################################################################
################## AD on RNFL thickness #########################

#### Kunkle  GWAS summary data extracted from the variants that showed best level of association after meta analysis of Stage 1 and Stage 2 Kunkle AD Gwas (from the list of published table 1)
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
head(ADx)


#### these columns are newly added to perform Steiger filtering for binary exposure (AD) 
ADx$exposure<-"AD"
ADx$prevalence.exposure<-0.1
ADx$ncase.exposure<-21982
ADx$ncontrol.exposure<-41944
ADx$units.exposure<-"log odds"


#### checking the pvalue: please note that 12 SNPs has pvalue<=5e-8 for the stage 1 kunkkle gwas
## however, they are found genome wide significant in the combined meta analysis of stage 1 and stage 2 ( link https://www.nature.com/articles/s41588-019-0358-2/tables/1)
length(which(ADx$pval.exposure<=5E-8)) # 12 snps

### no clumping performed as these variants were reported genome wide significant

write.csv(ADx,"./final kunkle outcome/New Kunkle data/R_AD_snps_for_MR.csv")

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

dim(RNFL)
head(RNFL)


RNFL$outcome<-"RNFL"
RNFL$units.outcome<-rep("SD", nrow(RNFL))
RNFL$samplesize.outcome<-rep(31434, nrow(RNFL))
head(RNFL)
length(which(RNFL$pvalue.outcome<=5E-8)) 

RNFL[, c("SNP","effect_allele.outcome", "other_allele.outcome", "beta.outcome")]
dat6 <- harmonise_data(
  exposure_dat = ADx, 
  outcome_dat = RNFL
)

dim(dat6)
dat6[, c("SNP","effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome", "beta.exposure", "beta.outcome")]

write.csv(dat6,"./final kunkle outcome/New Kunkle data/R_AD_RNFL_data.csv")
dat6<-read.csv("./final kunkle outcome/New Kunkle data/R_AD_RNFL_data.csv", header = TRUE, sep = ",")



##### steiger filtering:

dat6_steiger<-steiger_filtering(dat6)
table(dat6_steiger$steiger_dir)
dat6_RNFL_steiger<-subset(dat6_steiger, dat6_steiger$steiger_dir==TRUE)
dim(dat6_RNFL_steiger) 


#### MR analysis for entire dataset
mr_results<-mr(dat6, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

mr_results

## adjusting the beta by multiplying 0.693 for MR with binary exposure
mr_results$cbeta<-mr_results$b*0.693
mr_results$cse<-mr_results$se*0.693
mr_results$lcl<- mr_results$cbeta-1.96*mr_results$cse
mr_results$ucl<- mr_results$cbeta+1.96*mr_results$cse
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$cbeta,mr_results$lcl, mr_results$ucl,mr_results$pval)
results
names(results)<-c("Outcome", "Exposure", "N_SNPs", "MR_Methods", "Beta", "CI_Lower_Limit", "CI_Upper_Limit", "P_value")
results

write.csv(results,"./final kunkle outcome/New Kunkle data/R_AD_RNFL_results_all_snp.csv")





#### MR analysis for steiger filtered dataset
mr_results<-mr(dat6_RNFL_steiger, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

mr_results
## adjusting the beta by multiplying 0.693 for MR with binary exposure
mr_results$cbeta<-mr_results$b*0.693
mr_results$cse<-mr_results$se*0.693
mr_results$lcl<- mr_results$cbeta-1.96*mr_results$cse
mr_results$ucl<- mr_results$cbeta+1.96*mr_results$cse
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$cbeta,mr_results$lcl, mr_results$ucl,mr_results$pval)
results
names(results)<-c("Outcome", "Exposure", "N_SNPs", "MR_Methods", "Beta", "CI_Lower_Limit", "CI_Upper_Limit", "P_value")
results
# Export results
write.csv(results,"./final kunkle outcome/New Kunkle data/R_AD_RNFL_results_steiger.csv")


### MR analysis by excluding the snp rs429358

#dat6_exAPOE<-dat6_RNFL_steiger[ !dat6_RNFL_steiger$SNP == "rs429358", ]
#mr_results<-mr(dat6_exAPOE, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
#mr_results



# Run some sensitivity analyses
# c. Is there evidence of heterogeneity in the genetic effects?
het <- mr_heterogeneity(dat6)
het

# d. Is there evidence of horizontal pleiotropy?
pleio <- mr_pleiotropy_test(dat6)
pleio

pleio$egger_intercept-1.96*pleio$se
pleio$egger_intercept+1.96*pleio$se

res_single <- mr_singlesnp(dat6)
res_single



# Generate a scatter plot comparing the different methods

p1<-mr_scatter_plot(mr_results,dat6)
p1[[1]]
ggsave(p1[[1]], file="./final kunkle outcome/New Kunkle plot/AD_RNFL_main_scatter_plot.pdf", width=7, height=7)

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods

p2<-mr_forest_plot(res_single)
p2[[1]]
ggsave(p2[[1]], file="./final kunkle outcome/New Kunkle plot/AD_RNFL_main_forest_plot.pdf", width=7, height=7)

# Generate a funnel plot to check asymmetry

p3<-mr_funnel_plot(res_single)
p3[[1]]
ggsave(p3[[1]], file="./final kunkle outcome/New Kunkle plot/AD_RNFL_main_funnel_plot.pdf", width=7, height=7)

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat6)

p4<-mr_leaveoneout_plot(res_loo)
p4[[1]]
ggsave(p4[[1]], file="./final kunkle outcome/New Kunkle plot/AD_RNFL_main_leaveoneout_plot.pdf", width=7, height=7)

