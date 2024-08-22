---
title: "Combining Bidirectional Mendelian Randomization and Experimental Approaches for the Identification of miRNAs related to Major Depressive Disorder
"
author: "linlin Cao"
date: 2023.11.20 - 2024.8.23

---

# 1. Purpose 

There are over 280 million people suffer from Major depressive disorder(MDD) in 2023 around the world. 
A study conducted in Sweden revealed that around 15% of primary health care center attendees suffered from mental disorders, 
with depression being particularly prevalent. The clinical presentation of depression will lead to not only the poor quality life 
but also the burden to economy. Studies found the expression levels of miRNAs were associated with depression, given the possibility 
that miRNAs could served as biomarkers to detect depression. However, current researches mainly focused on miRNAs association with MDD, 
the casual relationship remains not fully understood. Mendelian Randomization (MR) method will help us to figure out the casual relationship 
between miRNAs and MDD. Furthermore, we applied the experiments to validate our identified miRNAs from ME analysis that have potential causal 
relationship with MDD based on blood samples from WHILA cohort.

MR method applied SNPs as instruments to investigate the casual relationship between exposure data(miRNAs) and outcome data(MDD).
For MR method explanation check article: <https://doi.org/10.1136/bmj.k601>

For MR related noun explanation: <https://mr-dictionary.mrcieu.ac.uk/>



# 2. Steps:

1. Download exposure and outcome data used in MR method
2. SNPs selection and harmonzing the exposure and outcome data
3. MR analysis
4. Experimental validation
5. MR results visualization (forest plot)
6. Power calculation for identified miRNAs (r^2 < 0.001)
7. Pathway analysis, result visualization (network)
8. Bi-directional MR analysis 
9. Experimental results analysis
10. RERI caculation for experimental results
11. phenome-wide scan analysis and visualization (phewas plot)

All analysis were conducted in R (version = 4.3.1). 

## Install and introduce pacakges we need
install.packages("remotes") * if you don't have this package

remotes::install_github("MRCIEU/TwoSampleMR") * install TwoSampleMR package from github 

* install_version("readxl", "1.4.3") 

```
if (!require("devtools")) { install.packages("devtools") } else {}

devtools::install_github("rondolab/MR-PRESSO")

library(remotes) # version 2.4.2.1

library(TwoSampleMR) # version 0.5.8

library(readxl) # version 1.4.3

library(dplyr) # Version 1.1.2

library(data.table) # version 1.14.8

library(MRPRESSO) # version 1.0
```



## 2.1 Download exposure and outcome data

Exposure data from Huan et al article (https://pubmed.ncbi.nlm.nih.gov/25791433/ ), can been downloaded from "Supplementary 4".

Exposure data contains 5,269 identified cis-eQTL associated with th expression level of miRNAs at FDR < 0.1.


Outcome data from Howar DM et al article (https://pubmed.ncbi.nlm.nih.gov/29662059/ ), can been downloaded from PGC website (https://pgc.unc.edu/for-researchers/download-results/). The filename is “Genome-wide Summary statistics from a meta-analysis of PGC and UK-Biobank (347.4Mb)”.

Outcome data contains 848,3301 SNPs assocaited Major depressive disorder, this is a GWAS summary result based on meta-analysis for PGC and UK Biobank cohort study from Wary et al for UK Biobank and Howard DM et al for PGC cohort. 

Download above exposure and outcome data and read into R project "binp52 project"



```
gwas_exposure_data <- readxl::read_excel("exposure_data.xlsx", skip = 1)
unique_miRNA <- unique(gwas_exposure_data$miRNA_FHS)
outcome_data <- fread("outcome_data.txt")
```

## 2.2 SNPs selection and harmonzing the exposure and outcome data
The careful selection of IVs (instrumental valiables which is SNPs in exposure data here) is essential before performing MR analysis, 

guided by three primary criteria: 

(1). IVs must exhibit a strong association (F-statistic > 10) with the exposure; 

(2). IVs should be independent of potential confounders; 

(3). IVs shouldn’t directly relate to the outcome without being linked to the exposure.

To meet with above assumptions, we performed filteration criteria like below : MAF > 0.01 & F > 10, 

phenoscanner to remove SNPs that may have strong potential relationship with confounders:
BMI, physical actiivities, drinking, and smoking
SNP has strong assocaition with above confounders could been found from phenoscanner website.
The dependent SNPs were removed using the clump_data () function in TwosampleMR package when r^2 > 0.001 withing the window size 10Mb, 
and the reference genotype data is Europeans 1000 Genomes project.


After read the exposure and outcome data using the TwosampleMR package, we performed the harmonzation to adjust the effect allele into the same direction 
for exposure and outcome data for later MR analysis step.


```
gwas_exposure_data <- gwas_exposure_data[gwas_exposure_data$MAF > 0.01 &
                      (gwas_exposure_data$Estimate^2) / (gwas_exposure_data$Std.Error^2) > 10, ]

# create a function get_clumped_data to get the harmonized data for each miRNA
get_allclumped_data <- function(miRNA, r_square) {
  test1 <- subset(gwas_exposure_data, miRNA_FHS == miRNA)
  test1 <- format_data(test1, type = "exposure", snp_col = "snpID", 
                       beta_col = "Estimate", se_col = "Std.Error", 
                       effect_allele_col = "effect", other_allele_col = "noneffect",
                       eaf_col = "obs_eaf", pval_col = "Pval")
  test_clumped <- clump_data(test1, clump_kb = 10000, clump_r2 = as.numeric(r_square))
  
  if (nrow(test_clumped) == 0) {
    return(data.frame(
      miRNA = miRNA, SNP = NA,
      effect_allele.exposure = NA, other_allele.exposure = NA,
      effect_allele.outcome = NA, other_allele.outcome = NA,
      beta.exposure = NA, beta.outcome = NA,
      eaf.exposure = NA, eaf.outcome = NA,
      remove = NA, palindromic = NA, ambiguous = NA,
      id.outcome = NA, se.outcome = NA,
      pval.outcome = NA, outcome = NA,
      mr_keep.outcome = NA, pval_origin.outcome = NA,
      se.exposure = NA, pval.exposure = NA, exposure = NA, mr_keep.exposure = NA, pval_origin.exposure = NA,
      id.exposure = NA, action = NA, SNP_index = NA, mr_keep = NA, samplesize.outcome = NA
    ))
  }
  
  test2 <- merge(test_clumped, outcome_data, by.x = "SNP", by.y = "MarkerName")
  if (nrow(test2) == 0){
    return(data.frame(
      miRNA = miRNA, SNP = NA,
      effect_allele.exposure = NA, other_allele.exposure = NA,
      effect_allele.outcome = NA, other_allele.outcome = NA,
      beta.exposure = NA, beta.outcome = NA,
      eaf.exposure = NA, eaf.outcome = NA,
      remove = NA, palindromic = NA, ambiguous = NA,
      id.outcome = NA, se.outcome = NA,
      pval.outcome = NA, outcome = NA,
      mr_keep.outcome = NA, pval_origin.outcome = NA,
      se.exposure = NA, pval.exposure = NA, exposure = NA, mr_keep.exposure = NA, pval_origin.exposure = NA,
      id.exposure = NA, action = NA, SNP_index = NA, mr_keep = NA, samplesize.outcome = NA
    ))
  }



# get the harmonized results of LD pruning criteria ( r^2 = 0.001) 

miRNA_results_list_0.001 <- lapply(unique(gwas_exposure_data$miRNA_FHS), 
                             function(miRNA) get_allclumped_data(miRNA, 0.001))
harmonized_df_0.001 <- do.call(rbind, miRNA_results_list_0.001)
write.table(harmonized_df, "harmonized_data_0.001.csv", sep = ",", row.names = F)



# remove SNPs that may have strong potential relationship with confounders :
# BMI, phydical actiivities, drinking and smoking
# SNP has strong association with above confounders could been found from phenoscanner website (P < 5*10(-8)).

# No selected SNP shown the strong association with above confounders on the phenoscanner website.


```

This step will generate harmonized files containing MR results based on LD pruning criteria,

"harmonized_data_0.001.csv"

Theharmonized data file was used in later mr-presso and mr analysis.



## 2.3 MR analysis 
### including MRPRESSO
MRPRESSO global test was performed to identify the possible bias from the horizontal pleiotrophy.

MRPRESSO can only been performed for those miRNAs which have over 3 SNPs.

We removed SNP which MRPRESSO p < 0.05 (considered as has horizontal pleitrophy) for each miRNA

Since our six identified miRNAs have only one corresponding IV remained, we cannot conduct MRPRESSO for them, so we didn't show it in the material and method part

After SNP selection, MR analysis will apply Wald ratio method for miRNA has only one SNP, and IVW method for miRNA has over one SNP.


 
```
# This function will perform MRPRESSO global test and MR analysis miRNA.

add_mrpresso <- function(harmonized_data){
  all_mr_results <- data.frame()
  for (RNA in unique(harmonized_data$miRNA)) {
    harmonized_subset <- harmonized_data[harmonized_data$miRNA == RNA, ]
    harmonized_subset <- harmonized_subset[, -1]
    # Perform MR analysis
    if (nrow(harmonized_subset) > 3) {
      # mr presso global test to remove outliers for over 3 IVs miRNA, 
      # 3 IVs is based on test, since miR_668 has 3 snps but cannot pass the mr-presso package process instrumental variables requirement number
      mr_presso_result <- suppressWarnings({
        mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                  data = harmonized_subset, NbDistribution = 1000,  
                  SignifThreshold = 0.05)
      })
      
      # Check for outliers and perform MR analysis
      if (mr_presso_result[["MR-PRESSO results"]]$`Global Test`$Pvalue < 0.05) {
        mr_presso_filtresu <- subset(harmonized_subset, !mr_presso_result$global_outliers)
        
        mr_test_resu <- mr(mr_presso_filtresu, method_list = "mr_ivw")
        mr_test_resu <- cbind(miRNA = rep(RNA, nrow(mr_test_resu)), mr_test_resu)
        all_mr_results <- rbind(all_mr_results, mr_test_resu)} else {
          
          mr_test_resu <- mr(harmonized_subset, method_list = "mr_ivw")
          mr_test_resu <- cbind(miRNA = rep(RNA, nrow(mr_test_resu)), mr_test_resu)
          all_mr_results <- rbind(all_mr_results, mr_test_resu)}
      
    }
    else {
      mr_test_resu <- mr(harmonized_subset, method_list = c("mr_wald_ratio", "mr_ivw"))
      mr_test_resu <- cbind(miRNA = rep(RNA, nrow(mr_test_resu)), mr_test_resu)
      all_mr_results <- rbind(all_mr_results, mr_test_resu)}
  }
  return (all_mr_results)
}
all_mr_results_0.001 <- add_mrpresso(harmonized_df_0.001)

# we performed the BH multiple test corresction to avoid the fasle positive error.
adjusted_mr_result_0.001 <- all_mr_results_0.001 %>%
  group_by( method) %>%
  mutate(adj.P = p.adjust(pval, method = "BH"))


OR_CI_result_0.001 <- generate_odds_ratios(mr_res = adjusted_mr_result_0.001)

OR_CI_result_0.001 <- subset(OR_CI_result_0.001, method %in% c("Wald ratio", "Inverse variance weighted"))

write.csv(OR_CI_result_0.001, "OR_CI_result_0.001.csv", row.names = FALSE)

```
This step generated three MR results files based on distinct LD pruning criteria (r^2 = 0.001)

"OR_CI_result_0.001.csv"



The significant identified miRNAs which related with MDD were retrieved based on P < 0.05. 

The results after BH corresction was not presented in our result, since we want keep as amny as candidate miRNAs to further experimental validation.




## 2.4 Experimental validation

We got 6 miRNAs are significantly associated with MDD. Based on previous study, We removed 2 miRNAs which CT value > 38, 
indicating these miRNAs expression were undetectable in our WHILA plasma samples. 

power caculation decides the number we needed to do experimental validation avoiding the false positive and negative situdation (power = 80%, significance 0.05)
```
# Set parameters
effect_size <- 1.49  # miR-132 expression level mean difference (case - control) based on previous study, divided by pooled sd
alpha <- 0.05     
power <- 0.8       


# perform the pwr test 
power.t.test(delta = 1.4, sd = 1.49, sig.level = 0.05, power = 0.8, type = "two")
pwr.t.test(d= 1.4/1.49, sig.level = 0.05, power = 0.8, type = "two.sample", alternative = "two.sided")

# two-side 62 individuals, one group 30 people.
```
Then we performed qRT-PCR to do validation for above 4 miRNAs based on 104 samples from WHILA population. 



## 2.5 MR results visualization (forest plot) 

Three Forest plots were used to present  MR results for identified miRNAs, and interval was 95%CIs for each miRNA.


```
# forest figure for identified miRNA 0.001

OR_CI_result_0.001_identified <- OR_CI_result_0.001[OR_CI_result_0.001$pval < 0.05 , ]
forest_plot_0.001 <- ggplot(OR_CI_result_0.001_identified, 
                            mapping = aes(x = or, y = reorder(miRNA, or_uci95))) +
                            geom_linerange(aes(xmin = or_lci95, xmax = or_uci95)) +
                            geom_point(aes(x = or, y = reorder(miRNA, -or_uci95), color = pval))+
                            geom_vline(xintercept = 1.00, linetype = "dashed") +
                            scale_color_gradient(low = "black", high = "gray") +
                            xlab("OR") +
                            ylab("Exposure") +
                            theme_minimal(base_size = 15) 
forest_plot_0.001

ggsave("forest_plot_0.001.png", plot = forest_plot_0.001, width = 10, height = 15)

```
The forest plot is in "forest_plot_0.001.png" file.


## 2.6  Power calculation for identified miRNAs (r^2 < 0.001)
Due to the substantial difference in sample size between the exposure and outcome datasets (exposure: 5,239 individuals; outcome: 500,199 individuals), 
it is necessary to assess the power of our study. Power calculation helps to assess whether our study can detect causal effects for the exposure on the 
outcome based on current sample size, effect size of IVs on the exposure, and significance level. When the power is low, it may cause a false negative result. 
We assessed the study power of IVs for identified miRNAs. 

First, we calculated the proportion of variance explained for the association between the IV and the miRNA using the formula R2=2×β2×MAF×(1-MAF),
where β represents the effect size of the IV on the exposure and MAF represents the minor allele frequency.

The study power was further examined on the website (https://shiny.cnsgenomics.com/mRnd/) using the samples size of exposure, outcome, R2, 
our causal effect size of the exposure-outcome, and significance level typically used 0.05.



## 2.7  Pathway analysis, result visualization (network)


To get the biological significance for the miRNAs and understand their target genes and invovled pathways, we performed the in-silicon search based on miRsystem, miRwalk, and miRnet website for identified six miRNAs (except the miR_635 which is not found on website).
The in-silicon results are shown in the file pathway analysis -> 
We combined the target genes come from three platforms above, and performed the pathway enrichemnt analysis (over-representative analysis) based on KEGG and GO databases (biological process). Significantly enriched KEGG pathways (P < 0.05) for all the combined target genes for each identified miRNA were selected, the significant pathways for all identified miRNAs are presented in Supplmentary files "Table S4a". Specifically, pathways containing terms like “pathway” or “Pathways” in their names were selected, and genes targeted by at least four miRNAs were included. The remaining miRNAs, target genes, and involved pathways are presented in Supplmentary files "Table S4b". A network illustrating the interactions between miRNAs, their filtered target genes and pathways was constructed using Cytoscape (version = 3.10.2). After filtering target genes from three websites for the KEGG results for five miRNAs (excluding miR-635-3p, which was not recorded), we identified 30 genes, 44 pathways, and five identified miRNAs involved in the interaction network. 


The result for the KEGG pathway analysis were shown in cytoscape to show all the miRNAs and their target genes and enriched pathways.

The script for pathway analysis and in-silicon results are shown in pathway analysis filefolder above. 

The three_web_all_genes.R will use the miRwalk, miRsystem, and miRnet in-silicon results which shown the target genes for identified 5 miRNAs in each corresponding plateform, 
the files "edge_pathname_interbig3_file.xlsx" and "type_pathname_interbig3_file.xlsx" are used to generate network using cytoscape for kegg result, and the
file "combined_pathname_all_pathways_kegg_p0.05_interaction3.xlsx" which is the tableS4a described above.

And also the file "combined_common_validated_gene_GO_top50_neuro_related_reorder.xlsx" (TableS4c) will show the results for GO analysis containing top50 significantly enriched and nurorelated pathways for each miRNAs validated target genes been predicted in three plateforms above, the important neuro related pathways that may related to MDD are selected from "combind_GO_top50pathway_5miRNA_allgenes_3web_genesymbol.xlsx" files, their  miRNAs, target genes, and related pathways are shown in "neuro_related_GO_pathway_fortop50combined.xlsx" file we uploaded.







## 2.8  Bi-directional MR analysis 
```
outcome_data <- fread("outcome_data.txt")
test <- select(gwas_exposure_data, c("snpID","miRNA_FHS","effect", "noneffect",
                                     "obs_eaf", "Estimate", "Std.Error", "Pval" ))
outcome_data$LogOR <- outcome_data$LogOR / (log(exp(1)))

# reverse analysis
bi_exposure_data <- outcome_data
bi_exposure_data$F <- (bi_exposure_data$LogOR^2) / (bi_exposure_data$StdErrLogOR^2)
bi_exposure_data <- bi_exposure_data[bi_exposure_data$Freq > 0.01 & bi_exposure_data$F > 10, ]
bi_exposure_data <- subset(bi_exposure_data, P < 5e-5)

colnames(bi_exposure_data) <- c("snpID", "effect_allele", "other_allele",
                                "eaf", "beta", "se", "pval","F")
bi_exposure_data <- data.frame(bi_exposure_data)
bi_exposure_data <- format_data(bi_exposure_data, snp_col = "snpID",
            effect_allele_col = "effect_allele", other_allele_col = "other_allele",
            beta_col = "beta", se_col = "se", 
            pval_col = "pval", eaf_col = "eaf")
dim(bi_exposure_data) # remians 4625 snps for p < 5*10(-8), 25741 snps for p < 5*10(-5)
# clump data
bi_clump_0.001 <- clump_data(bi_exposure_data, clump_kb = 10000, clump_r2 = 0.001) # 50 snps remains for p < 5*10(-8), 332 snps remians for p < 5*10(-5)
bi_clump_0.01 <- clump_data(bi_exposure_data, clump_kb = 10000, clump_r2 = 0.01) # 53 snps remains for p < 5*10(-8), 577 snps remians for p < 5*10(-5)
bi_clump_0.1 <- clump_data(bi_exposure_data, clump_kb = 10000, clump_r2 = 0.1) # 66 snps remains for p < 5*10(-8), 782 snps remians for p < 5*10(-5)

# select same snp in gwas_exposure data
bi_clump_0.001.test <- merge(bi_clump_0.001, gwas_exposure_data, by.x = "SNP", by.y = "snpID")
bi_clump_0.01.test <- merge(bi_clump_0.01, gwas_exposure_data, by.x = "SNP", by.y = "snpID")
bi_clump_0.1.test <- merge(bi_clump_0.1, gwas_exposure_data, by.x = "SNP", by.y = "snpID")
dim(bi_clump_0.1.test) # there is no snps in bi_outcome data for r^2 0.1,0.01, 0.001 for p < 5*10(-8)
                        # there is no snps in bi_outcome data (gwas_exposure_data) for r^2 0.001 for P < 5*10(-5)
                        # there is 1 snp in bi_outcome data (gwas_exposure_data) for r^2 0.01 for P < 5*10(-5)
# there is 1 snp in bi_outcome data (gwas_exposure_data) for r^2 0.1 for P < 5*10(-5)
```
Above results suggesting there is no corresponding SNP remaining after performing the same criteria as the forward MR analysis (LD pruning r^2 < 0.001)
Thus, no evidence showed there is no reverse causation for miRNA and MDD.

## 2.9 Experimental results analysis 

We reported the clinical chracteristics of WHILA population, baseline age and body mass index (BMI) were reported as mean and standard deviation (SD), while smoking and alcohol consumption were expressed as percentages. An independent sample t-test was employed to compare age and BMI between cases and controls, whereas smoking and alcohol use were assessed using a chi-square test. To investigate the association between cases and controls for miRNA expression levels at baseline, logistic regression was conducted, with adjustments made for age. Statistical significance was considered at p < 0.05. 


The file "mirna_mental.csv" contains clinical chracteristics of WHILA population for 104 sampes we used, if you need to get these data, pls ask the permission from CPF. 

The file "new_experiment_results.xlsx" contains our experimental results for 5 identified miRNAs from MR analysis results (LD pruning r^2 = 0.001)
"mirna_mental.csv" and "new_experiment_results.xlsx" are accessible after requesting for permission.

```
script:  experimental_results_analysis.R

```

All statistical results can been found in code with the beginning of "#check"

## 2.10 RERI calculation
Please check RERI calculation fileforlder for script and experimental results.
The new_mirna2.xlsx file is avaliable after requesting for permisson.

```
script:  RERI_calculation.R

```

The result are shown after the code "#check" in each code line.

## 2.11 phenome-wide scan analysis and visualization (phewas plot)
phenome-wide scan analysis file folder contians phenome-wide scan analysis.R script and phewas results for identified miRNAs corresponding IV in GWASATLAS.
'''
script phenome-wide scan analysis.R
'''
The script will generate a file forlder names phewas plot revised y axis and this file contains figure showing the phewas plots names rs685149_version1.jpeg, rs3803809_version1.jpeg, rs12479469_version1.jpeg.
