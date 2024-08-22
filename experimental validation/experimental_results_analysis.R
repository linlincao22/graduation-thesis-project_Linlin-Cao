library(dplyr)
library(readxl)
library(ggplot2)
mirna_mental.csv <- read.table("mirna_mental.csv", sep = ";", header = T)
mirna_mental.csv$bmi[mirna_mental.csv$bmi ==  " " ] <- "NA"

View(mirna_mental.csv)


# randomly selected samples after power calculation
# select 10 samples which got MDD withing three years after collection
case_10 <- subset(mirna_mental.csv, mirna_mental.csv$case_control == "Case" & mirna_mental.csv$time_MDD_years < 4  )
View(case_10)

case_10.test<- case_10[case_10$pkod %in% sample(case_10$pkod, 10),]
View(case_10.test)

control_10 <- subset(mirna_mental.csv, mirna_mental.csv$case_control == "Control"  )
control_10.test<- control_10[control_10$pkod %in% sample(control_10$pkod, 10),]
View(control_10.test)


### table 1 

##  calculate mean and SD for age and BMI, 52 case and 52 control and 104 total
mean_age_total <- mean(mirna_mental.csv$age_baseline) # 56.00962
mean_age_case <- mean(mirna_mental.csv$age_baseline[mirna_mental.csv$case_control == "Case"]) # 55.71154
mean_age_control <- mean(mirna_mental.csv$age_baseline[mirna_mental.csv$case_control == "Control"]) # 56.30769

sd_age_total <- sd(mirna_mental.csv$age_baseline) # 2.781685
sd_age_case <- sd(mirna_mental.csv$age_baseline[mirna_mental.csv$case_control == "Case"]) # 2.569325
sd_age_control <- sd(mirna_mental.csv$age_baseline[mirna_mental.csv$case_control == "Control"]) # 2,973995



# change blank cell in BMI column into NA and change "," in the middle of each value into "."to calculate mean and sd
mean_BMI_total <- mean(as.numeric(lapply(mirna_mental.csv$bmi[mirna_mental.csv$bmi != "NA"],
                                         function(x) gsub(",", ".", x)))) # 24.50875
mean_BMI_case <- mean(as.numeric(lapply(mirna_mental.csv$bmi[mirna_mental.csv$bmi != "NA" & mirna_mental.csv$case_control == "Case"],
                                        function(x) gsub(",", ".", x)))) # 24.96427
mean_BMI_control <- mean(as.numeric(lapply(mirna_mental.csv$bmi[mirna_mental.csv$bmi != "NA" & mirna_mental.csv$case_control == "Control"],
                                            function(x) gsub(",", ".", x)))) # 24.04333

sd_BMI_total <- sd(as.numeric(lapply(mirna_mental.csv$bmi[mirna_mental.csv$bmi != "NA"],
                                    function(x) gsub(",", ".", x)))) # 3.924952
sd_BMI_case <- sd(as.numeric(lapply(mirna_mental.csv$bmi[mirna_mental.csv$bmi != "NA" & mirna_mental.csv$case_control == "Case"],
                                       function(x) gsub(",", ".", x)))) # 3.856467
sd_BMI_control <- sd(as.numeric(lapply(mirna_mental.csv$bmi[mirna_mental.csv$bmi != "NA" & mirna_mental.csv$case_control == "Control"],
                                       function(x) gsub(",", ".", x)))) # 3.981731
length(mirna_mental.csv$bmi[mirna_mental.csv$bmi != "NA" & mirna_mental.csv$case_control == "Case"])
# check 
mean_age_total # 56.00962
mean_age_case #  55.71154
mean_age_control # 56.30769

sd_age_total # 2.781685
sd_age_case # 2.569325
sd_age_control # 2.973995



mean_BMI_total # 24.50875
mean_BMI_case # 24.96427
mean_BMI_control # 24.04333

sd_BMI_total # 3.924952
sd_BMI_case # 3.856467
sd_BMI_control # 3.981731


## current or past smoking (value = 0,1) and severe-alcohol (value = 2) percentage for total, case and control
observed_total_number_sm <- sum(!is.na(mirna_mental.csv$smoke_cat)) # 101 
smoking_number <- sum(mirna_mental.csv$smoke_cat %in% c(0,1)) # 26
smoking_percen <- smoking_number / observed_total_number_sm 
smoking_case_percen <- sum(mirna_mental.csv$smoke_cat %in% c(0,1) & mirna_mental.csv$case_control == "Case" & !is.na(mirna_mental.csv$smoke_cat)) /
                           sum(!is.na(mirna_mental.csv$smoke_cat) & mirna_mental.csv$case_control == "Case")
smoking_control_percen <- sum(mirna_mental.csv$smoke_cat %in% c(0,1) & mirna_mental.csv$case_control == "Control" & !is.na(mirna_mental.csv$smoke_cat)) /
                              sum(!is.na(mirna_mental.csv$smoke_cat) & mirna_mental.csv$case_control == "Control")

sev_alc_number <- sum(mirna_mental.csv$alc_cat %in% 2)
sev_alc_percen <- sev_alc_number / 104
sev_alc_case_percen <- sum(mirna_mental.csv$alc_cat %in% 2 & mirna_mental.csv$case_control == "Case") / 52
sev_alc_control_percen <- sum(mirna_mental.csv$alc_cat %in% 2 & mirna_mental.csv$case_control == "Control") / 52

# check
observed_total_number_sm # 101
smoking_number # 26
smoking_percen # 25.74257%
smoking_case_percen # 25.49021%
smoking_control_percen # 26%

# observed alcohol people 104
sev_alc_number # 20
sev_alc_percen # 19.2%
sev_alc_case_percen # 17.30769%
sev_alc_control_percen # 21.15385




## parameter test and non-parameter test for age, bmi, smoking, alcohol
t.test_age <- t.test(mirna_mental.csv$age_baseline[mirna_mental.csv$case_control == "Case"], 
       mirna_mental.csv$age_baseline[mirna_mental.csv$case_control == "Control"]) # p = 0.2767
t.test_bmi <- t.test(as.numeric(lapply(mirna_mental.csv$bmi[mirna_mental.csv$bmi != "NA" & mirna_mental.csv$case_control == "Case"],
                                       function(x) gsub(",", ".", x))),
                     as.numeric(lapply(mirna_mental.csv$bmi[mirna_mental.csv$bmi != "NA" & mirna_mental.csv$case_control == "Control"],
                                       function(x) gsub(",", ".", x))))

# check 
t.test_age # p = 0.2767
t.test_bmi # p = 0.2603


## non-parameter test for smoking (0,1) or not (2)
mirna_mental.csv$smoke_cat[mirna_mental.csv$smoke_cat == 0 ] <- "smokers"
mirna_mental.csv$smoke_cat[mirna_mental.csv$smoke_cat == 1 ] <- "smokers"
mirna_mental.csv$smoke_cat[mirna_mental.csv$smoke_cat == 2] <- "non smokers"

data <- data.frame(group = mirna_mental.csv$case_control,
                   smoking_status =  mirna_mental.csv$smoke_cat)

data <- na.omit(data)

test <- table(data$group, data$smoking_status)

smoking_chisq_res <- chisq.test(test)


# non-parameter test for alcohol (2) and alcohol (0,1) in case and control
mirna_mental.csv$alc_cat[mirna_mental.csv$alc_cat == 0] <- "non alcohol"
mirna_mental.csv$alc_cat[mirna_mental.csv$alc_cat == 1] <- "non alcohol"
mirna_mental.csv$alc_cat[mirna_mental.csv$alc_cat == 2] <- "sev_alcohol"

data_alc <- data.frame(group = mirna_mental.csv$case_control,
                       alc_status =  mirna_mental.csv$alc_cat)

data_alc <- na.omit(data_alc)

test_alc <- table(data_alc$group, data_alc$alc_status)
test_alc
alc_chisq_res <- chisq.test(test_alc) # p = 0.8035

# check 
alc_chisq_res # p = 0.8035
smoking_chisq_res # p = 1


## introducing our experimental results and remove old results from others
new_mirna <- mirna_mental.csv[, !names(mirna_mental.csv) %in% c("dCt_miR_17_5p", "dCt_miR_134_5p" , "dCt_miR_144_5p" ,"dCt_let_7b_5p",  "dCt_let_7c_5p")]
experi_results <- read_excel("new_experiment results.xlsx")
new_mirna_test1 <- merge(experi_results, new_mirna, by = "pkod")
new_mirna <- new_mirna_test1[, !names(new_mirna_test1) %in% c("Cq Mean_425","Cq Mean_133",  "Cq Mean_132" , "Cq Mean_598", "Cq Mean_130", "Cq Mean_629" )]
new_mirna$`delta 598`<- as.numeric(new_mirna$`delta 598`)
new_mirna$case_control <- as.factor(ifelse(new_mirna$case_control == "Case", 1, 0))

## table 2: test delta mirna expression level statistical difference between case and control
## t test for normal distribution, wilxcon test for skewed distribution

# have a look at the distribution of each delta mirna
hist(new_mirna$`delta 132`, main = "mirna_132 distribution histogram", xlab = "value", ylab = "count")
hist(new_mirna$`delta 133`, main = "mirna_133 distribution histogram", xlab = "value", ylab = "count")
hist(as.numeric(new_mirna$`delta 598`), main = "mirna_598 distribution histogram", xlab = "value", ylab = "count")
hist(new_mirna$`delta 130`, main = "mirna_130 distribution histogram", xlab = "value", ylab = "count")
hist(new_mirna$`delta 629`, main = "mirna_629 distribution histogram", xlab = "value", ylab = "count")

# check
# all five mirna data are in normal distribution

# mean and median and sd of each mirna
mean_132_case <- mean(new_mirna$`delta 132`[!is.na(new_mirna$`delta 132`) & new_mirna$case_control == 1 ])
mean_132_control <- mean(new_mirna$`delta 132`[!is.na(new_mirna$`delta 132`) & new_mirna$case_control == 0 ])
median_132_case <- median(new_mirna$`delta 132`[!is.na(new_mirna$`delta 132`) & new_mirna$case_control == 1 ])
median_132_control <- median(new_mirna$`delta 132`[!is.na(new_mirna$`delta 132`) & new_mirna$case_control == 0 ])
sd_132_case <- sd(new_mirna$`delta 132`[!is.na(new_mirna$`delta 132`) & new_mirna$case_control == 1 ])
sd_132_control <- sd(new_mirna$`delta 132`[!is.na(new_mirna$`delta 132`) & new_mirna$case_control == 0 ])
t.test_132 <-  t.test(new_mirna$`delta 132`[!is.na(new_mirna$`delta 132`) & new_mirna$case_control == 1 ],
                      new_mirna$`delta 132`[!is.na(new_mirna$`delta 132`) & new_mirna$case_control == 0 ])

mean_133_case <- mean(new_mirna$`delta 133`[!is.na(new_mirna$`delta 133`) & new_mirna$case_control == 1 ])
mean_133_control <- mean(new_mirna$`delta 133`[!is.na(new_mirna$`delta 133`) & new_mirna$case_control == 0 ])
median_133_case <- median(new_mirna$`delta 133`[!is.na(new_mirna$`delta 133`) & new_mirna$case_control == 1 ])
median_133_control <- median(new_mirna$`delta 133`[!is.na(new_mirna$`delta 133`) & new_mirna$case_control == 0 ])
sd_133_case <- sd(new_mirna$`delta 133`[!is.na(new_mirna$`delta 133`) & new_mirna$case_control == 1 ])
sd_133_control <- sd(new_mirna$`delta 133`[!is.na(new_mirna$`delta 133`) & new_mirna$case_control == 0 ])
t.test_133 <-  t.test(new_mirna$`delta 133`[!is.na(new_mirna$`delta 133`) & new_mirna$case_control == 1 ],
                      new_mirna$`delta 133`[!is.na(new_mirna$`delta 133`) & new_mirna$case_control == 0 ])


mean_130_case <- mean(new_mirna$`delta 130`[!is.na(new_mirna$`delta 130`) & new_mirna$case_control == 1 ])
mean_130_control <- mean(new_mirna$`delta 130`[!is.na(new_mirna$`delta 130`) & new_mirna$case_control == 0 ])
median_130_case <- median(new_mirna$`delta 130`[!is.na(new_mirna$`delta 130`) & new_mirna$case_control == 1 ])
median_130_control <- median(new_mirna$`delta 130`[!is.na(new_mirna$`delta 130`) & new_mirna$case_control == 0 ])
sd_130_case <- sd(new_mirna$`delta 130`[!is.na(new_mirna$`delta 130`) & new_mirna$case_control == 1 ])
sd_130_control <- sd(new_mirna$`delta 130`[!is.na(new_mirna$`delta 130`) & new_mirna$case_control == 0 ])
t.test_130 <-  t.test(new_mirna$`delta 130`[!is.na(new_mirna$`delta 130`) & new_mirna$case_control == 1 ],
                      new_mirna$`delta 130`[!is.na(new_mirna$`delta 130`) & new_mirna$case_control == 0 ])


mean_598_case <- mean(new_mirna$`delta 598`[!is.na(new_mirna$`delta 598`) & new_mirna$case_control == 1 ])
mean_598_control <- mean(new_mirna$`delta 598`[!is.na(new_mirna$`delta 598`) & new_mirna$case_control == 0 ])
median_598_case <- median(new_mirna$`delta 598`[!is.na(new_mirna$`delta 598`) & new_mirna$case_control == 1 ])
median_598_control <- median(new_mirna$`delta 598`[!is.na(new_mirna$`delta 598`) & new_mirna$case_control == 0 ])
sd_598_case <- sd(new_mirna$`delta 598`[!is.na(new_mirna$`delta 598`) & new_mirna$case_control == 1 ])
sd_598_control <- sd(new_mirna$`delta 598`[!is.na(new_mirna$`delta 598`) & new_mirna$case_control == 0 ])
t.test_598 <-  t.test(new_mirna$`delta 598`[!is.na(new_mirna$`delta 598`) & new_mirna$case_control == 1 ],
                      new_mirna$`delta 598`[!is.na(new_mirna$`delta 598`) & new_mirna$case_control == 0 ])


mean_629_case <- mean(new_mirna$`delta 629`[!is.na(new_mirna$`delta 629`) & new_mirna$case_control == 1 ])
mean_629_control <- mean(new_mirna$`delta 629`[!is.na(new_mirna$`delta 629`) & new_mirna$case_control == 0 ])
median_629_case <- median(new_mirna$`delta 629`[!is.na(new_mirna$`delta 629`) & new_mirna$case_control == 1 ])
median_629_control <- median(new_mirna$`delta 629`[!is.na(new_mirna$`delta 629`) & new_mirna$case_control == 0 ])
sd_629_case <- sd(new_mirna$`delta 629`[!is.na(new_mirna$`delta 629`) & new_mirna$case_control == 1 ])
sd_629_control <- sd(new_mirna$`delta 629`[!is.na(new_mirna$`delta 629`) & new_mirna$case_control == 0 ])
t.test_629 <-  t.test(new_mirna$`delta 629`[!is.na(new_mirna$`delta 629`) & new_mirna$case_control == 1 ],
                      new_mirna$`delta 629`[!is.na(new_mirna$`delta 629`) & new_mirna$case_control == 0 ])


# check 
new_test_idx <- order(new_mirna$case_control)
new_test_idx
new_ordered_mirna <- new_mirna[new_test_idx,]
View(new_ordered_mirna)


mean_132_case # -4.187755
mean_132_control # -3.819796
median_132_case # -4.21
median_132_control # -3.92
sd_132_case # 0.6774285
sd_132_control # 0.92921
t.test_132 # p value = 0.02762

mean_133_case # -2.8544
mean_133_control # -3.188958
median_133_case # -2.795
median_133_control # -3.32
sd_133_case # 1.137362
sd_133_control # 1.297856
t.test_133 # p value = 0.1787

mean_130_case # 0.01705882
mean_130_control # -0.256
median_130_case # -0.15
median_130_control # -0.29
sd_130_case # 0.659722
sd_130_control # 0.6240716
t.test_130 # p value = 0.03505

mean_598_case # -3.916327
mean_598_control # -4.033404
median_598_case # -4.05
median_598_control # -4.27
sd_598_case # 0.6529858
sd_598_control # 1.12238
t.test_598 # p value = 0.5363


mean_629_case # -3.945294
mean_629_control # -3.801429
median_629_case # -4.03
median_629_control # -4.07
sd_629_case # 1.164325
sd_629_control # 1.045498
t.test_629 # p value = 0.5168



## logistic regression

# logistic regression without age, bmi adjusted
log_132 <- glm(case_control ~ `delta 132`, 
               data = new_mirna, family = binomial, na.action = na.exclude)

log_133 <- glm(case_control ~ `delta 133`, 
               data = new_mirna, family = binomial, na.action = na.exclude)

log_130 <- glm(case_control ~ `delta 130`, 
               data = new_mirna, family = binomial, na.action = na.exclude)

log_598 <- glm(case_control ~ `delta 598`, 
               data = new_mirna, family = binomial, na.action = na.exclude)
 
log_629 <- glm(case_control ~ `delta 629`, 
               data = new_mirna, family = binomial, na.action = na.exclude)


# 95 CI interval, interval dont including 0, statistical significant.
confic_95ci_132 <- exp(confint(log_132)[2,])
confic_95ci_133 <- exp(confint(log_133)[2,])
confic_95ci_130 <- exp(confint(log_130)[2,])
confic_95ci_598 <- exp(confint(log_598)[2,])
confic_95ci_629 <- exp(confint(log_629)[2,])

# log(OR) = coef(log_132)[2]
OR_132 <- exp(coef(log_132)[2])
OR_133 <- exp(coef(log_133)[2])
OR_130 <- exp(coef(log_130)[2])
OR_598 <- exp(coef(log_598)[2])
OR_629 <- exp(coef(log_629)[2])



# check 
summary(log_132) # p = 0.0327*
OR_132 # 0.5597074 , means 1 unit of miRNA 132 increasing,  the possibility of getting case decreased 0.559
confic_95ci_132 #  0.3176055 0.9309509 

summary(log_133) # P = 0.181
OR_133 # 1.263692
confic_95ci_133 # 0.9064936 1.8120590

summary(log_130) # p = 0.0419 *
OR_130 #  2.011052
confic_95ci_130 # 1.061790 4.153805

summary(log_598) # p = 0.528
OR_598 #  1.154412
confic_95ci_598 # 0.7407306 1.8264510

summary(log_629) # p = 0.515
OR_629 #  0.8869728 
confic_95ci_629 # 0.6085727 1.2702380

# logistic regression with age adjusted
# adjust bmi value , turn into numeric and change the "," into "." 
new_mirna$bmi <-  as.numeric(lapply(new_mirna$bmi, function(x) gsub(",", ".", x)))

log_132_adj_age <- glm(case_control ~ `delta 132` + age_baseline   , 
                        data = new_mirna, family = binomial, na.action = na.exclude)

log_133_adj_age <- glm(case_control ~ `delta 133` + age_baseline , 
                        data = new_mirna, family = binomial, na.action = na.exclude)

log_130_adj_age  <- glm(case_control ~ `delta 130` + age_baseline , 
                         data = new_mirna, family = binomial, na.action = na.exclude)

log_598_adj_age  <- glm(case_control ~ `delta 598` + age_baseline  , 
                         data = new_mirna, family = binomial, na.action = na.exclude)

log_629_adj_age  <- glm(case_control ~ `delta 629` + age_baseline  , 
                         data = new_mirna, family = binomial, na.action = na.exclude)

# 95 CI interval, interval dont including 0, statistical significant.
confic_95ci_132_adj_age <- exp(confint(log_132_adj_age)[2,])
confic_95ci_133_adj_age <- exp(confint(log_133_adj_age)[2,])
confic_95ci_130_adj_age <- exp(confint(log_130_adj_age)[2,])
confic_95ci_598_adj_age <- exp(confint(log_598_adj_age)[2,])
confic_95ci_629_adj_age <- exp(confint(log_629_adj_age)[2,])

# log(OR) = coef(log_132)[2]
OR_132_adj_age <- exp(coef(log_132_adj_age)[2])
OR_133_adj_age <- exp(coef(log_133_adj_age)[2])
OR_130_adj_age <- exp(coef(log_130_adj_age)[2])
OR_598_adj_age <- exp(coef(log_598_adj_age)[2])
OR_629_adj_age <- exp(coef(log_629_adj_age)[2])


# check
summary(log_132_adj_age) # p = 0.0216 *
summary(log_133_adj_age) # p = 0.194
summary(log_130_adj_age) # p = 0.037 *
summary(log_598_adj_age) # p = 0.605
summary(log_629_adj_age) # p = 0.585

confic_95ci_132_adj_age #  0.2888053 0.8835155 
confic_95ci_133_adj_age #  0.8998505 1.8007906 
confic_95ci_130_adj_age #  1.082140 4.281113 
confic_95ci_598_adj_age #  0.7182494 1.7875313 
confic_95ci_629_adj_age #  0.6194637 1.2989638 

OR_132_adj #  0.506546 
OR_133_adj #  1.19547 
OR_130_adj #  1.843425 
OR_598_adj #  1.011382 
OR_629_adj #  0.8648989 







