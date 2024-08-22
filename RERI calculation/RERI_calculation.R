mydata <- read.xlsx("RERI calculation/new_mirna2.xlsx") # experimental result
colnames(mydata) <- c("pkod","miR132", "miR130", "MDD")
mydata <- mydata[,c("miR132", "miR130", "MDD")]
testn <- mydata
testn$miR132_miR130 <- testn$miR132  * testn$miR130


testnrun_data <- testn[, c("miR132", "miR130",  "miR132_miR130", "MDD")]
View(testnrun_data)
testnrun_data <- na.omit(testnrun_data)

# calculate reri for 98 samples for miR132, miR130, don't scaled 
lm_98sample_test2n <- glm(as.numeric(MDD) ~ ., data = testnrun_data, family = binomial ) # logistic regression
threebeta_98sample_test2n <- lm_98sample_test2n$coefficients[c("miR132", "miR130",  "miR132_miR130")]
reri_98sample_test2n <- exp(sum(threebeta_98sample_test2n)) - exp(threebeta_98sample_test2n[1]) - exp(threebeta_98sample_test2n[2])+1 
reri_98sample_test2n # 0.10141

# randomly select the 90 samples from 98 samples, 
reri_run1000 <- matrix(0,1000,1)  
set.seed(1000)

for (i in 1:1000){
  ind <- sort(sample(1:98, 90, replace = F)) 
  log_fit <- glm(as.numeric(MDD)~., data = testnrun_data[ind,], family = binomial)
  coef <- log_fit$coefficients[c("miR132", "miR130",  "miR132_miR130")]
  reri_run1000[i] <- exp(sum(coef))-exp(coef[1])-exp(coef[2])+1
}
refi_run1000_s <- sort(reri_run1000)
refi_run1000_ci <- c(refi_run1000_s[25], refi_run1000_s[975])
refi_run1000_ci # -0.3428663  0.2802874
# unscaled samples, randomly select 90 samples for 1000 times without duplicate sampling, 
# reri (95%): -0.3428663  0.2802874

