#reviewer analysis 
#first calculate models and then work out AUC and bootstrp confidence intervals 
#setwd("~/re_gecip/renal/oalavijeh/projects/stones/prs/ishan/")

library(pscl, lib.loc = "~/re_gecip/renal/Rpackages/")
library(pwr, lib.loc = "~/re_gecip/renal/Rpackages/")
library(msm, lib.loc = "~/re_gecip/renal/Rpackages/")
library(ggplot2)
library(dplyr, lib.loc = "~/re_gecip/renal/Rpackages/")
library(logistf)
library(data.table, lib.loc = "~/re_gecip/renal/Rpackages/")
library(FSA)
library(oddsratio, lib.loc = "~/re_gecip/renal/Rpackages/")
library(msm)
library(ggpubr)
options(bitmapType="cairo")

options(scipen=999)

# Calculate liability threshold heritability
#K=pop prevalence
#P=proportion of cases in study 374/24930
#hsq=Heritability estimate (on observed scale)
#bigT = liability threshold
#tau = density of gaussian

h2_liab <- function(x){
  K=0.10
  P=0.01500201
  zv <- dnorm(qnorm(K))
  x * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
}

# Read in the phenotype files 
score <- read.table("~/re_gecip/renal/oalavijeh/projects/stones/prs/ishan/ishan_stone_prs.sscore", header = T, comment.char = '')
cohort <- read.table("~/re_gecip/renal/oalavijeh/projects/stones/cohorts_pheno/final_hes_stone_rdcontrol_ancestry_matched_pc.ped", header = T)
solved_and_slc34a3 <- read.table("~/re_gecip/renal/oalavijeh/projects/stones/cohorts_pheno/solved_slc34a3_platekeys.txt")
solved_only <- fread("~/re_gecip/renal/oalavijeh/projects/stones/analysis/solved_stones.txt")
solved_only <- as.data.frame(unique(solved_only$plate_key))
colnames(solved_only) <- "V1"
slc34a3_only <- fread("~/re_gecip/renal/oalavijeh/projects/stones/analysis/SLC34A3_hes_stones_maf0.001_missense_cadd_variant_GTs_cons_hpo_concat.csv")
slc34a3_only <- as.data.frame(unique(slc34a3_only$plate_key))
colnames(slc34a3_only) <- "V1"
control_slc34a3 <- fread("~/re_gecip/renal/oalavijeh/projects/stones/analysis/SLC34A3_hes_stones_missense_controls_GTs.txt")
control_slc34a3 <- as.data.frame(control_slc34a3$V5)
colnames(control_slc34a3) <- "V1"
#make pheno.prs which is the base pheno file 
pheno.prs <- merge(cohort, score, by.x = 'IID', by.y = 'X.IID', all.x = T, all.y = F)
pheno.prs <- as.data.frame(pheno.prs)
#standardise
control.pheno.prs <- subset(pheno.prs, pheno.prs$PHENO==0)
pheno.prs$prs_standardised <- (pheno.prs$ldpred_beta_SUM-mean(control.pheno.prs$ldpred_beta_SUM))/sd(control.pheno.prs$ldpred_beta_SUM)
#make pheno files for analysis
pheno.prs <- pheno.prs[,c(1:4,6,19,5,7:16)]
pheno.prs$prs_standardised <- as.numeric(pheno.prs$prs_standardised)

#make unsolved which has no solved, no case or control slc34a3 
unsolved.prs <- subset(pheno.prs, !pheno.prs$IID %in% control_slc34a3$V1)
unsolved.prs <- subset(unsolved.prs, !unsolved.prs$IID %in% solved_and_slc34a3$V1)

gel_unsolved.prs <- subset(pheno.prs, !pheno.prs$IID %in% control_slc34a3$V1)
gel_unsolved.prs <- subset(gel_unsolved.prs, !gel_unsolved.prs$IID %in% solved_only$V1)

#modelling, first with all cohort then with unsolved 
# We can then calculate the null model (model with PRS) using logistic regression 
formodel.pheno.prs <- pheno.prs[,c(5,6,7,8:17)]
for_null_model.pheno.prs <- formodel.pheno.prs[,c(1,3:13)]
null.model <- glm(PHENO~., family=binomial(link='logit'),data=for_null_model.pheno.prs, maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(null.model)
prs.result <- NULL
#fit all patients model
model <- glm(PHENO~., family=binomial(link='logit'),data=formodel.pheno.prs, maxit=1000)

# model R2 is obtained as 
model.r2 <- pR2(model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
prs.h2_liab <- h2_liab(prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.coef <- summary(model)$coeff["prs_standardised",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

# Best result is:
prs.result[which.max(prs.result$R2),]
#
#modelling, now with unsolved 
# We can then calculate the null model (model with PRS) using logistic regression 


for_null_unsolved <- unsolved.prs[,c(5,7,8:17)]
unsolved_null.model <- glm(PHENO~., family=binomial(link='logit'),
                  data=for_null_unsolved,
                  maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(unsolved_null.model)
prs.result <- NULL
# Now perform a linear regression on phenotype with PRS and the covariates
# ignoring the FID and IID from our model 
formodel.unsolved.prs <- unsolved.prs[,c(5,6,7,8:17)]
unsolved_model <- glm(PHENO~., family=binomial(link='logit'),data=formodel.unsolved.prs, maxit=1000)
# model R2 is obtained as 
unsolved_model.r2 <- pR2(unsolved_model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
unsolved.prs.r2 <- unsolved_model.r2-null.r2
unsolved.prs.h2_liab <- h2_liab(unsolved.prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.coef <- summary(unsolved_model)$coeff["prs_standardised",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
unsolved_prs.result <- rbind(prs.result, data.frame(R2=unsolved.prs.r2, H2=unsolved.prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

# Best result is:
unsolved_prs.result[which.max(unsolved_prs.result$R2),]


options(scipen=999)


#AUC 
library(DescTools)
library(boot)
Cstat(slc34a3_model)
#bootstrap confidence intervals 
FUN <- function(x,i) {
  r.glm <- glm(PHENO ~ ., data=x[i,], family=binomial)
  Cstat(r.glm)
}
boot.res <- boot(formodel.unsolved.prs, FUN, R=999) 
boot.ci(boot.res, type="perc")

# the percentile confidence intervals
boot.ci(boot.res, type="perc")

#add in whether SLC34A3 monogenic issue 0/1, and another column SLC34AA3*PRS
independence.prs <- pheno.prs 
#standardize scores 
#control.unsolved.prs <- subset(unsolved.prs, unsolved.prs$PHENO==0)
#independence.prs$prs_standardised <- (independence.prs$ldpred_beta_SUM-mean(control.unsolved.prs$ldpred_beta_SUM))/sd(control.unsolved.prs$ldpred_beta_SUM)

independence.prs$slc34a3 <- 0
independence.prs$slc34a3[independence.prs$FID %in% slc34a3_only$V1] <- 1
independence.prs$slc34a3[independence.prs$FID %in% control_slc34a3$V1] <- 1
independence.prs$combined <- independence.prs$prs_standardised*independence.prs$slc34a3
#regression analysis 
fornullmodel.independence.prs <- independence.prs[,c(5,7,8:17)]
formodel.independence.prs <- independence.prs[,c(5,6,18,19,7,8:17)]

#null model 
independence_null.model <- glm(PHENO~., family=binomial(link='logit'),
                           data=fornullmodel.independence.prs,
                           maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(independence_null.model)
#prs.result <- NULL


independence.model <- glm(PHENO~., family=binomial(link='logit'),data=formodel.independence.prs, maxit=10000)
independence.model.r2 <- pR2(independence.model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
independence_prs.r2 <- independence.model.r2-null.r2
independence_prs.h2_liab <- h2_liab(independence_prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
independence_prs.coef <- summary(independence.model)$coeff["prs_standardised",]
prs.beta <- as.numeric(independence_prs.coef[1])
prs.se <- as.numeric(independence_prs.coef[2])
prs.p <- as.numeric(independence_prs.coef[4])
# We can then store the results
prs.result_independence <- rbind(prs.result, data.frame(R2=independence_prs.r2, H2=independence_prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

#AUC
Cstat(independence.model)
FUN <- function(x,i) {
  r.glm <- glm(PHENO ~ ., data=x[i,], family=binomial)
  Cstat(r.glm)
}
boot.res <- boot(independence.prs, FUN, R=10)
boot.ci(boot.res, type="perc")

# the percentile confidence intervals
boot.ci(boot.res, type="perc")



#AUC 
library(DescTools)
library(boot)
Cstat(slc34a3_model)
#bootstrap confidence intervals 
FUN <- function(x,i) {
  r.glm <- glm(PHENO ~ ., data=x[i,], family=binomial)
  Cstat(r.glm)
}
boot.res <- boot(formodel.unsolved.prs, FUN, R=999) 
boot.ci(boot.res, type="perc")

# the percentile confidence intervals
boot.ci(boot.res, type="perc")

#calculate SD increase in PRS
sd_increase_in_prs <- exp(unsolved_model$coefficients["ldpred_beta_SUM"]*sd(unsolved.prs$ldpred_beta_SUM))
sd_increase_in_prs
#95% CI
five_percent <- exp(confint(unsolved_model)["ldpred_beta_SUM",1]*sd(unsolved.prs$ldpred_beta_SUM))
ninety_five_percent <- exp(confint(unsolved_model)["ldpred_beta_SUM",2]*sd(unsolved.prs$ldpred_beta_SUM))

#plot centile of model as prediction of getting USD
library(jtools)
library(OddsPlotty)
library(dplyr)
unsolved_plotting <- formodel.unsolved.prs
unsolved_plotting$rank <- rank(unsolved_plotting$prs_standardised)
unsolved_plotting$rank <- round(unsolved_plotting$rank/length(unsolved_plotting$rank) * 100)
x <- unsolved_plotting %>% group_by(rank) %>% mutate(prevelance = (sum(PHENO==1)/length(PHENO)*100))
l <- x[,c(14,15)]
l <- unique(l)
l <- l[order(l$rank,decreasing = F),]
ggplot(l, aes(x=rank, y=prevelance)) + geom_smooth()
#fit model with rank 
unsolved_plotting_centile <- unsolved_plotting[,c(1,14,3:13)]
unsolved_centile <- glm(PHENO~., family=binomial(link='logit'),data=unsolved_plotting_centile, maxit=10000)
effect_plot(data = unsolved_plotting_centile, unsolved_centile, pred=rank, interval = T, x.label = "PRS Centile",
            y.label = "Prevelance of USD")
#plot_summs(unsolved_model, model, coefs=coef_names)

#older method
#plot by decile of PRS
slc34a3_case_only <- subset(pheno.prs, pheno.prs$IID %in% slc34a3_only$V1)
decile.prs <- formodel.unsolved.prs %>%
  mutate(
    decile = ntile(`prs_standardised`, 100),
    decile = as.factor(decile),
    `PHENO` = as.factor(`PHENO`)
  )

#turn pheno.prs into covariates and rank only
decile.prs <- decile.prs[,c(1,14,3:13)]
# Fit regression model
prs_glm <- glm(`PHENO` ~.,
               data = decile.prs,
               family = 'binomial')

# Put results in data.frame
summs <- prs_glm %>% summary()

# Get point estimates and SEs
results <- bind_cols(coef(prs_glm),
                     summs$coefficients[, 2]) %>%
  setNames(c("estimate", "se"))  %>%
  #slice(1:10) %>%
  mutate(decile = 1:111)

# Your coefficients are on the log odds scale, such that a coefficient is
# log(odds_y==1 / odds_y == 0). We can exponentiate to get odds instead.
results_odds <- results %>% mutate(across(.cols = -decile,
                                          ~ exp(.x)))

# Need SEs on the odds scale too
results_odds <- results_odds %>%
  mutate(var_diag = diag(vcov(prs_glm)),
         se = sqrt(estimate ^ 2 * var_diag))

# Plot with +/- 1 SE
ggplot(results_odds, aes(x = as.factor(decile), y = estimate, color)) +
  geom_point(stat = "identity", shape = 15) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width = 0.4) +
  xlab("PRS Decile") +
  ylab("Odds")

ggplot(results_odds[c(1:100),], aes(x = decile, y = estimate, color)) +
  geom_smooth(method="lm", alpha=.15) +
  xlab("PRS Decile") +
  ylab("Odds")

#Plotting ratios of coefficients on odds scale

estmean <- coef(prs_glm)
estvar <- vcov(prs_glm)
se_dec1 <- deltamethod(~ exp(x1 / x5), estmean, estvar)
se_dec2 <- deltamethod(~ exp(x2 / x5), estmean, estvar)

#But for now we'll omit the SEs since it's not clear what "should" be used:


# Get 5th decile's odds
odds_dec5 <- results_odds %>%
  filter(decile == 5) %>%
  pull(estimate)

# Calculate odds ratios
results_or <- results_odds %>%
  mutate(or_dec5 = estimate / odds_dec5)
#95% CI of each OR
results_or$centile95th <- results_or$estimate+(1.96*results_or$se)
results_or$centile5th <- results_or$estimate-(1.96*results_or$se)

# Plot OR with no SEs
ggplot(results_or[c(1:100),], aes(x = as.factor(decile), y = or_dec5, color)) +
  geom_point(stat = "identity", shape = 15) +
  geom_smooth() + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = centile5th, ymax = centile95th), width = 0.4) +
  xlab("PRS Decile") +
  ylab("Odds Ratio")

#plot smooth line through data
ggplot(results_or[c(1:100),], aes(x=decile, y= or_dec5)) + 
  geom_smooth(method="glm", alpha=.15) +
  #geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  xlab("PRS Decile") +
  ylab("Odds Ratio") +
  expand_limits(x=0, y=0) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,100,by=10))

#stratify PRS by ethnicity 
ancestry <- fread("/re_gecip/renal/oalavijeh/projects/phenotypes/final_phenotypes/all_agg_ancestry.tsv")
ancestry.prs <- merge(pheno.prs, mysql[,c(2,5:9)], by.x = "IID", by.y = "plate_key", all.x = T, all.y = F)
# Create a new column called "ancestry"
ancestry.prs$ancestry <- ifelse(ancestry.prs$pred_european_ancestries >= 0.9, "european",
                         ifelse(ancestry.prs$pred_african_ancestries >= 0.9, "african",
                                ifelse(ancestry.prs$pred_american_ancestries >= 0.9, "american",
                                       ifelse(ancestry.prs$pred_east_asian_ancestries >= 0.9, "east_asian", 
                                              ifelse(ancestry.prs$pred_south_asian_ancestries >= 0.9, "south_asian", NA)))))

ancestry.prs$ancestry[is.na(ancestry.prs$ancestry)] <- "mixed"
table(ancestry.prs$ancestry)

ancestry.prs<- ancestry.prs[,c(1:7,24,8:17)]
write.table(ancestry.prs, "../../../../../cohorts_pheno/ancestry_identified_total_stone_cohort.ped", col.names = T, row.names = F, quote = F, sep = '\t')

#do PRS per ethnic group 
ancestry.prs <- fread("~/re_gecip/renal/oalavijeh/projects/stones/cohorts_pheno/ancestry_identified_total_stone_cohort.ped")
european.prs <- subset(ancestry.prs, ancestry.prs$ancestry=="european")
east_asian.prs <- subset(ancestry.prs, ancestry.prs$ancestry=="east_asian")
south_asian.prs <- subset(ancestry.prs, ancestry.prs$ancestry=="south_asian")
mixed.prs <- subset(ancestry.prs, ancestry.prs$ancestry=="mixed") 
african.prs <- subset(ancestry.prs, ancestry.prs$ancestry=="african") 
european.prs.formodel <- european.prs[,c(7,5,6,9:18)]
east_asian.prs.formodel <- east_asian.prs[,c(7,5,6,9:18)]
south_asian.prs.formodel <- south_asian.prs[,c(7,5,6,9:18)]
mixed.prs.formodel <- mixed.prs[,c(7,5,6,9:18)]
african.prs.formodel <- african.prs[,c(7,5,6,9:18)]

european.prs.fornullmodel <- european.prs[,c(7,6,9:18)]
east_asian.prs.fornullmodel <- east_asian.prs[,c(7,6,9:18)]
south_asian.prs.fornullmodel <- south_asian.prs[,c(7,6,9:18)]
mixed.prs.fornullmodel <- mixed.prs[,c(7,6,9:18)]
african.prs.fornullmodel <- african.prs[,c(7,6,9:18)]

#run PRS models 
european.null.model <- glm(PHENO~., family=binomial(link='logit'),data=european.prs.fornullmodel, maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(european.null.model)
prs.result <- NULL
#fit all patients model
european.model <- glm(PHENO~., family=binomial(link='logit'),data=european.prs.formodel, maxit=1000)

# model R2 is obtained as 
model.r2 <- pR2(european.model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
prs.h2_liab <- h2_liab(prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.coef <- summary(european.model)$coeff["ldpred_beta_SUM",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
european.prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

#run PRS models 
south_asian.null.model <- glm(PHENO~., family=binomial(link='logit'),data=south_asian.prs.fornullmodel, maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(south_asian.null.model)
prs.result <- NULL
#fit all patients model
south_asian.model <- glm(PHENO~., family=binomial(link='logit'),data=south_asian.prs.formodel, maxit=1000)

# model R2 is obtained as 
model.r2 <- pR2(south_asian.model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
prs.h2_liab <- h2_liab(prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.coef <- summary(south_asian.model)$coeff["ldpred_beta_SUM",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
south_asian.prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

#run PRS models 
east_asian.null.model <- glm(PHENO~., family=binomial(link='logit'),data=south_asian.prs.fornullmodel, maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(east_asian.null.model)
prs.result <- NULL
#fit all patients model
east_asian.model <- glm(PHENO~., family=binomial(link='logit'),data=east_asian.prs.formodel, maxit=1000)

# model R2 is obtained as 
model.r2 <- pR2(east_asian.model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
prs.h2_liab <- h2_liab(prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.coef <- summary(east_asian.model)$coeff["ldpred_beta_SUM",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
east_asian.prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

#run PRS models 
african.null.model <- glm(PHENO~., family=binomial(link='logit'),data=african.prs.fornullmodel, maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(african.null.model)
prs.result <- NULL
#fit all patients model
african.model <- glm(PHENO~., family=binomial(link='logit'),data=african.prs.formodel, maxit=1000)

# model R2 is obtained as 
model.r2 <- pR2(african.model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
prs.h2_liab <- h2_liab(prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.coef <- summary(african.model)$coeff["ldpred_beta_SUM",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
african.prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

#run PRS models 
mixed.null.model <- glm(PHENO~., family=binomial(link='logit'),data=mixed.prs.fornullmodel, maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(mixed.null.model)
prs.result <- NULL
#fit all patients model
mixed.model <- glm(PHENO~., family=binomial(link='logit'),data=mixed.prs.formodel, maxit=1000)

# model R2 is obtained as 
model.r2 <- pR2(mixed.model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
prs.h2_liab <- h2_liab(prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.coef <- summary(mixed.model)$coeff["ldpred_beta_SUM",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
mixed.prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

#run PRS models for both asian 
asian.prs <- rbind(east_asian.prs, south_asian.prs)
asian.prs.formodel <- asian.prs[,c(7,5,6,9:18)]
asian.prs.fornullmodel <- asian.prs[,c(7,6,9:18)]
asian.null.model <- glm(PHENO~., family=binomial(link='logit'),data=asian.prs.fornullmodel, maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(asian.null.model)
prs.result <- NULL
#fit all patients model
asian.model <- glm(PHENO~., family=binomial(link='logit'),data=asian.prs.formodel, maxit=1000)

# model R2 is obtained as 
model.r2 <- pR2(asian.model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
prs.h2_liab <- h2_liab(prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.coef <- summary(asian.model)$coeff["ldpred_beta_SUM",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
asian.prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

#african without PCs 
african.prs.formodel.nopc <- african.prs.formodel[,c(1:3)]
african.nopc.model <- glm(PHENO~., family=binomial(link='logit'),data=african.prs.formodel.nopc, maxit=1000)



#AUC 
library(DescTools)
library(boot)
Cstat(slc34a3_model)
#bootstrap confidence intervals 
FUN <- function(x,i) {
  r.glm <- glm(PHENO ~ ., data=x[i,], family=binomial)
  Cstat(r.glm)
}
boot.res <- boot(formodel.unsolved.prs, FUN, R=999) 
boot.ci(boot.res, type="perc")

# the percentile confidence intervals
boot.ci(boot.res, type="perc")

#calculate SD increase in PRS
sd_increase_in_prs <- exp(african.model$coefficients["prs.standardised"]*sd(african.prs.formodel$prs.standardised))
sd_increase_in_prs
#95% CI
five_percent <- exp(confint(african.nopc.model)["ldpred_beta_SUM",1]*sd(african.prs.formodel.nopc$ldpred_beta_SUM))
five_percent
ninety_five_percent <- exp(confint(african.nopc.model)["ldpred_beta_SUM",2]*sd(african.prs.formodel.nopc$ldpred_beta_SUM))
ninety_five_percent

#metanalysis of odds ratios from ancestry 
library(metafor)
#meta anlaysis of PRS per ethnicity
coef1 <- coef(european.model)["ldpred_beta_SUM"]
SE1 <- sqrt(vcov(european.model)["ldpred_beta_SUM", "ldpred_beta_SUM"])

coef2 <- coef(african.model)["ldpred_beta_SUM"]
SE2 <- sqrt(vcov(african.model)["ldpred_beta_SUM", "ldpred_beta_SUM"])

coef3 <- coef(asian.model)["ldpred_beta_SUM"]
SE3 <- sqrt(vcov(asian.model)["ldpred_beta_SUM", "ldpred_beta_SUM"])

coef4 <- coef(mixed.model)["ldpred_beta_SUM"]
SE4 <- sqrt(vcov(mixed.model)["ldpred_beta_SUM", "ldpred_beta_SUM"])

# Calculate weight for each study based on sample size
n1 <- sum(length(european.prs$ldpred_beta_SUM))
n2 <- sum(!is.na(african.prs$ldpred_beta_SUM))
n3 <- sum(!is.na(asian.prs$ldpred_beta_SUM))
n4 <- sum(!is.na(mixed.prs$ldpred_beta_SUM))


w1 <- n1 / (n1 + n2 + n3 + n4)
w2 <- n2 / (n1 + n2 + n3 + n4)
w3 <- n3 / (n1 + n2 + n3 + n4)
w4 <- n4 / (n1 + n2 + n3 + n4)

# Calculate standard deviation of predictor variable for each dataset
sd_predictor1 <- sd(european.prs$ldpred_beta_SUM)
sd_predictor2 <- sd(african.prs$ldpred_beta_SUM)
sd_predictor3 <- sd(asian.prs$ldpred_beta_SUM)
sd_predictor4 <- sd(mixed.prs$ldpred_beta_SUM)

# Calculate one standard deviation increase in predictor variable for each dataset
one_sd_increase1 <- sd_predictor1
one_sd_increase2 <- sd_predictor2
one_sd_increase3 <- sd_predictor3
one_sd_increase4 <- sd_predictor4

# Calculate change in log odds for one standard deviation increase in predictor for each dataset
change_in_log_odds1 <- coef1 * one_sd_increase1
change_in_log_odds2 <- coef2 * one_sd_increase2
change_in_log_odds3 <- coef3 * one_sd_increase3
change_in_log_odds4 <- coef4 * one_sd_increase4

# Calculate odds ratio for one standard deviation increase in predictor for each dataset
odds_ratio1 <- exp(change_in_log_odds1)
odds_ratio2 <- exp(change_in_log_odds2)
odds_ratio3 <- exp(change_in_log_odds3)
odds_ratio4 <- exp(change_in_log_odds4)



# Calculate standard error of odds ratio for one standard deviation increase in predictor for each dataset
SE_odds_ratio1 <- abs(log(exp(1))) * SE1 * one_sd_increase1
SE_odds_ratio2 <- abs(log(exp(1))) * SE2 * one_sd_increase2
SE_odds_ratio3 <- abs(log(exp(1))) * SE3 * one_sd_increase3
SE_odds_ratio4 <- abs(log(exp(1))) * SE4 * one_sd_increase4

# Create data frame with effect size, standard error, and weight for each study
dat <- data.frame(effect = c(odds_ratio1, odds_ratio2, odds_ratio3, odds_ratio4), SE = c(SE_odds_ratio1, SE_odds_ratio2, SE_odds_ratio3, SE_odds_ratio4), w = c(w1, w2,w3,w4))
dat$study <- c("European (n=19420)", "African (n=559)", "Asian (n=3035)","Mixed (n=2221)")
# Meta-analyze odds ratio using random-effects model
res <- rma(yi = effect, sei = SE, weights = w, data = dat, )

# Print summary of meta-analysis results and forest plot
summary(res)
labs <- dat$study
options(scipen = 0)
forest(res, olim = c(0,2), refline = 1,xlab = "Odds Ratio", slab=labs, mlab = "Meta-analysis")
mtext(paste("Association p-value=",signif(summary(res)$pval,digits = 3)),side=3, line=-1)
mtext(paste("Heterogeneity p-value=",signif(summary(res)$QEp,digits=3)),side=3, line=-2.25)

#top 5% using one SD as reference
ldpred_beta_SUM_sd <- sd(formodel.unsolved.prs$ldpred_beta_SUM)
ldpred_beta_SUM_95th <- quantile(formodel.unsolved.prs$ldpred_beta_SUM, 0.95)
formodel.unsolved.prs$top_5_percent <- ifelse(formodel.unsolved.prs$ldpred_beta_SUM > ldpred_beta_SUM_95th, 1, 0)
for_fitting <- formodel.unsolved.prs[,c(1:13)]
five_percent_model <- glm(PHENO~.,data=for_fitting,family = "binomial")
OR_top_5_percent_vs_1SD <- exp((coef(five_percent_model)["ldpred_beta_SUM"] * ldpred_beta_SUM_95th) - (0.5 * coef(five_percent_model)["ldpred_beta_SUM"] * ldpred_beta_SUM_sd))

#top1%  using one SD as reference
ldpred_beta_SUM_sd <- sd(formodel.unsolved.prs$ldpred_beta_SUM)
ldpred_beta_SUM_99th <- quantile(formodel.unsolved.prs$ldpred_beta_SUM, 0.99)
formodel.unsolved.prs$top_1_percent <- ifelse(formodel.unsolved.prs$ldpred_beta_SUM > ldpred_beta_SUM_99th, 1, 0)
for_fitting <- formodel.unsolved.prs[,c(1:13)]
one_percent_model <- glm(PHENO~.,data=for_fitting,family = "binomial")
OR_top_1_percent_vs_1SD <- exp((coef(one_percent_model)["ldpred_beta_SUM"] * ldpred_beta_SUM_99th) - (0.5 * coef(one_percent_model)["ldpred_beta_SUM"] * ldpred_beta_SUM_sd))

#top 1% compared to other 99%
top1percent.prs<- formodel.unsolved.prs
top1percent.prs$top_centile <- ifelse(top1percent.prs$ldpred_beta_SUM >= quantile(top1percent.prs$ldpred_beta_SUM, probs = 0.99), 1, 0)
top1percent.prs <- top1percent.prs[,c(1,14,3:13)]
model <- glm(PHENO ~., data = top1percent.prs, family = binomial)
odds_ratio <- exp(coef(model)["top_centile"])

#top 5% compared to other 95%
top5percent.prs<- formodel.unsolved.prs
top5percent.prs$top_centile <- ifelse(top5percent.prs$ldpred_beta_SUM >= quantile(top5percent.prs$ldpred_beta_SUM, probs = 0.95), 1, 0)
top5percent.prs <- top5percent.prs[,c(1,14,3:13)]
model <- glm(PHENO ~., data = top5percent.prs, family = binomial)
odds_ratio <- exp(coef(model)["top_centile"])

#standardise data to control data
unsolved.prs.standardised <- unsolved.prs
control.unsolved.prs <- subset(unsolved.prs, unsolved.prs$PHENO==0)
unsolved.prs.standardised$prs_standardised <- (unsolved.prs.standardised$ldpred_beta_SUM-mean(control.unsolved.prs$ldpred_beta_SUM))/sd(control.unsolved.prs$ldpred_beta_SUM)
hist(unsolved.prs.standardised$prs_standardised)
#modelling, now with unsolved 
# We can then calculate the null model (model with PRS) using logistic regression 
for_null_unsolved_standardised <- unsolved.prs.standardised[,c(7,6,8:17)]
unsolved_null.model <- glm(PHENO~., family=binomial(link='logit'),
                           data=for_null_unsolved_standardised,
                           maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(unsolved_null.model)
prs.result <- NULL
# Now perform a linear regression on phenotype with PRS and the covariates
# ignoring the FID and IID from our model 
formodel.unsolved.prs.standardised <- unsolved.prs.standardised[,c(7,19,6,8:17)]
unsolved_model <- glm(PHENO~., family=binomial(link='logit'),data=formodel.unsolved.prs.standardised, maxit=1000)
# model R2 is obtained as 
unsolved_model.r2 <- pR2(unsolved_model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
unsolved.prs.r2 <- unsolved_model.r2-null.r2
unsolved.prs.h2_liab <- h2_liab(unsolved.prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.coef <- summary(unsolved_model)$coeff["prs_standardised",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
unsolved_prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))

# Best result is:
unsolved_prs.result[which.max(unsolved_prs.result$R2),]


options(scipen=999)


#AUC 
library(DescTools)
library(boot)
Cstat(slc34a3_model)
#bootstrap confidence intervals 
FUN <- function(x,i) {
  r.glm <- glm(PHENO ~ ., data=x[i,], family=binomial)
  Cstat(r.glm)
}
boot.res <- boot(formodel.unsolved.prs, FUN, R=999) 
boot.ci(boot.res, type="perc")

# the percentile confidence intervals
boot.ci(boot.res, type="perc")

#calculate SD increase in PRS
sd_increase_in_prs <- exp(unsolved_model$coefficients["prs_standardised"]*sd(unsolved.prs.standardised$prs_standardised))
sd_increase_in_prs
#95% CI
five_percent <- exp(confint(unsolved_model)["ldpred_beta_SUM",1]*sd(unsolved.prs$ldpred_beta_SUM))
ninety_five_percent <- exp(confint(unsolved_model)["ldpred_beta_SUM",2]*sd(unsolved.prs$ldpred_beta_SUM))

#plot standardies graphs and distributions 
#create another pheno cohort to subdivide solved, slc34a3 and controls that have slc34a3 
subgroup.pheno.prs <- pheno.prs
#standardize scores 
control.unsolved.prs <- subset(unsolved.prs, unsolved.prs$PHENO==0)
subgroup.pheno.prs$prs_standardised <- (subgroup.pheno.prs$ldpred_beta_SUM-mean(control.unsolved.prs$ldpred_beta_SUM))/sd(control.unsolved.prs$ldpred_beta_SUM)
hist(subgroup.pheno.prs$prs_standardised)

subgroup.pheno.prs$PHENO <- as.character(subgroup.pheno.prs$PHENO)
subgroup.pheno.prs$PHENO[subgroup.pheno.prs$IID %in% slc34a3_only$V1] <- "Case SLC34A3"
subgroup.pheno.prs$PHENO[subgroup.pheno.prs$IID %in% control_slc34a3$V1] <- "Control SLC34A3"
subgroup.pheno.prs$PHENO[subgroup.pheno.prs$PHENO =="0"] <- "Control"
subgroup.pheno.prs$PHENO[subgroup.pheno.prs$PHENO =="1"] <- "Unsolved"
subgroup.pheno.prs$PHENO <- as.factor(subgroup.pheno.prs$PHENO)
subgroup.pheno.prs$PHENO <- reorder(subgroup.pheno.prs$PHENO, subgroup.pheno.prs$prs_standardised, mean)

#plotting for subgroup 
# Start plotting
ggplot(subgroup.pheno.prs, aes(x=prs_standardised, fill=PHENO))+
  geom_density(alpha=0.3)+
  theme_classic()+
  labs(x="Polygenic Score", y="Density", fill="Cohort") 
#ggsave("slc34a3_subgroup_prs_stones.density.png", height = 7, width = 7)

# Violin plot to compare phenotypes
ggplot(subgroup.pheno.prs, aes(x=PHENO, y=prs_standardised, fill=PHENO)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") +
  theme(axis.text.x = element_text(angle=30, vjust=0.6, hjust=0.5))
#ggsave("slc34a3_subgroup_prs_stones.violin.png", height = 7, width = 10)
#stats on subgroup 
group_by(subgroup.pheno.prs, PHENO) %>%
  summarise(
    count = n(),
    mean = mean(prs_standardised, na.rm = TRUE),
    median = median(prs_standardised, na.rm = TRUE),
    sd = sd(prs_standardised, na.rm = TRUE)
  )


#comparing two distributions 
kruskal.test(prs_standardised ~ PHENO, data=subgroup.pheno.prs)
pairwise.wilcox.test(subgroup.pheno.prs$prs_standardised, subgroup.pheno.prs$PHENO, p.adjust.method = "BH")
#pop <- aov(ldpred_beta_SUM ~ PHENO, data=subgroup.pheno.prs)
#summary(pop)
#t.test(x$ldpred_beta_SUM, y$ldpred_beta_SUM)

#plot with stats 
my_comparisons <- list(c("Unsolved", "Control"), c("Case SLC34A3", "Control"), c("Case SLC34A3","Unsolved"),
                       c("Case SLC34A3", "Control SLC34A3"))
ggplot(subgroup.pheno.prs, aes(x=PHENO, y=prs_standardised, fill=PHENO)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  scale_fill_brewer(palette = "Dark2") +
  ylim(-4,7) + 
  theme_classic() +
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") +
  theme(axis.text.x = element_text(angle=30, vjust=0.6, hjust=0.5)) +
  stat_compare_means(comparisons = my_comparisons) + 
  stat_compare_means(label.y =50)


ggsave("slc34a3_subgroup_prs_stones.violin.png", height = 7, width = 10)
options(scipen = 0)
ggviolin(subgroup.pheno.prs, x = "PHENO", y = "prs_standardised", fill = "PHENO",
         palette = c("#00AFBB", "#28def0", "#FC4E07","#f2241f"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 8) + 
  stat_cor(p.accuracy = 0.01) +
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") 
ggsave("slc34a3_subgroup_prs_stones_ggpubr_4groups.violin.png", height = 7, width = 10)

#plot PRS distributions per ethnicity 
#first normalise each to specific ancestry controls 
asian.prs.controls <- subset(asian.prs,asian.prs$PHENO==0)
asian.prs$prs.standardised <- (asian.prs$ldpred_beta_SUM-mean(asian.prs.controls$ldpred_beta_SUM))/sd(asian.prs.controls$ldpred_beta_SUM)
asian.prs$cohort <- "Asian"
african.prs.controls <- subset(african.prs,african.prs$PHENO==0)
african.prs$prs.standardised <- (african.prs$ldpred_beta_SUM-mean(african.prs.controls$ldpred_beta_SUM))/sd(african.prs.controls$ldpred_beta_SUM)
african.prs$cohort <- "African"
european.prs.controls <- subset(european.prs,european.prs$PHENO==0)
european.prs$prs.standardised <- (european.prs$ldpred_beta_SUM-mean(european.prs.controls$ldpred_beta_SUM))/sd(european.prs.controls$ldpred_beta_SUM)
european.prs$cohort <- "European"
mixed.prs.controls <- subset(mixed.prs,mixed.prs$PHENO==0)
mixed.prs$prs.standardised <- (mixed.prs$ldpred_beta_SUM-mean(mixed.prs.controls$ldpred_beta_SUM))/sd(mixed.prs.controls$ldpred_beta_SUM)
mixed.prs$cohort <- "Mixed "
#then rebind and plot 
ancestry_defined <- rbind(asian.prs, african.prs, european.prs, mixed.prs)
ggplot(ancestry_defined, aes(x=prs.standardised, fill=cohort))+
  geom_density(alpha=0.3)+
  theme_classic()+
  labs(x="Polygenic Score", y="Density", fill="Cohort") 
#try and plot as facet 
ggplot(ancestry_defined, aes(x=ldpred_beta_SUM, fill=cohort))+
  geom_density(alpha=0.3)+
  theme_classic()+
  labs(x="Polygenic Score", y="Density", fill="Cohort") +
  facet_grid(rows=vars(cohort)) 

#plot PRS distributions per ethnicity without normalising 
#first normalise each to specific ancestry controls 
asian.prs.controls <- subset(asian.prs,asian.prs$PHENO==0)
asian.prs$prs.standardised <- (asian.prs$ldpred_beta_SUM-mean(asian.prs.controls$ldpred_beta_SUM))/sd(asian.prs.controls$ldpred_beta_SUM)
asian.prs$cohort <- "Asian"
african.prs.controls <- subset(african.prs,african.prs$PHENO==0)
african.prs$prs.standardised <- (african.prs$ldpred_beta_SUM-mean(african.prs.controls$ldpred_beta_SUM))/sd(african.prs.controls$ldpred_beta_SUM)
african.prs$cohort <- "African"
european.prs.controls <- subset(european.prs,european.prs$PHENO==0)
european.prs$prs.standardised <- (european.prs$ldpred_beta_SUM-mean(european.prs.controls$ldpred_beta_SUM))/sd(european.prs.controls$ldpred_beta_SUM)
european.prs$cohort <- "European"
mixed.prs.controls <- subset(mixed.prs,mixed.prs$PHENO==0)
mixed.prs$prs.standardised <- (mixed.prs$ldpred_beta_SUM-mean(mixed.prs.controls$ldpred_beta_SUM))/sd(mixed.prs.controls$ldpred_beta_SUM)
mixed.prs$cohort <- "Mixed "
#then rebind and plot 
ancestry_defined <- rbind(asian.prs, african.prs, european.prs, mixed.prs)
ggplot(ancestry_defined, aes(x=prs.standardised, fill=cohort))+
  geom_density(alpha=0.3)+
  theme_classic()+
  labs(x="Polygenic Score", y="Density", fill="Cohort") 
#try and plot as facet 
ggplot(ancestry_defined, aes(x=prs.standardised, fill=cohort))+
  geom_density(alpha=0.3)+
  theme_classic()+
  labs(x="Polygenic Score", y="Density", fill="Cohort") +
  facet_grid(rows=vars(cohort))

#power calculations 
n <- 25304
alpha <- 0.05
OR_interact <- 2
power <- pwr.f2.test(u=log(OR_interact), v1 =15,v2=25289, sig.level = alpha, power=NULL)

#pairwise calculations 
library(emmeans)
options(scipen = 0)
pairs <- emmeans(independence.model, pairwise ~ prs_standardised + slc34a3 + combined)
pdf("test.pdf")
plot(pairs, type="response")
dev.off()
#calculate SD increase in PRS
sd_increase_in_prs <- exp(independence.model$coefficients["prs_standardised"]*sd(independence.prs$prs_standardised))
sd_increase_in_prs
#95% CI
five_percent <- exp(confint(independence.model)["prs_standardised",1]*sd(independence.prs$prs_standardised))
ninety_five_percent <- exp(confint(independence.model)["prs_standardised",2]*sd(independence.prs$prs_standardised))


#preveleance per cohort plotting
library(jtools)
library(OddsPlotty)
library(dplyr)
all_slc34a3 <- rbind(control_slc34a3, slc34a3_only)
slc34a3_carriers_for_model <- subset(pheno.prs, pheno.prs$IID %in% all_slc34a3$V1)
slc34a3_plotting <- slc34a3_carriers_for_model
slc34a3_plotting$rank <- rank(slc34a3_plotting$prs_standardised)
slc34a3_plotting$rank <- round(slc34a3_plotting$rank/length(slc34a3_plotting$rank) * 100)
y <- slc34a3_plotting %>% group_by(rank) %>% mutate(prevelance = (sum(PHENO==1)/length(PHENO)*100))
p<- y[,c(18,19)]
p <- unique(p)
p <- p[order(p$rank,decreasing = F),]
ggplot(p, aes(x=rank, y=prevelance)) + geom_smooth(method="glm", alpha=.15)
ggplot(p, aes(x=rank, y=prevelance)) + geom_smooth()

#fit model with rank 
slc34a3_plotting_centile <- slc34a3_plotting[,c(5,18,7:17)]
slc34a3_centile <- glm(PHENO~., family=binomial(link='logit'),data=slc34a3_plotting_centile, maxit=10000)
slc34a3_plot<-effect_plot(data = slc34a3_plotting_centile, slc34a3_centile, pred=rank, interval = F, x.label = "PRS Centile",
                          y.label = "Prevelance of USD") + scale_y_continuous(labels = scales::percent)
###all commerts 
unsolved_plotting <- gel_unsolved.prs[,c(5:17)]
unsolved_plotting$rank <- rank(unsolved_plotting$prs_standardised)
unsolved_plotting$rank <- round(unsolved_plotting$rank/length(unsolved_plotting$rank) * 100)
x <- unsolved_plotting %>% group_by(rank) %>% mutate(prevelance = (sum(PHENO==1)/length(PHENO)*100))
l <- x[,c(14,15)]
l <- unique(l)
l <- l[order(l$rank,decreasing = F),]
ggplot(l, aes(x=rank, y=prevelance)) + geom_smooth(method="glm", alpha=.15)
#fit model with rank 
unsolved_plotting_centile <- unsolved_plotting[,c(1,14,3:13)]
unsolved_centile <- glm(PHENO~., family=binomial(link='logit'),data=unsolved_plotting_centile, maxit=10000)
unsolved_plot<-effect_plot(data = unsolved_plotting_centile, unsolved_centile, pred=rank, interval = F, x.label = "PRS Centile",
                           y.label = "Prevelance of USD") + scale_y_continuous(labels = scales::percent)

#combine data for ggplot
data1 <- make_predictions(unsolved_centile, pred="rank")
data1$model <- "SLC34A3 non-carrier"
data2 <- make_predictions(slc34a3_centile, pred="rank")
data2$model <- "SLC34A3 carrier"
cdata <- bind_rows(data1,data2)

ggplot(data = cdata,aes(x=rank,y=PHENO,color=model)) +
  geom_smooth(method = "glm", alpha=0.15, se = T) + 
  labs(x = "Percentile of Polygenic risk score", y = "Prevelance of USD", color = "Model") +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.02)) +
  theme_minimal() + 
  scale_color_discrete(name = "Cohort") 


ggplot() +
  geom_line(data=data1, aes(x=rank, y=PHENO), color='green') + geom_smooth(method="glm", alpha=.15) +
  geom_line(data=data2, aes(x=rank, y=PHENO), color='red') + geom_smooth(method="glm", alpha=.15) +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.02)) 


###########probability plot#######
# Create a new column for prevalence by cohort
df <- pheno.prs
df$cohort <- NULL
df$cohort[df$IID %in% all_slc34a3$V1] <- "SLC34A3 carrier"
df$cohort[!df$IID %in% all_slc34a3$V1] <- "SLC34A3 non-carrier"
df$prevalence <- ifelse(df$PHENO == 1, 1, 0)
prevalence_by_cohort <- aggregate(prevalence ~ prs_standardised + cohort, df, mean)
prevalence_by_cohort$prevalence <- prevalence_by_cohort$prevalence * 100
# Create the plot - facet
ggplot(df, aes(x = prs_standardised, y = prevalence, color = cohort)) + 
  #geom_point() + 
  geom_smooth(method = "glm", se = FALSE, span = 0.5) +
  xlab("PRS") + 
  ylab("Prevalence") + 
  ggtitle("Prevalence of phenotype by PRS and cohort") +
  facet_wrap(~cohort, scales = "free_x")

#single plot
ggplot(df, aes(x = prs_standardised, y = prevalence, color = cohort, group = cohort)) + 
  #geom_point() + 
  geom_smooth(method = "glm", se = FALSE, span = 0.5) +
  xlab("PRS") + 
  ylab("Proportion with USD in 100KGP") + 
  labs(color="Carrier Status")+
  scale_y_continuous(labels = scales::percent) + 
  theme_minimal()
ggsave("analysis/proportion_with_usd_plot_100kgp.png", height = 7, width = 7)

#distrbutions for presentation 
dist <- pheno.prs
dist$PHENO <- ifelse(!dist$IID %in% solved_only$V1 & dist$PHENO==1, "Unsolved Cases", dist$PHENO)
dist$PHENO[dist$PHENO==0] <- "Control"
dist$PHENO[dist$PHENO==1] <- "Case"
dist$PHENO[dist$IID %in% slc34a3_only$V1] <- "Case SLC34A3"
dist$PHENO[dist$IID %in% control_slc34a3$V1] <- "Control SLC34A3"


no_control_slc34a3 <- dist[!dist$PHENO=="Control SLC34A3",]
no_control_slc34a3 <- no_control_slc34a3[!no_control_slc34a3$PHENO=="Case",]
no_control_slc34a3$PHENO <- factor(no_control_slc34a3$PHENO, levels =c("Case SLC34A3", "Unsolved Cases","Control"))
#plot
ggplot(no_control_slc34a3, aes(x=prs_standardised, fill=PHENO))+
  geom_density(alpha=0.2)+
  theme_classic()+
  labs(x="Polygenic Score", y="Density", fill="Cohort") + facet_grid(rows = vars(PHENO)) +
  geom_vline(xintercept = 0, linetype = "dotted")


ggsave("densities_slc34a3_case_unsolved_cases_controls.density.png", height = 7, width = 7)

no_control_slc34a3_mean <- no_control_slc34a3 %>%
  group_by(PHENO) %>%
  summarize(mean_prs = mean(prs_standardised))

no_control_slc34a3_median <- no_control_slc34a3 %>%
  group_by(PHENO) %>%
  summarize(median_prs = median(prs_standardised))

ggplot(no_control_slc34a3, aes(x = prs_standardised, fill = PHENO)) +
  geom_density(alpha = 0.2) +
  geom_vline(data = no_control_slc34a3_median, 
             aes(xintercept = median_prs, group = PHENO), 
             color = "black", 
             size = 1, linetype = "dotted") +
  theme_classic() +
  labs(x = "PRS", y = "Density", fill = "Cohort") +
  facet_grid(rows = vars(PHENO))

ggsave("analysis/densities_slc34a3_case_unsolved_cases_controls.density_median_line.png", height = 7, width = 7)

