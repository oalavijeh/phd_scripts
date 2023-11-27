#reviewer analysis 
#first calculate models and then work out AUC and bootstrp confidence intervals 
#setwd("~/re_gecip/renal/oalavijeh/projects/stones/prs/ishan/")

library(pscl, lib.loc = "~/re_gecip/renal/Rpackages/")
library(pwr, lib.loc = "~/re_gecip/renal/Rpackages/")
library(ggplot2)
library(dplyr, lib.loc = "~/re_gecip/renal/Rpackages/")
library(logistf)
library(data.table, lib.loc = "~/re_gecip/renal/Rpackages/")
library(FSA)
library(oddsratio, lib.loc = "~/re_gecip/renal/Rpackages/")
library(msm)
library(ggpubr)
#defo need this one for multinomial logistic regression 
library(nnet)
options(bitmapType="cairo")

options(scipen=0)

# Calculate liability threshold heritability
#K=pop prevalence
#P=proportion of cases in study 1209/26096
#hsq=Heritability estimate (on observed scale)
#bigT = liability threshold
#tau = density of gaussian

h2_liab <- function(x){
  K=0.001
  P=0.04632894
  zv <- dnorm(qnorm(K))
  x * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
}
#read in prs
setwd("/re_gecip/renal/oalavijeh/projects/unexplained_esrf/prs/")
# Read in the prs file
ckd_prs <- fread("unexplained_esrf_controls_prs_apol1_status_covariates.txt")
colnames(ckd_prs)[14] <- "ckd_prs"
iga_prs <- fread("unexplained_iga_77snp_prs_scores.sscore")
colnames(iga_prs)[1] <- "id"
colnames(iga_prs)[3] <- "iga_prs"
ssns_prs <- fread("ssns/unexplained_ssns_prs_scores.sscore")
colnames(ssns_prs)[1] <- "id"
colnames(ssns_prs)[3] <- "ssns_prs"
membranous_prs <- fread("membranous/unexplained_membranous_prs_scores.sscore")
colnames(membranous_prs)[1] <- "id"
colnames(membranous_prs)[3] <- "membranous_prs"

#merge all 
prs <- Reduce(merge,list(ckd_prs, iga_prs[,c(1,3)],membranous_prs[,c(1,3)],ssns_prs[,c(1,3)]))
#get rid of monogenic 
monogenic <- fread("/re_gecip/renal/oalavijeh/projects/unexplained_esrf/apol1/all_solved_cases_with_apol1.txt")
prs <- subset(prs, !id %in% monogenic$plate_key)

#standardise for each score
control.prs <- subset(prs, pheno==0)
prs$iga_prs_standardised <- (prs$iga_prs-mean(control.prs$iga_prs))/sd(control.prs$iga_prs)
prs$ssns_prs_standardised <- (prs$ssns_prs-mean(control.prs$ssns_prs))/sd(control.prs$ssns_prs)
prs$membranous_prs_standardised <- (prs$membranous_prs-mean(control.prs$membranous_prs))/sd(control.prs$membranous_prs)
prs$ckd_prs_standardised <- (prs$ckd_prs-mean(control.prs$ckd_prs))/sd(control.prs$ckd_prs)

#make apol1 column 
prs$apol1 <- ifelse(prs$type=="Case - High Rx APOL1"| prs$type=="Control - High Rx APOL1",1,0)
write.table(prs, "../prs/unexplained_case_rdcontrol_iga77_mem_ssns_ckd_prs_standardised_pcs_sex_apol1.txt", col.names = F, row.names = F, quote = F, sep = '\t')
#make pheno files for modelling with variants analysis
formodel.prs <- prs[,c(3,22,19:21,23,2,4:13)]

#modelling - first do cases versus controls only ten multinomial lr
# We can then calculate the null model (model with PRS) using logistic regression 
formodel.null.prs <- formodel.prs[,c(1,7:17)]
null.model <- glm(pheno~., family=binomial(link='logit'),data=formodel.null.prs, maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(null.model)
prs.result <- NULL
#fit all patients model
model <- glm(pheno~., family=binomial(link='logit'),data=formodel.prs, maxit=1000)

# model R2 is obtained as 
model.r2 <- pR2(model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
#prs.h2_liab <- h2_liab(prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.ckd.coef <- summary(model)$coeff["ckd_prs_standardised",]
prs.beta <- as.numeric(prs.ckd.coef[1])
prs.se <- as.numeric(prs.ckd.coef[2])
prs.p <- as.numeric(prs.ckd.coef[4])
# We can then store the results
#prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))
prs.result <- rbind(prs.result, data.frame(R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))

# Best result is:
prs.result[which.max(prs.result$R2),]
#get model data
summary(model)

#increase in various things and OR
#calculate SD increase in apol1
sd_increase_in_prs <- exp(model$coefficients["apol1"]*sd(prs$apol1))
sd_increase_in_prs
#95% CI
five_percent <- exp(confint(model)["apol1",1]*sd(prs$apol1))
ninety_five_percent <- exp(confint(model)["apol1",2]*sd(prs$apol1))

#plot difference between controls and cases by apol1 and other PRS
to_plot <- prs

#plot membranous 
ggplot(to_plot, aes(x = membranous_prs_standardised, fill = factor(pheno))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of PRS in cases and controls",
       x = "PRS",
       y = NULL) +
  theme_minimal() + 
  labs(fill = "pheno")

#plot ckd versus apol1 type 
ggplot(to_plot, aes(x = ckd_prs_standardised, fill = factor(type))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of PRS in cases and controls",
       x = "PRS",
       y = NULL) +
  theme_minimal() + 
  labs(fill = "type")

#plot ssns and iga
ggplot(to_plot, aes(x = membranous_prs_standardised, fill = factor(type))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of PRS in cases and controls",
       x = "PRS",
       y = NULL) +
  theme_minimal() + 
  labs(fill = "type")

#plot box and violin plot 
#comparing two distributions - ckd
kruskal.test(ckd_prs_standardised ~ type, data=prs)
pairwise.wilcox.test(prs$ckd_prs_standardised, prs$type, p.adjust.method = "BH")
kruskal.test(apol1 ~ pheno, data=prs)
#iga
kruskal.test(iga_prs_standardised ~ type, data=prs)
pairwise.wilcox.test(prs$iga_prs_standardised, prs$type, p.adjust.method = "BH")
kruskal.test(iga_prs_standardised ~ pheno, data=prs)
#ssns
kruskal.test(ssns_prs_standardised ~ type, data=prs)
pairwise.wilcox.test(prs$ssns_prs_standardised, prs$type, p.adjust.method = "BH")
kruskal.test(ssns_prs_standardised ~ pheno, data=prs)
#membranous
kruskal.test(membranous_prs_standardised ~ type, data=prs)
pairwise.wilcox.test(prs$membranous_prs, prs$type, p.adjust.method = "BH")
kruskal.test(membranous_prs_standardised ~ pheno, data=prs)
#pop <- aov(ldpred_beta_SUM ~ PHENO, data=subgroup.pheno.prs)
#summary(pop)
#t.test(x$ldpred_beta_SUM, y$ldpred_beta_SUM)

#plot with stats 
my_comparisons <- list(c("Case - High Rx APOL1", "Case - Non Rx APOL1 Control",
                         "Control - High Rx APOL1", "Control - Non Rx APOL1"))
ggplot(to_plot, aes(x=type, y=ckd_prs_standardised, fill=type)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  scale_fill_brewer(palette = "Dark2") +
  ylim(-4,5) + 
  theme_classic() +
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") +
  theme(axis.text.x = element_text(angle=30, vjust=0.6, hjust=0.5)) +
  stat_compare_means(comparisons = my_comparisons) + 
  #scale_x_discrete(labels = c("Control", "Solved Case")) + 
  stat_compare_means(label.y =50)


ggviolin(to_plot, x = "type", y = "ckd_prs_standardised", fill = "pheno",
         palette = c("#00AFBB", "#f2241f"),
         add = "boxplot", add.params = list(fill = "white"))+
  #facet_wrap(~type)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 8) + 
  stat_cor(p.accuracy = 0.01) +
  #scale_x_discrete(labels = c("Control", "EEHTN", "Primary HTN")) + 
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") 
#ggsave("slc34a3_subgroup_prs_stones_ggpubr_4groups.violin.png", height = 7, width = 10)
ggviolin(to_plot, x = "pheno", y = "ckd_prs_standardised", fill = "pheno",
         palette = c("#00AFBB", "#f2241f"),
         add = "boxplot", add.params = list(fill = "white"))+
  #facet_wrap(~type)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 8) + 
  stat_cor(p.accuracy = 0.01) +
  #scale_x_discrete(labels = c("Control", "EEHTN", "Primary HTN")) + 
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") 


#compare distibutions 

prs$prs_standardised <- as.numeric(prs$prs_standardised)
comparison_results <- prs %>%
  group_by(type) %>%
  kruskal.test(prs_standardised)
comparison_results 
kruskal.test(prs_standardised ~ pheno, data=pheno.prs)

# Create faceted density plots of PRS by type and pheno_esrf using ggplot2
ggplot(prs, aes(x = prs_standardised, fill = factor(type))) +
  geom_density(alpha = 0.5) +
  labs(title = "Faceted Density Plot of PRS by Type and pheno_esrf",
       x = "PRS",
       y = "Density") +
  scale_fill_discrete(name = "pheno_esrf") +
  facet_grid(type ~ ., scales = "free_y") 

#compare
library(dplyr)
library(stats)

# Assuming your_data_frame is your data frame
# Perform t-test or Wilcoxon rank-sum test for each combination of type and pheno_esrf
comparison_results <- prs %>%
  group_by(factor(type)) %>%
  summarize(p_value = wilcox.test(prs_standardised ~ factor(type))$p.value)

print(comparison_results)

library(dplyr)
library(stats)

# Assuming your_data_frame is your data frame
# Perform Kruskal-Wallis test for each type
comparison_results <- prs %>%
  group_by(type) %>%
  summarise(p_value = kruskal.test(prs_standardised ~ factor(pheno))$p.value)

print(comparison_results)

## Assuming your_data_frame is your data frame
# Perform one-way ANOVA
anova_result <- aov(ckd_prs_standardised ~ type, data = prs)
summary(anova_result)
# Perform post-hoc tests (Tukey's HSD)
posthoc_result <- TukeyHSD(anova_result)

# View ANOVA table
print(summary(anova_result))

# View post-hoc test results
print(posthoc_result)



#AUC 
library(DescTools)
library(boot)
Cstat(model)
#bootstrap confidence intervals 
FUN <- function(x,i) {
  r.glm <- glm(pheno ~ ., data=x[i,], family=binomial)
  Cstat(r.glm)
}
boot.res <- boot(formodel.pheno.prs, FUN, R=999) 
boot.ci(boot.res, type="perc")

# the percentile confidence intervals
boot.ci(boot.res, type="perc")


#calculate SD increase in PRS
sd_increase_in_prs <- exp(model$coefficients["prs_standardised"]*sd(prs$prs_standardised))
sd_increase_in_prs
#95% CI
five_percent <- exp(confint(model)["prs_standardised",1]*sd(prs$prs_standardised))
ninety_five_percent <- exp(confint(model)["prs_standardised",2]*sd(prs$prs_standardised))

#plot centile of model as prediction of getting USD
library(jtools)
library(OddsPlotty)
library(dplyr)
unsolved_plotting <- formodel.prs
unsolved_plotting$rank <- rank(unsolved_plotting$prs_standardised)
unsolved_plotting$rank <- round(unsolved_plotting$rank/length(unsolved_plotting$rank) * 100)
x <- unsolved_plotting %>% group_by(rank) %>% mutate(prevelance = (sum(pheno==1)/length(pheno)*100))
l <- x[,c(15,16)]
l <- unique(l)
l <- l[order(l$rank,decreasing = F),]
ggplot(l, aes(x=rank, y=prevelance)) + geom_smooth()
#fit model with rank 
unsolved_plotting_centile <- unsolved_plotting[,c(1,15,4:13)]
unsolved_centile <- glm(pheno~., family=binomial(link='logit'),data=unsolved_plotting_centile, maxit=10000)
effect_plot(data = unsolved_plotting_centile, unsolved_centile, pred=rank, interval = T, x.label = "PRS Centile",
            y.label = "Prevelance of USD")
#plot_summs(unsolved_model, model, coefs=coef_names)


#plot standardies graphs and distributions 
#create another pheno cohort to subdivide solved, slc34a3 and controls that have slc34a3 
subgroup.pheno.prs <- pheno.prs

subgroup.pheno.prs$PHENO <- as.factor(subgroup.pheno.prs$PHENO)
subgroup.pheno.prs$PHENO <- reorder(subgroup.pheno.prs$PHENO, subgroup.pheno.prs$prs_standardised, mean)

#stats on subgroup 
group_by(prs, type) %>%
  summarise(
    count = n(),
    mean = mean(prs_standardised, na.rm = TRUE),
    median = median(prs_standardised, na.rm = TRUE),
    sd = sd(prs_standardised, na.rm = TRUE)
  )



#plot with stats 
my_comparisons <- list(c("Case - High Rx APOL1", "Control - High Rx APOL1"), c("Case - Non Rx APOL1", "Control - Non Rx APOL1"), c("Case - High Rx APOL1","Case - Non Rx APOL1 Control"))

ggviolin(prs, x = "pheno", y = "membranous_prs_standardised", fill = "pheno",
         palette = c("#00AFBB", "#FC4E07","#f2241f", "black"),
         add = "boxplot", add.params = list(fill = "white"))+
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  #stat_compare_means(label.y = 8) + 
  #stat_cor(p.accuracy = 0.01) +
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") 
ggsave("slc34a3_subgroup_prs_stones_ggpubr_4groups.violin.png", height = 7, width = 10)


#power calculations 
n <- 25304
alpha <- 0.05
OR_interact <- 2
power <- pwr.f2.test(u=log(OR_interact), v1 =15,v2=25289, sig.level = alpha, power=NULL)

###all commerts 
prs_plotting <- prs
prs_plotting$rank <- rank(prs_plotting$prs_standardised)
prs_plotting$rank <- round(prs_plotting$rank/length(prs_plotting$rank) * 100)
x <- prs_plotting %>% group_by(rank) %>% mutate(prevelance = (sum(pheno==1)/length(pheno)*100))
l <- x[,c(18,19)]
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
df <- prs
df$prevalence <- ifelse(df$pheno == 1, 1, 0)
prevalence_by_cohort <- aggregate(prevalence ~ prs_standardised + type, df, mean)
prevalence_by_cohort$prevalence <- prevalence_by_cohort$prevalence * 100
# Create the plot - facet
ggplot(df, aes(x = prs_standardised, y = prevalence, color = type)) + 
  #geom_point() + 
  geom_smooth(method = "glm", se = FALSE, span = 0.5) +
  xlab("PRS") + 
  ylab("Prevalence") + 
  ggtitle("Prevalence of phenotype by PRS and cohort") +
  facet_wrap(~pheno, scales = "free_x")

#single plot
ggplot(df, aes(x = prs_standardised, y = prevalence, color = pheno, group = pheno)) + 
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

#plot
ggplot(dist, aes(x=prs_standardised, fill=pheno))+
  geom_density(alpha=0.2)+
  theme_classic()+
  labs(x="Polygenic Score", y="Density", fill="Cohort") + facet_grid(rows = vars(pheno)) +
  geom_vline(xintercept = 0, linetype = "dotted")


ggsave("densities_slc34a3_case_unsolved_cases_controls.density.png", height = 7, width = 7)

ggplot(dist, aes(x = prs_standardised, fill = pheno)) +
  geom_density(alpha = 0.2) +
  geom_vline(data = pheno.prs, 
             aes(xintercept = median_prs, group = pheno), 
             color = "black", 
             size = 1, linetype = "dotted") +
  theme_classic() +
  labs(x = "PRS", y = "Density", fill = "Cohort") +
  facet_grid(rows = vars(pheno))

ggsave("analysis/densities_slc34a3_case_unsolved_cases_controls.density_median_line.png", height = 7, width = 7)


#reshape
library(reshape2)
to_melt <- prs[,c(3,19:22)]
melted_df <- melt(to_melt, id.vars="pheno")
ggplot(melted_df, aes(x=variable, y=value, fill=factor(pheno))) +
  geom_violin(scale="width", draw_quantiles = c(0.25,0.5,0.75)) +
  geom_boxplot(width=0.1, position = position_dodge(0.9), alpha=0.5) +
  labs(x="Polygenic score type", y="PRS", fill = "pheno") + 
  theme_minimal() + 
  theme(legend.position = "top") +
  scale_fill_manual(values = c("0" = "blue", "1" = "red"))

to_melt <- prs[,c(3,19:22)]
melted_df <- melt(to_melt, id.vars="pheno")
ggplot(melted_df, aes(x=variable, y=value, color=factor(pheno))) +
  geom_line() +
  labs(x="Polygenic score type", y="PRS") + 
  theme_minimal() 

phenotype_0 <- melted_df[melted_df$pheno==0,]
phenotype_1 <- melted_df[melted_df$pheno==1,]


plot_0 <- ggplot(phenotype_0, aes(x=value, factor(variable))) +
  geom_density(alpha=0.5) +
  labs(x="PRS", y="Density", fill = "PRS") + 
  theme_minimal() 

ggplot(melted_df, aes(x=variable,y=value, group=variable))+
  geom_density_ridges2()
