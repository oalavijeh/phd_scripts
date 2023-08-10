#turn CKD grades into 0 or 1
library(data.table)
pkd_kaplan <- fread("phenotypes/final_phenotypes/all_cystic_tte_v16.ped")
pkd_kaplan <- pkd_kaplan[,c(2,3,4)]
names(pkd_kaplan)[names(pkd_kaplan)=="pheno_esrf"] <- "status"
names(pkd_kaplan)[names(pkd_kaplan)=="tte"] <- "time"
pkd_kaplan <- as.data.frame(pkd_kaplan)
pkd_kaplan$time <- as.numeric(pkd_kaplan$time)
#write.table(pkd_kaplan, "~/re_gecip/renal/oalavijeh/projects/phenotypes/stratified/pkd_kaplan_bioinformaticallyv16_svs_as_truncating.txt", col.names = T, row.names = F, quote = F, sep = '\t')
library(survival)
library(ggfortify)
#library(survminer, "~/re_gecip/renal/Rpackages")

fit <- survfit(Surv(time, status) ~ type, data = pkd_kaplan2)
autoplot(fit, conf.int = FALSE, censor = TRUE, main = "renal survival in 100K cystic renal disease cohort (n=419)", xlab ="Age at ESRF", ylab = "survival", surv.geom = 'step')
#no censor (crosses)
autoplot(fit, conf.int = FALSE, censor = FALSE, grid = FALSE,
         xlab ="Age at ESRF", ylab = "Percentage without ESRF", surv.geom = 'step') + theme_classic()




#see total graph using surviner 
ggsurvplot(fit, data = pkd_kaplan)

#customise graph in survminer 
pp = ggsurvplot(
  fit, 
  data = pkd_kaplan, 
  size = 1,                  # change line size
  risk.table = TRUE,         # Add risk table
  risk.table.col = "strata", # Risk table color by groups
  risk.table.height = 0.31,  # Useful to change when you have multiple groups
  ggtheme = theme_bw(),      # Change ggplot2 theme
  ylab = "Renal Survival",
  xlab= "Age/years",
  title = "100K cystic disease cohort", 
)
