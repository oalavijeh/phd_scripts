library(genpwr)

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=27305, Case.Rate=0.04, k=NULL,
                  MAF=seq(0.001, 0.01, 0.001), OR=c(5),Alpha=0.0000000482,
                  True.Model=c("Additive"), 
                  Test.Model=c("Additive"))

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=3200, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.001), OR=c(11),Alpha=0.00000258,
                  True.Model=c("Dominant", "Additive"), 
                  Test.Model=c("Dominant", "Recessive", "Additive", "2df"))


pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=3200, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0005), OR=c(10),Alpha=0.00000258,
                  True.Model=c("Dominant", "Recessive", "Additive"), 
                  Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=3200, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0005), OR=c(5),Alpha=0.00000258,
                  True.Model=c("Dominant", "Recessive", "Additive"), 
                  Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=3200, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0005), OR=c(2),Alpha=0.00000258,
                  True.Model=c("Dominant", "Recessive", "Additive"), 
                  Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=3200, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0005), OR=c(5),Alpha=0.00000258,
                  True.Model="All", 
                  Test.Model=c("Dominant", "Recessive"))

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=3200, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0005), OR=c(5),Alpha=0.00000258,
                  True.Model="Dominant", 
                  Test.Model="Dominant")

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=3200, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0005), OR=c(5),Alpha=0.00000258,
                  True.Model="Recessive", 
                  Test.Model="Recessive")

ss <- genpwr.calc(calc = "ss", model = "logistic", ge.interaction = NULL,
                  OR=11, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.001), Power=0.8, Alpha=0.00000258,
                  True.Model=c("Dominant", "Additive"), 
                  Test.Model=c("Dominant","Additive"))

ss <- genpwr.calc(calc = "ss", model = "logistic", ge.interaction = NULL,
                  OR=5, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0001), Power=0.8, Alpha=0.00000258,
                  True.Model="Dominant", 
                  Test.Model="Dominant")

ss <- genpwr.calc(calc = "ss", model = "logistic", ge.interaction = NULL,
                  OR=5, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0001), Power=0.8, Alpha=0.00000258,
                  True.Model=c("Dominant", "Recessive", "Additive"), 
                  Test.Model=c("Dominant", "Additive"))

ss <- genpwr.calc(calc = "ss", model = "logistic", ge.interaction = NULL,
                  OR=20, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0001), Power=0.8, Alpha=0.00000258,
                  True.Model=c("Dominant", "Recessive", "Additive"), 
                  Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

ss <- genpwr.calc(calc = "ss", model = "logistic", ge.interaction = NULL,
                  OR=5, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.0001), Power=0.8, Alpha=0.00000258,
                  True.Model="Dominant", 
                  Test.Model="Dominant")

or <- genpwr.calc(calc = "es", model = "logistic", ge.interaction = NULL,
                  N=3200, Case.Rate=0.07, k=NULL,
                  MAF=seq(0.001, 0.01, 0.001), Power=0.8, Alpha=0.00000258,
                  True.Model="All", Test.Model="All")


power.plot(pw)
ss.plot(ss)
  