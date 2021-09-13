
library(Epi)
library(survival)
library(haven)
library(foreign)
library(broom)
library(popEpi)
library(dplyr)

setwd("R:/LRWE_Proj26/shared/SL/dataset/analysis") 
db <-read_dta("db.dta")
nrow(db)
attach(db)


##############################################################################################################################################

#####CRM-mortality####

##Diab & Male
db11 <-Lexis(entry = list(period  = yearin,
                          age  = agein),
             exit = list(period  = outm),
             exit.status = crm,
             id = patid,
             data = subset(db, DM == 1 & gender ==1))

dbs11  <-splitMulti(db11, age = seq(30,100,1), period= seq(1998,2019,1))

a.kn <- with(subset(dbs11, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
p.kn <- with(subset(dbs11, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))

r11 <- glm((lex.Xst==1)~Ns(age, knots = a.kn)*Ns(period, knots = p.kn),
           family = poisson,
           offset = log(lex.dur),
           data   = dbs11)

##Diab & Female
db12 <-Lexis(entry = list(period  = yearin,
                          age  = agein),
             exit = list(period  = outm),
             exit.status = crm,
             id = patid,
             data = subset(db, DM == 1 & gender ==2))

dbs12  <-splitMulti(db12, age = seq(30,100,1), period= seq(1998,2019,1))

a.kn <- with(subset(dbs12, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
p.kn <- with(subset(dbs12, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))

r12 <- glm((lex.Xst==1)~Ns(age, knots = a.kn)*Ns(period, knots = p.kn),
           family = poisson,
           offset = log(lex.dur),
           data   = dbs12)



##No-Diab & Male
db01 <-Lexis(entry = list(period = yearin,
                          age = agein),
             exit = list(period = outm),
             exit.status = crm,
             id = patid,
             data = subset(db, DM == 0  & gender ==1))

dbs01  <-splitMulti(db01, age = seq(30,100,1), period= seq(1998,2019,1))

a.kn <- with(subset(dbs01, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
p.kn <- with(subset(dbs01, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))

r01 <- glm((lex.Xst==1)~Ns(age, knots = a.kn)*Ns(period, knots = p.kn),
           family = poisson,
           offset = log(lex.dur),
           data   = dbs01)

##No-Diab & Female
db02 <-Lexis(entry = list(period = yearin,
                          age = agein),
             exit = list(period = outm),
             exit.status = crm,
             id = patid,
             data = subset(db, DM == 0  & gender ==2))

dbs02  <-splitMulti(db02, age = seq(30,100,1), period= seq(1998,2019,1))

a.kn <- with(subset(dbs02, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
p.kn <- with(subset(dbs02, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))

r02 <- glm((lex.Xst==1)~Ns(age, knots = a.kn)*Ns(period, knots = p.kn),
           family = poisson,
           offset = log(lex.dur),
           data   = dbs02)

##Combined DM+NonDM
age          <- c(35:100)
period       <- seq(1998,2019,0.1)
nd           <- expand.grid(age, period)
colnames(nd) <- c("age","period")
nd           <- cbind(nd, lex.dur=1000)

p11           <- ci.pred(r11, newdata = nd, Exp = T)
colnames(p11) <- c("rate", "lb", "ub")
p11           <- cbind(p11, nd, gender = 1, DM = 1, out="crm")

p12           <- ci.pred(r12, newdata = nd, Exp = T)
colnames(p12)  <- c("rate", "lb", "ub")
p12           <- cbind(p12, nd, gender = 2, DM = 1, out="crm")

p01           <- ci.pred(r01, newdata = nd, Exp = T)
colnames(p01) <- c("rate", "lb", "ub")
p01           <- cbind(p01, nd, gender = 1, DM =0, out="crm")

p02           <- ci.pred(r02, newdata = nd, Exp = T)
colnames(p02) <- c("rate", "lb", "ub")
p02           <- cbind(p02, nd, gender=2, DM = 0, out="crm")
res_crm      <- rbind(p11, p12, p01, p02)

###############################################################################################

#####ALL-cause-mortality####

##Diab & Male
db11 <-Lexis(entry = list(period  = yearin,
                          age  = agein),
             exit = list(period  = outm),
             exit.status = acm,
             id = patid,
             data = subset(db, DM == 1 & gender ==1))

dbs11  <-splitMulti(db11, age = seq(30,100,1), period= seq(1998,2019,1))

a.kn <- with(subset(dbs11, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
p.kn <- with(subset(dbs11, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))

r11 <- glm((lex.Xst==1)~Ns(age, knots = a.kn)*Ns(period, knots = p.kn),
           family = poisson,
           offset = log(lex.dur),
           data   = dbs11)

##Diab & Female
db12 <-Lexis(entry = list(period  = yearin,
                          age  = agein),
             exit = list(period  = outm),
             exit.status = acm,
             id = patid,
             data = subset(db, DM == 1 & gender ==2))

dbs12  <-splitMulti(db12, age = seq(30,100,1), period= seq(1998,2019,1))

a.kn <- with(subset(dbs12, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
p.kn <- with(subset(dbs12, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))

r12 <- glm((lex.Xst==1)~Ns(age, knots = a.kn)*Ns(period, knots = p.kn),
           family = poisson,
           offset = log(lex.dur),
           data   = dbs12)



##No-Diab & Male
db01 <-Lexis(entry = list(period = yearin,
                          age = agein),
             exit = list(period = outm),
             exit.status = acm,
             id = patid,
             data = subset(db, DM == 0  & gender ==1))

dbs01  <-splitMulti(db01, age = seq(30,100,1), period= seq(1998,2019,1))

a.kn <- with(subset(dbs01, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
p.kn <- with(subset(dbs01, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))

r01 <- glm((lex.Xst==1)~Ns(age, knots = a.kn)*Ns(period, knots = p.kn),
           family = poisson,
           offset = log(lex.dur),
           data   = dbs01)

##No-Diab & Female
db02 <-Lexis(entry = list(period = yearin,
                          age = agein),
             exit = list(period = outm),
             exit.status = acm,
             id = patid,
             data = subset(db, DM == 0  & gender ==2))

dbs02  <-splitMulti(db02, age = seq(30,100,1), period= seq(1998,2019,1))

a.kn <- with(subset(dbs02, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
p.kn <- with(subset(dbs02, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))

r02 <- glm((lex.Xst==1)~Ns(age, knots = a.kn)*Ns(period, knots = p.kn),
           family = poisson,
           offset = log(lex.dur),
           data   = dbs02)

##Combined DM+NonDM
age          <- c(35:100)
period       <- seq(1998,2019,0.1)
nd           <- expand.grid(age, period)
colnames(nd) <- c("age","period")
nd           <- cbind(nd, lex.dur=1000)

p11           <- ci.pred(r11, newdata = nd, Exp = T)
colnames(p11) <- c("rate", "lb", "ub")
p11           <- cbind(p11, nd, gender = 1, DM = 1, out="acm")

p12           <- ci.pred(r12, newdata = nd, Exp = T)
colnames(p12)  <- c("rate", "lb", "ub")
p12           <- cbind(p12, nd, gender = 2, DM = 1, out="acm")

p01           <- ci.pred(r01, newdata = nd, Exp = T)
colnames(p01) <- c("rate", "lb", "ub")
p01           <- cbind(p01, nd, gender = 1, DM =0, out="acm")

p02           <- ci.pred(r02, newdata = nd, Exp = T)
colnames(p02) <- c("rate", "lb", "ub")
p02           <- cbind(p02, nd, gender=2, DM = 0, out="acm")
res_acm      <- rbind(p11, p12, p01, p02)
res          <-rbind(res_crm, res_acm)
write.dta(res, file = "res.dta")






