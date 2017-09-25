### R code from vignette source 'msm-manual.Rnw'

###################################################
### code chunk number 1: msm-manual.Rnw:41-43
###################################################
version <- gsub("Version: +", "",
               packageDescription("msm", lib.loc=c("../..",.libPaths()))$Version)


###################################################
### code chunk number 2: msm-manual.Rnw:48-49
###################################################
cat(version)


###################################################
### code chunk number 3: msm-manual.Rnw:52-53
###################################################
cat(format(Sys.time(), "%d %B, %Y"))


###################################################
### code chunk number 4: msm-manual.Rnw:861-862
###################################################
options(width = 60)


###################################################
### code chunk number 5: msm-manual.Rnw:897-898
###################################################
library(msm)


###################################################
### code chunk number 6: msm-manual.Rnw:957-958
###################################################
cav[1:21,]


###################################################
### code chunk number 7: msm-manual.Rnw:968-969
###################################################
statetable.msm(state, PTNUM, data=cav)


###################################################
### code chunk number 8: msm-manual.Rnw:1020-1024
###################################################
Q  <-  rbind ( c(0, 0.25, 0, 0.25),
               c(0.166, 0, 0.166, 0.166),
               c(0, 0.25, 0, 0.25),
               c(0, 0, 0, 0) )


###################################################
### code chunk number 9: msm-manual.Rnw:1058-1060
###################################################
Q.crude  <- crudeinits.msm(state ~ years, PTNUM, data=cav,
                                   qmatrix=Q)


###################################################
### code chunk number 10: msm-manual.Rnw:1084-1086
###################################################
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                    qmatrix = Q, deathexact = 4)


###################################################
### code chunk number 11: msm-manual.Rnw:1114-1115 (eval = FALSE)
###################################################
## help(optim)


###################################################
### code chunk number 12: msm-manual.Rnw:1138-1139
###################################################
cav.msm


###################################################
### code chunk number 13: msm-manual.Rnw:1178-1180
###################################################
cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     qmatrix = Q, deathexact = 4, covariates = ~ sex)


###################################################
### code chunk number 14: msm-manual.Rnw:1188-1189
###################################################
cavsex.msm


###################################################
### code chunk number 15: msm-manual.Rnw:1205-1207
###################################################
qmatrix.msm(cavsex.msm, covariates=list(sex=0)) # Male
qmatrix.msm(cavsex.msm, covariates=list(sex=1)) # Female


###################################################
### code chunk number 16: msm-manual.Rnw:1219-1222 (eval = FALSE)
###################################################
## cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                   qmatrix = Q, deathexact = 4,
##                   covariates = list("1-2" = ~ sex, "1-4" = ~sex) )


###################################################
### code chunk number 17: msm-manual.Rnw:1233-1237 (eval = FALSE)
###################################################
## cav3.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                 qmatrix = Q, deathexact = 4,
##                 covariates = ~ sex,
##                 constraint = list(sex=c(1,2,3,1,2,3,2)) )


###################################################
### code chunk number 18: msm-manual.Rnw:1273-1277 (eval = FALSE)
###################################################
## cav4.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                 qmatrix = Q, deathexact = 4,
##                 control = list(trace=2, REPORT=1),
##                 fixedpars = c(6, 7) )


###################################################
### code chunk number 19: msm-manual.Rnw:1316-1317
###################################################
pmatrix.msm(cav.msm, t=10)


###################################################
### code chunk number 20: msm-manual.Rnw:1346-1347
###################################################
sojourn.msm(cav.msm)


###################################################
### code chunk number 21: msm-manual.Rnw:1359-1360
###################################################
pnext.msm(cav.msm)


###################################################
### code chunk number 22: msm-manual.Rnw:1385-1386
###################################################
totlos.msm(cav.msm)


###################################################
### code chunk number 23: msm-manual.Rnw:1408-1409
###################################################
qratio.msm(cav.msm, ind1=c(2,1), ind2=c(1,2))


###################################################
### code chunk number 24: msm-manual.Rnw:1417-1418
###################################################
hazard.msm(cavsex.msm)


###################################################
### code chunk number 25: msm-manual.Rnw:1427-1428 (eval = FALSE)
###################################################
## qmatrix.msm(cav.msm)


###################################################
### code chunk number 26: msm-manual.Rnw:1437-1438 (eval = FALSE)
###################################################
## qmatrix.msm(cavsex.msm, covariates = 0)


###################################################
### code chunk number 27: msm-manual.Rnw:1443-1444 (eval = FALSE)
###################################################
## qmatrix.msm(cavsex.msm, covariates = list(sex = 1))


###################################################
### code chunk number 28: msm-manual.Rnw:1470-1471
###################################################
plot(cav.msm, legend.pos=c(8, 1))


###################################################
### code chunk number 29: msm-manual.Rnw:1713-1715
###################################################
options(digits=3)
prevalence.msm(cav.msm, times=seq(0,20,2))


###################################################
### code chunk number 30: msm-manual.Rnw:1717-1718
###################################################
plot.prevalence.msm(cav.msm, mintime=0, maxtime=20)


###################################################
### code chunk number 31: msm-manual.Rnw:1852-1855
###################################################
options(digits=2)
pearson.msm(cav.msm, timegroups=2,
            transitions=c(1,2,3,4,5,6,7,8,9,9,9,10))


###################################################
### code chunk number 32: msm-manual.Rnw:1980-1992
###################################################
Qm <- rbind(c(0, 0.148, 0, 0.0171),
            c(0, 0, 0.202, 0.081),
            c(0, 0, 0, 0.126),
            c(0, 0, 0, 0))
ematrix <- rbind(c(0, 0.1, 0, 0),
                 c(0.1, 0, 0.1, 0),
                 c(0, 0.1, 0, 0),
                 c(0, 0, 0, 0))
cavmisc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                   qmatrix = Qm, ematrix = ematrix, deathexact = 4,
                   obstrue = firstobs)
cavmisc.msm


###################################################
### code chunk number 33: msm-manual.Rnw:2020-2024
###################################################
cavmiscsex.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                      qmatrix = Qm, ematrix = ematrix,
                      deathexact = 4, misccovariates = ~sex,
                      obstrue=firstobs)


###################################################
### code chunk number 34: msm-manual.Rnw:2026-2027
###################################################
cavmiscsex.msm


###################################################
### code chunk number 35: msm-manual.Rnw:2047-2049
###################################################
ematrix.msm(cavmiscsex.msm, covariates=list(sex=0))
ematrix.msm(cavmiscsex.msm, covariates=list(sex=1))


###################################################
### code chunk number 36: msm-manual.Rnw:2096-2098
###################################################
pearson.msm(cavmisc.msm, timegroups=2,
            transitions=c(1,2,3,4,5,6,7,8,9,9,9,10))


###################################################
### code chunk number 37: msm-manual.Rnw:2144-2146
###################################################
vit <- viterbi.msm(cavmisc.msm)
vit[vit$subject==100103,]


###################################################
### code chunk number 38: msm-manual.Rnw:2344-2345
###################################################
three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))


###################################################
### code chunk number 39: msm-manual.Rnw:2363-2375
###################################################
hmodel1 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18),
                hmmIdent(999))

fev1.msm <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q,
                deathexact=3, hmodel=hmodel1,
                hcovariates=list(~acute, ~acute, NULL),
                hcovinits = list(-8, -8, NULL),
                hconstraint = list(acute = c(1,1)))

fev1.msm

sojourn.msm(fev1.msm)


###################################################
### code chunk number 40: msm-manual.Rnw:2642-2643 (eval = FALSE)
###################################################
## help(msm)


###################################################
### code chunk number 41: msm-manual.Rnw:2651-2652 (eval = FALSE)
###################################################
## help.start()


