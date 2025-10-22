# ---- Multistate Model Selection ----
library(mstate)
library(data.table)
library(survival)
library(Hmisc)
library(RColorBrewer)

# ---- Load data and correct for last observations ----
load("Data/data_mstate_ready_updated.Rda")
tmat <- attr(data, "trans")
data <- subset(data, Tstart_from_int1 != Tstop_from_int1)
attr(data, "trans") <- tmat

# ---- Preliminary operations ----
### modify data a bit
# 1. remove transitions 2 and 4
data <- data[!(data$trans %in% c(2, 4)), ]

# 2. prepare for proportionality into state kid
data$trans2 <- ifelse(data$trans %in% c(7, 9), 7, data$trans)

#table(data$trans, data$trans2)
#table(data$trans2, data$status)

# make variable to distinguish between having moved in or married
data$kcohab <- ifelse(data$trans == 7, 1, 0)
data$kmarr <- ifelse(data$trans == 9, 1, 0)

# new transition specific covariates
data$agreement_intentions1.7.new <- data$agreement_intentions1.7 + data$agreement_intentions1.9
data$agreement_intentions2.7.new <- data$agreement_intentions2.7 + data$agreement_intentions2.9
data$agreement_intentions3.7.new <- data$agreement_intentions3.7 + data$agreement_intentions3.9

data$reldur_base_cat.7.new <- data$reldur_base_cat.7 + data$reldur_base_cat.9

data$educdiff1.7.new <- data$educdiff1.7 + data$educdiff1.9
data$educdiff2.7.new <- data$educdiff2.7 + data$educdiff2.9
data$educdiff3.7.new <- data$educdiff3.7 + data$educdiff3.9

data$agea.7.new <- data$agea.7 + data$agea.9

data$relsatf.7.new <- data$relsatf.7 + data$relsatf.9
data$relsatm.7.new <- data$relsatm.7 + data$relsatm.9

# 3. prepare for proportionality into state dissolution
data$trans3 <- ifelse(data$trans2 %in% c(6, 8), 6, data$trans2)

# make variable to distinguish between having moved in or married
data$dismarr <- ifelse(data$trans2 == 8, 1, 0)

# new transition specific covariates
data$agreement_intentions1.6.new <- data$agreement_intentions1.6 + data$agreement_intentions1.8
data$agreement_intentions2.6.new <- data$agreement_intentions2.6 + data$agreement_intentions2.8
data$agreement_intentions3.6.new <- data$agreement_intentions3.6 + data$agreement_intentions3.8

data$reldur_base_cat.6.new <- data$reldur_base_cat.6 + data$reldur_base_cat.8

data$educdiff1.6.new <- data$educdiff1.6 + data$educdiff1.8
data$educdiff2.6.new <- data$educdiff2.6 + data$educdiff2.8
data$educdiff3.6.new <- data$educdiff3.6 + data$educdiff3.8

data$agea.6.new <- data$agea.6 + data$agea.8

data$relsatf.6.new <- data$relsatf.6 + data$relsatf.8
data$relsatm.6.new <- data$relsatm.6 + data$relsatm.8

# 4: Recenter age at baseline, for transition-specific covariates
mean(data[data$agea.1 > 0, ]$agea.1)
mean(data[data$agea.3 > 0, ]$agea.3)
mean(data[data$agea.5 > 0, ]$agea.5)
mean(data[data$agea.6.new > 0, ]$agea.6.new)
mean(data[data$agea.7.new > 0, ]$agea.7.new)

data$agea_rep.1 <- ifelse(data$agea.1 == 0, 24.74341, data$agea.1)
data$agea_rep.3 <- ifelse(data$agea.3 == 0, 24.74341, data$agea.3)
data$agea_rep.5 <- ifelse(data$agea.5 == 0, 27.0226, data$agea.5)
data$agea_rep.6 <- ifelse(data$agea.6.new == 0, 27.9974, data$agea.6.new)
data$agea_rep.7 <- ifelse(data$agea.7.new == 0, 27.9974, data$agea.7.new)

# ---- Model selection ----
#* Model 1: ALL TRANSITION-SPECIFIC COVARIATES, ALL TRANSITIONS BUT 2,4
cox_full <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                    agreement_intentions1.1 + agreement_intentions1.3 +
                    agreement_intentions1.5 + agreement_intentions1.6 +
                    agreement_intentions1.7 + agreement_intentions1.8 +
                    agreement_intentions1.9 +
                    agreement_intentions2.1 + agreement_intentions2.3 +
                    agreement_intentions2.5 + agreement_intentions2.6 +
                    agreement_intentions2.7 + agreement_intentions2.8 +
                    agreement_intentions2.9 +
                    agreement_intentions3.1 + agreement_intentions3.3 +
                    agreement_intentions3.5 + agreement_intentions3.6 +
                    agreement_intentions3.7 + agreement_intentions3.8 +
                    agreement_intentions3.9 +
                    sex.1 + sex.3 + sex.5 + sex.6 + sex.7 +
                    sex.8 + sex.9 +
                    livewithpar.1 + livewithpar.3 + livewithpar.5 +
                    livewithpar.6 + livewithpar.7 + livewithpar.8 +
                    livewithpar.9 +
                    educdiff1.1 + educdiff1.3 + educdiff1.5 +
                    educdiff1.6 + educdiff1.7 + educdiff1.8 + educdiff1.9 +
                    educdiff2.1 + educdiff2.3 + educdiff2.5 +
                    educdiff2.6 + educdiff2.7 + educdiff2.8 + educdiff2.9 +
                    educdiff3.1 + educdiff3.3 + educdiff3.5 +
                    educdiff3.6 + educdiff3.7 + educdiff3.8 + educdiff3.9 +
                    agea.1 + agea.3 +
                    agea.5 + agea.6 + agea.7 + agea.8 +
                    agea.9 +
                    agegap1.1 + agegap1.3 + agegap1.5 +
                    agegap1.6 + agegap1.7 + agegap1.8 + agegap1.9 +
                    agegap2.1 + agegap2.3 + agegap2.5 +
                    agegap2.6 + agegap2.7 + agegap2.8 + agegap2.9 +
                    reldur_base.1 + reldur_base.3 +
                    reldur_base.5 + reldur_base.6 + reldur_base.7 + reldur_base.8 +
                    reldur_base.9 +
                    relsatm.1 + relsatm.3 + relsatm.5 +
                    relsatm.6 + relsatm.7 + relsatm.8 + relsatm.9 +
                    relsatf.1 + relsatf.3 + relsatf.5 +
                    relsatf.6 + relsatf.7 + relsatf.8 + relsatf.9 +
                    strata(trans), data = data)

s_coxfull <- summary(cox_full)
testfull <- cox.zph(cox_full)

#* Model 2: Interaction between agreement variable and duration of the relationship
#* categorical
cox_full_reldurint <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                              agreement_intentions1.1 * reldur_base_cat.1 +
                              agreement_intentions1.3 * reldur_base_cat.3 +
                              agreement_intentions1.5 * reldur_base_cat.5 +
                              agreement_intentions1.6 * reldur_base_cat.6 +
                              agreement_intentions1.7 * reldur_base_cat.7 +
                              agreement_intentions1.8 * reldur_base_cat.8 +
                              agreement_intentions1.9 * reldur_base_cat.9 +
                              agreement_intentions2.1 * reldur_base_cat.1 +
                              agreement_intentions2.3 * reldur_base_cat.3 +
                              agreement_intentions2.5 * reldur_base_cat.5 +
                              agreement_intentions2.6 * reldur_base_cat.6 +
                              agreement_intentions2.7 * reldur_base_cat.7 +
                              agreement_intentions2.8 * reldur_base_cat.8 +
                              agreement_intentions2.9 * reldur_base_cat.9 +
                              agreement_intentions3.1 * reldur_base_cat.1 +
                              agreement_intentions3.3 * reldur_base_cat.3 +
                              agreement_intentions3.5 * reldur_base_cat.5 +
                              agreement_intentions3.6 * reldur_base_cat.6 +
                              agreement_intentions3.7 * reldur_base_cat.7 +
                              agreement_intentions3.8 * reldur_base_cat.8 +
                              agreement_intentions3.9 * reldur_base_cat.9 +
                              sex.1 + sex.3 + sex.5 + sex.6 + sex.7 + sex.8 + sex.9 +
                              livewithpar.1 + livewithpar.3 + livewithpar.5 +
                              livewithpar.6 + livewithpar.7 + livewithpar.8 +
                              livewithpar.9 +
                              educdiff1.1 + educdiff1.3 + educdiff1.5 +
                              educdiff1.6 + educdiff1.7 + educdiff1.8 + educdiff1.9 +
                              educdiff2.1 + educdiff2.3 + educdiff2.5 +
                              educdiff2.6 + educdiff2.7 + educdiff2.8 + educdiff2.9 +
                              educdiff3.1 + educdiff3.3 + educdiff3.5 +
                              educdiff3.6 + educdiff3.7 + educdiff3.8 + educdiff3.9 +
                              agegap1.1 + agegap1.3 + agegap1.5 +
                              agegap1.6 + agegap1.7 + agegap1.8 + agegap1.9 +
                              agegap2.1 + agegap2.3 + agegap2.5 +
                              agegap2.6 + agegap2.7 + agegap2.8 + agegap2.9 +
                              relsatm.1 + relsatm.3 + relsatm.5 +
                              relsatm.6 + relsatm.7 + relsatm.8 + relsatm.9 +
                              relsatf.1 + relsatf.3 + relsatf.5 +
                              relsatf.6 + relsatf.7 + relsatf.8 + relsatf.9 +
                              agea.1 + agea.3 +
                              agea.5 + agea.6 + agea.7 +
                              agea.8 + agea.9 +
                              strata(trans), data = data)

# compare with previous model
AIC(cox_full)
AIC(cox_full_reldurint)
# Even though the AIC would prefer the simpler model, we have theoretical reasons 
# to keep the interaction in the model.

#* Model 3: remove livewithpar.3-livewithpar.9
cox_full_reldurint2 <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                               agreement_intentions1.1 * reldur_base_cat.1 +
                               agreement_intentions1.3 * reldur_base_cat.3 +
                               agreement_intentions1.5 * reldur_base_cat.5 +
                               agreement_intentions1.6 * reldur_base_cat.6 +
                               agreement_intentions1.7 * reldur_base_cat.7 +
                               agreement_intentions1.8 * reldur_base_cat.8 +
                               agreement_intentions1.9 * reldur_base_cat.9 +
                               agreement_intentions2.1 * reldur_base_cat.1 +
                               agreement_intentions2.3 * reldur_base_cat.3 +
                               agreement_intentions2.5 * reldur_base_cat.5 +
                               agreement_intentions2.6 * reldur_base_cat.6 +
                               agreement_intentions2.7 * reldur_base_cat.7 +
                               agreement_intentions2.8 * reldur_base_cat.8 +
                               agreement_intentions2.9 * reldur_base_cat.9 +
                               agreement_intentions3.1 * reldur_base_cat.1 +
                               agreement_intentions3.3 * reldur_base_cat.3 +
                               agreement_intentions3.5 * reldur_base_cat.5 +
                               agreement_intentions3.6 * reldur_base_cat.6 +
                               agreement_intentions3.7 * reldur_base_cat.7 +
                               agreement_intentions3.8 * reldur_base_cat.8 +
                               agreement_intentions3.9 * reldur_base_cat.9 +
                               sex.1 + sex.3 + sex.5 + sex.6 + sex.7 + sex.8 + sex.9 +
                               livewithpar.1 +
                               educdiff1.1 + educdiff1.3 + educdiff1.5 +
                               educdiff1.6 + educdiff1.7 + educdiff1.8 + educdiff1.9 +
                               educdiff2.1 + educdiff2.3 + educdiff2.5 +
                               educdiff2.6 + educdiff2.7 + educdiff2.8 + educdiff2.9 +
                               educdiff3.1 + educdiff3.3 + educdiff3.5 +
                               educdiff3.6 + educdiff3.7 + educdiff3.8 + educdiff3.9 +
                               agegap1.1 + agegap1.3 + agegap1.5 +
                               agegap1.6 + agegap1.7 + agegap1.8 + agegap1.9 +
                               agegap2.1 + agegap2.3 + agegap2.5 +
                               agegap2.6 + agegap2.7 + agegap2.8 + agegap2.9 +
                               relsatm.1 + relsatm.3 + relsatm.5 +
                               relsatm.6 + relsatm.7 + relsatm.8 + relsatm.9 +
                               relsatf.1 + relsatf.3 + relsatf.5 +
                               relsatf.6 + relsatf.7 + relsatf.8 + relsatf.9 +
                               agea.1 + agea.3 +
                               agea.5 + agea.6 + agea.7 +
                               agea.8 + agea.9 +
                               strata(trans), data = data)

s_coxfull_reldurint2 <- summary(cox_full_reldurint2)
anova(cox_full_reldurint, cox_full_reldurint2)

#* Model 4: Test sex as global variable
cox_full_testsex <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                            agreement_intentions1.1 * reldur_base_cat.1 +
                            agreement_intentions1.3 * reldur_base_cat.3 +
                            agreement_intentions1.5 * reldur_base_cat.5 +
                            agreement_intentions1.6 * reldur_base_cat.6 +
                            agreement_intentions1.7 * reldur_base_cat.7 +
                            agreement_intentions1.8 * reldur_base_cat.8 +
                            agreement_intentions1.9 * reldur_base_cat.9 +
                            agreement_intentions2.1 * reldur_base_cat.1 +
                            agreement_intentions2.3 * reldur_base_cat.3 +
                            agreement_intentions2.5 * reldur_base_cat.5 +
                            agreement_intentions2.6 * reldur_base_cat.6 +
                            agreement_intentions2.7 * reldur_base_cat.7 +
                            agreement_intentions2.8 * reldur_base_cat.8 +
                            agreement_intentions2.9 * reldur_base_cat.9 +
                            agreement_intentions3.1 * reldur_base_cat.1 +
                            agreement_intentions3.3 * reldur_base_cat.3 +
                            agreement_intentions3.5 * reldur_base_cat.5 +
                            agreement_intentions3.6 * reldur_base_cat.6 +
                            agreement_intentions3.7 * reldur_base_cat.7 +
                            agreement_intentions3.8 * reldur_base_cat.8 +
                            agreement_intentions3.9 * reldur_base_cat.9 +
                            sex +
                            livewithpar.1 +
                            educdiff1.1 + educdiff1.3 + educdiff1.5 +
                            educdiff1.6 + educdiff1.7 + educdiff1.8 + educdiff1.9 +
                            educdiff2.1 + educdiff2.3 + educdiff2.5 +
                            educdiff2.6 + educdiff2.7 + educdiff2.8 + educdiff2.9 +
                            educdiff3.1 + educdiff3.3 + educdiff3.5 +
                            educdiff3.6 + educdiff3.7 + educdiff3.8 + educdiff3.9 +
                            agegap1.1 + agegap1.3 + agegap1.5 +
                            agegap1.6 + agegap1.7 + agegap1.8 + agegap1.9 +
                            agegap2.1 + agegap2.3 + agegap2.5 +
                            agegap2.6 + agegap2.7 + agegap2.8 + agegap2.9 +
                            relsatm.1 + relsatm.3 + relsatm.5 +
                            relsatm.6 + relsatm.7 + relsatm.8 + relsatm.9 +
                            relsatf.1 + relsatf.3 + relsatf.5 +
                            relsatf.6 + relsatf.7 + relsatf.8 + relsatf.9 +
                            agea.1 + agea.3 +
                            agea.5 + agea.6 + agea.7 +
                            agea.8 + agea.9 +
                            strata(trans), data = data)
anova(cox_full_reldurint2, cox_full_testsex)
# keep sex as global

#* Model 5: test educational differences as global
cox_full_testeduc <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                             agreement_intentions1.1 * reldur_base_cat +
                             agreement_intentions1.3 * reldur_base_cat +
                             agreement_intentions1.5 * reldur_base_cat +
                             agreement_intentions1.6 * reldur_base_cat +
                             agreement_intentions1.7 * reldur_base_cat +
                             agreement_intentions1.8 * reldur_base_cat +
                             agreement_intentions1.9 * reldur_base_cat +
                             agreement_intentions2.1 * reldur_base_cat +
                             agreement_intentions2.3 * reldur_base_cat +
                             agreement_intentions2.5 * reldur_base_cat +
                             agreement_intentions2.6 * reldur_base_cat +
                             agreement_intentions2.7 * reldur_base_cat +
                             agreement_intentions2.8 * reldur_base_cat +
                             agreement_intentions2.9 * reldur_base_cat +
                             agreement_intentions3.1 * reldur_base_cat +
                             agreement_intentions3.3 * reldur_base_cat +
                             agreement_intentions3.5 * reldur_base_cat +
                             agreement_intentions3.6 * reldur_base_cat +
                             agreement_intentions3.7 * reldur_base_cat +
                             agreement_intentions3.8 * reldur_base_cat +
                             agreement_intentions3.9 * reldur_base_cat +
                             sex + livewithpar.1 +
                             educdiff +
                             agegap1.1 + agegap1.3 + agegap1.5 +
                             agegap1.6 + agegap1.7 + agegap1.8 + agegap1.9 +
                             agegap2.1 + agegap2.3 + agegap2.5 +
                             agegap2.6 + agegap2.7 + agegap2.8 + agegap2.9 +
                             relsatm.1 + relsatm.3 + relsatm.5 +
                             relsatm.6 + relsatm.7 + relsatm.8 + relsatm.9 +
                             relsatf.1 + relsatf.3 + relsatf.5 +
                             relsatf.6 + relsatf.7 + relsatf.8 + relsatf.9 +
                             agea.1 + agea.3 +
                             agea.5 + agea.6 + agea.7 +
                             agea.8 + agea.9 +
                             strata(trans), data = data)
anova(cox_full_testsex, cox_full_testeduc)
# we should keep the education differences as transition-specific

#* Model 6: test for age difference as global
cox_full_testagediff <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                                agreement_intentions1.1 * reldur_base_cat.1 +
                                agreement_intentions1.3 * reldur_base_cat.3 +
                                agreement_intentions1.5 * reldur_base_cat.5 +
                                agreement_intentions1.6 * reldur_base_cat.6 +
                                agreement_intentions1.7 * reldur_base_cat.7 +
                                agreement_intentions1.8 * reldur_base_cat.8 +
                                agreement_intentions1.9 * reldur_base_cat.9 +
                                agreement_intentions2.1 * reldur_base_cat.1 +
                                agreement_intentions2.3 * reldur_base_cat.3 +
                                agreement_intentions2.5 * reldur_base_cat.5 +
                                agreement_intentions2.6 * reldur_base_cat.6 +
                                agreement_intentions2.7 * reldur_base_cat.7 +
                                agreement_intentions2.8 * reldur_base_cat.8 +
                                agreement_intentions2.9 * reldur_base_cat.9 +
                                agreement_intentions3.1 * reldur_base_cat.1 +
                                agreement_intentions3.3 * reldur_base_cat.3 +
                                agreement_intentions3.5 * reldur_base_cat.5 +
                                agreement_intentions3.6 * reldur_base_cat.6 +
                                agreement_intentions3.7 * reldur_base_cat.7 +
                                agreement_intentions3.8 * reldur_base_cat.8 +
                                agreement_intentions3.9 * reldur_base_cat.9 +
                                sex +
                                livewithpar.1 +
                                educdiff1.1 + educdiff1.3 + educdiff1.5 +
                                educdiff1.6 + educdiff1.7 + educdiff1.8 + educdiff1.9 +
                                educdiff2.1 + educdiff2.3 + educdiff2.5 +
                                educdiff2.6 + educdiff2.7 + educdiff2.8 + educdiff2.9 +
                                educdiff3.1 + educdiff3.3 + educdiff3.5 +
                                educdiff3.6 + educdiff3.7 + educdiff3.8 + educdiff3.9 +
                                agegap +
                                relsatm.1 + relsatm.3 + relsatm.5 +
                                relsatm.6 + relsatm.7 + relsatm.8 + relsatm.9 +
                                relsatf.1 + relsatf.3 + relsatf.5 +
                                relsatf.6 + relsatf.7 + relsatf.8 + relsatf.9 +
                                agea.1 + agea.3 +
                                agea.5 + agea.6 + agea.7 +
                                agea.8 + agea.9 +
                                strata(trans), data = data)
anova(cox_full_testsex, cox_full_testagediff)
# marginally significant. Let's compare the AIC
AIC(cox_full_testsex)
AIC(cox_full_testagediff)
# AIC is slightly lower, I would keep the simpler model with global effect of age gap

#* Model 7: test for relationship satisfaction male partner
cox_full_testrelsatm <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                                agreement_intentions1.1 * reldur_base_cat.1 +
                                agreement_intentions1.3 * reldur_base_cat.3 +
                                agreement_intentions1.5 * reldur_base_cat.5 +
                                agreement_intentions1.6 * reldur_base_cat.6 +
                                agreement_intentions1.7 * reldur_base_cat.7 +
                                agreement_intentions1.8 * reldur_base_cat.8 +
                                agreement_intentions1.9 * reldur_base_cat.9 +
                                agreement_intentions2.1 * reldur_base_cat.1 +
                                agreement_intentions2.3 * reldur_base_cat.3 +
                                agreement_intentions2.5 * reldur_base_cat.5 +
                                agreement_intentions2.6 * reldur_base_cat.6 +
                                agreement_intentions2.7 * reldur_base_cat.7 +
                                agreement_intentions2.8 * reldur_base_cat.8 +
                                agreement_intentions2.9 * reldur_base_cat.9 +
                                agreement_intentions3.1 * reldur_base_cat.1 +
                                agreement_intentions3.3 * reldur_base_cat.3 +
                                agreement_intentions3.5 * reldur_base_cat.5 +
                                agreement_intentions3.6 * reldur_base_cat.6 +
                                agreement_intentions3.7 * reldur_base_cat.7 +
                                agreement_intentions3.8 * reldur_base_cat.8 +
                                agreement_intentions3.9 * reldur_base_cat.9 +
                                sex +
                                livewithpar.1 +
                                educdiff1.1 + educdiff1.3 + educdiff1.5 +
                                educdiff1.6 + educdiff1.7 + educdiff1.8 + educdiff1.9 +
                                educdiff2.1 + educdiff2.3 + educdiff2.5 +
                                educdiff2.6 + educdiff2.7 + educdiff2.8 + educdiff2.9 +
                                educdiff3.1 + educdiff3.3 + educdiff3.5 +
                                educdiff3.6 + educdiff3.7 + educdiff3.8 + educdiff3.9 +
                                agegap +
                                relsatm +
                                relsatf.1 + relsatf.3 + relsatf.5 +
                                relsatf.6 + relsatf.7 + relsatf.8 + relsatf.9 +
                                agea.1 + agea.3 +
                                agea.5 + agea.6 + agea.7 +
                                agea.8 + agea.9 +
                                strata(trans), data = data)
anova(cox_full_testagediff, cox_full_testrelsatm)
# should have transition specific effect

#* Model 8: test for relationship satisfaction female partner
cox_full_testrelsatf <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                                agreement_intentions1.1 * reldur_base_cat +
                                agreement_intentions1.3 * reldur_base_cat +
                                agreement_intentions1.5 * reldur_base_cat +
                                agreement_intentions1.6 * reldur_base_cat +
                                agreement_intentions1.7 * reldur_base_cat +
                                agreement_intentions1.8 * reldur_base_cat +
                                agreement_intentions1.9 * reldur_base_cat +
                                agreement_intentions2.1 * reldur_base_cat +
                                agreement_intentions2.3 * reldur_base_cat +
                                agreement_intentions2.5 * reldur_base_cat +
                                agreement_intentions2.6 * reldur_base_cat +
                                agreement_intentions2.7 * reldur_base_cat +
                                agreement_intentions2.8 * reldur_base_cat +
                                agreement_intentions2.9 * reldur_base_cat +
                                agreement_intentions3.1 * reldur_base_cat +
                                agreement_intentions3.3 * reldur_base_cat +
                                agreement_intentions3.5 * reldur_base_cat +
                                agreement_intentions3.6 * reldur_base_cat +
                                agreement_intentions3.7 * reldur_base_cat +
                                agreement_intentions3.8 * reldur_base_cat +
                                agreement_intentions3.9 * reldur_base_cat +
                                sex + livewithpar.1 +
                                educdiff1.1 + educdiff1.3 + educdiff1.5 +
                                educdiff1.6 + educdiff1.7 + educdiff1.8 + educdiff1.9 +
                                educdiff2.1 + educdiff2.3 + educdiff2.5 +
                                educdiff2.6 + educdiff2.7 + educdiff2.8 + educdiff2.9 +
                                educdiff3.1 + educdiff3.3 + educdiff3.5 +
                                educdiff3.6 + educdiff3.7 + educdiff3.8 + educdiff3.9 +
                                agegap +
                                relsatm.1 + relsatm.3 + relsatm.5 +
                                relsatm.6 + relsatm.7 + relsatm.8 + relsatm.9 +
                                relsatf +
                                agea.1 + agea.3 +
                                agea.5 + agea.6 + agea.7 +
                                agea.8 + agea.9 +
                                strata(trans), data = data)
anova(cox_full_testagediff, cox_full_testrelsatf)
# should also stay separately

#* Model 9: Test for proportionality into state kid from cohabitation/marriage
cox_full_propkid <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                            agreement_intentions1.1 * reldur_base_cat.1 +
                            agreement_intentions1.3 * reldur_base_cat.3 +
                            agreement_intentions1.5 * reldur_base_cat.5 +
                            agreement_intentions1.6 * reldur_base_cat.6 +
                            agreement_intentions1.7.new * reldur_base_cat.7.new +
                            agreement_intentions1.8 * reldur_base_cat.8 +
                            agreement_intentions2.1 * reldur_base_cat.1 +
                            agreement_intentions2.3 * reldur_base_cat.3 +
                            agreement_intentions2.5 * reldur_base_cat.5 +
                            agreement_intentions2.6 * reldur_base_cat.6 +
                            agreement_intentions2.7.new * reldur_base_cat.7.new +
                            agreement_intentions2.8 * reldur_base_cat.8 +
                            agreement_intentions3.1 * reldur_base_cat.1 +
                            agreement_intentions3.3 * reldur_base_cat.3 +
                            agreement_intentions3.5 * reldur_base_cat.5 +
                            agreement_intentions3.6 * reldur_base_cat.6 +
                            agreement_intentions3.7.new * reldur_base_cat.7.new +
                            agreement_intentions3.8 * reldur_base_cat.8 +
                            sex +
                            livewithpar.1 +
                            educdiff1.1 + educdiff1.3 + educdiff1.5 +
                            educdiff1.6 + educdiff1.7.new + educdiff1.8 +
                            educdiff2.1 + educdiff2.3 + educdiff2.5 +
                            educdiff2.6 + educdiff2.7.new + educdiff2.8 +
                            educdiff3.1 + educdiff3.3 + educdiff3.5 +
                            educdiff3.6 + educdiff3.7.new + educdiff3.8 +
                            agea.1 + agea.3 + agea.5 + agea.6 + agea.7.new + agea.8 +
                            agegap +
                            relsatm.1 + relsatm.3 + relsatm.5 +
                            relsatm.6 + relsatm.7.new + relsatm.8 +
                            relsatf.1 + relsatf.3 + relsatf.5 +
                            relsatf.6 + relsatf.7.new + relsatf.8 +
                            kmarr +
                            strata(trans2), data = data)
testfull_propkid <- cox.zph(cox_full_propkid)
# we look at the results for the variable kmarr -> they indicate that proportionality
# is satisfied so we can keep this simpler model

#* Model 10: Proportionality of transitions to state dissolution from cohabitation/marriage
cox_full_propdiss <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                             agreement_intentions1.1 * reldur_base_cat.1 +
                             agreement_intentions1.3 * reldur_base_cat.3 +
                             agreement_intentions1.5 * reldur_base_cat.5 +
                             agreement_intentions1.6.new * reldur_base_cat.6.new +
                             agreement_intentions1.7.new * reldur_base_cat.7.new +
                             agreement_intentions2.1 * reldur_base_cat.1 +
                             agreement_intentions2.3 * reldur_base_cat.3 +
                             agreement_intentions2.5 * reldur_base_cat.5 +
                             agreement_intentions2.6.new * reldur_base_cat.6.new +
                             agreement_intentions2.7.new * reldur_base_cat.7.new +
                             agreement_intentions3.1 * reldur_base_cat.1 +
                             agreement_intentions3.3 * reldur_base_cat.3 +
                             agreement_intentions3.5 * reldur_base_cat.5 +
                             agreement_intentions3.6.new * reldur_base_cat.6.new +
                             agreement_intentions3.7.new * reldur_base_cat.7.new +
                             sex +
                             livewithpar.1 +
                             educdiff1.1 + educdiff1.3 + educdiff1.5 +
                             educdiff1.6.new + educdiff1.7.new +
                             educdiff2.1 + educdiff2.3 + educdiff2.5 +
                             educdiff2.6.new + educdiff2.7.new +
                             educdiff3.1 + educdiff3.3 + educdiff3.5 +
                             educdiff3.6.new + educdiff3.7.new +
                             agea.1 + agea.3 + agea.5 + agea.6.new + agea.7.new +
                             agegap +
                             relsatm.1 + relsatm.3 + relsatm.5 +
                             relsatm.6.new + relsatm.7.new +
                             relsatf.1 + relsatf.3 + relsatf.5 +
                             relsatf.6.new + relsatf.7.new +
                             kmarr + dismarr +
                             strata(trans3), data = data)
testfull_propdiss <- cox.zph(cox_full_propdiss)
# we look at the results for the dismarr variable -> we can keep the proportional model

#* Model 11: we specify age at baseline effects smooth
cox_agesmooth <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~
                         agreement_intentions1.1 * reldur_base_cat.1 +
                         agreement_intentions1.3 * reldur_base_cat.3 +
                         agreement_intentions1.5 * reldur_base_cat.5 +
                         agreement_intentions1.6.new * reldur_base_cat.6.new +
                         agreement_intentions1.7.new * reldur_base_cat.7.new +
                         agreement_intentions2.1 * reldur_base_cat.1 +
                         agreement_intentions2.3 * reldur_base_cat.3 +
                         agreement_intentions2.5 * reldur_base_cat.5 +
                         agreement_intentions2.6.new * reldur_base_cat.6.new +
                         agreement_intentions2.7.new * reldur_base_cat.7.new +
                         agreement_intentions3.1 * reldur_base_cat.1 +
                         agreement_intentions3.3 * reldur_base_cat.3 +
                         agreement_intentions3.5 * reldur_base_cat.5 +
                         agreement_intentions3.6.new * reldur_base_cat.6.new +
                         agreement_intentions3.7.new * reldur_base_cat.7.new +
                         sex +
                         livewithpar.1 +
                         educdiff1.1 + educdiff1.3 + educdiff1.5 +
                         educdiff1.6.new + educdiff1.7.new +
                         educdiff2.1 + educdiff2.3 + educdiff2.5 +
                         educdiff2.6.new + educdiff2.7.new +
                         educdiff3.1 + educdiff3.3 + educdiff3.5 +
                         educdiff3.6.new + educdiff3.7.new +
                         pspline(agea_rep.1, df = 4) +
                         pspline(agea_rep.3, df = 4) +
                         pspline(agea_rep.5, df = 4) +
                         pspline(agea_rep.6, df = 4) +
                         pspline(agea_rep.7, df = 4) +
                         agegap +
                         relsatm.1 + relsatm.3 + relsatm.5 +
                         relsatm.6.new + relsatm.7.new +
                         relsatf.1 + relsatf.3 + relsatf.5 +
                         relsatf.6.new + relsatf.7.new +
                         kmarr + dismarr +
                         strata(trans3), data = data)
AIC(cox_agesmooth)
AIC(cox_full_propdiss)

# We have theoretical reasons to want to estimate smooth age effects. Additionally,
# age was giving some indication of non-proportionality

# we select this as the best fitting model
bestmod <- cox_agesmooth
save(bestmod, file = "Results/bestmodel_updated.Rda")
save(data, file = "Results/datafitting_updated.Rda")
