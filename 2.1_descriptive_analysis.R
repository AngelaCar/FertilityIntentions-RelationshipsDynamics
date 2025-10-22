# ---- Descriptive Analysis ----
library(mstate)
library(data.table)
library(survival)
library(Hmisc)
library(RColorBrewer)

# ---- Load data and correct for last observations ----
load("Data/data_mstate_ready_updated.Rda")
events(data)

tmat <- attr(data, "trans")
datatb <- as.data.table(data)

# correction for last transitions that happen at the same time as the entry into
# the state. This does not affect the results, but it will make sure we can
# estimate the confidence intervals
datatb <- datatb[!(Tstart_from_int1 == Tstop_from_int1)]
data <- subset(data, Tstart_from_int1 != Tstop_from_int1)

describe(data$cid)
describe(data$id)

# ---- Descriptives stats ----
data_unique <- datatb[, head(.SD, 1), by = "cid"]

# descriptive tables
data_unique[, range(agea)]
data_unique[, range(agep)]

data_unique[, mean(agea), by = start_state] # 1 = cohab, 2 = dating, 3 = married
data_unique[, sd(agea), by = start_state]
data_unique[, range(agea), by = start_state]
data_unique[, mean(reldur_base), by = start_state] # 1 = cohab, 2 = dating, 3 = married
data_unique[, sd(reldur_base), by = start_state]
data_unique[, range(reldur_base), by = start_state]
data_unique[, mean(relsatm), by = start_state]
data_unique[, sd(relsatm), by = start_state]
data_unique[, range(relsatm), by = start_state]
data_unique[, mean(relsatf), by = start_state]
data_unique[, sd(relsatf), by = start_state]
data_unique[, range(relsatf), by = start_state]

(table(data_unique[start_state == 1, educdiff]) / nrow(data_unique[start_state == 1])) * 100
(table(data_unique[start_state == 2, educdiff]) / nrow(data_unique[start_state == 2])) * 100
(table(data_unique[start_state == 3, educdiff]) / nrow(data_unique[start_state == 3])) * 100

(table(data_unique[start_state == 1, livewithpar]) / nrow(data_unique[start_state == 1])) * 100
(table(data_unique[start_state == 2, livewithpar]) / nrow(data_unique[start_state == 2])) * 100
(table(data_unique[start_state == 3, livewithpar]) / nrow(data_unique[start_state == 3])) * 100

# distribution of agreement by starting state
(table(data_unique[start_state == 1, agreement_intentions]) / nrow(data_unique[start_state == 1])) * 100
(table(data_unique[start_state == 2, agreement_intentions]) / nrow(data_unique[start_state == 2])) * 100
(table(data_unique[start_state == 3, agreement_intentions]) / nrow(data_unique[start_state == 3])) * 100
datatb[, table(agreement_intentions, trans)]
datatb[status == 1, table(agreement_intentions, trans)]

# with other agreement variable
# distribution of agreement by starting state
(table(data_unique[start_state == 1, agreement_timing]) / nrow(data_unique[start_state == 1])) * 100
(table(data_unique[start_state == 2, agreement_timing]) / nrow(data_unique[start_state == 2])) * 100
(table(data_unique[start_state == 3, agreement_timing]) / nrow(data_unique[start_state == 3])) * 100
datatb[, table(agreement_timing, trans)]
data_unique[, table(agreement_timing)]
data_unique[, table(agreement_timing) / nrow(data_unique)] * 100
datatb[status == 1, table(agreement_timing, trans)]

# class(data) <- "msdata"
attr(data, "trans") <- tmat
events(data)



# ---- Simple model with only fertility variable ----
# To reproduce Figure 1

cox_agree <- coxph(Surv(Tstart_from_int1, Tstop_from_int1, status) ~ 
                     agreement_intentions1.1 +
                     agreement_intentions1.2 + agreement_intentions1.3 + 
                     agreement_intentions1.4 + agreement_intentions1.5 +
                     agreement_intentions1.6 + agreement_intentions1.7 +
                     agreement_intentions1.8 + agreement_intentions1.9 +
                     agreement_intentions2.1 + agreement_intentions2.2 + 
                     agreement_intentions2.3 + agreement_intentions2.4 +
                     agreement_intentions2.5 + agreement_intentions2.6 + 
                     agreement_intentions2.7 + agreement_intentions2.8 +
                     agreement_intentions2.9 +
                     agreement_intentions3.1 + agreement_intentions3.2 + 
                     agreement_intentions3.3 + agreement_intentions3.4 +
                     agreement_intentions3.5 + agreement_intentions3.6 + 
                     agreement_intentions3.7 + agreement_intentions3.8 +
                     agreement_intentions3.9 +
                     strata(trans), 
                   data = data)

summary(cox_agree)

# Prediction from the model
newd_bothy <- data.frame(
  agreement_intentions = as.factor(rep("both_yes", 9)),
  trans = 1:9, strata = 1:9
)

newd_bothy$agreement_intentions <- factor(newd_bothy$agreement_intentions,
                                          levels = levels(data$agreement_intentions)
)
attr(newd_bothy, "trans") <- tmat
class(newd_bothy) <- c("msdata", "data.frame")
newd_bothy <- expand.covs(newd_bothy, "agreement_intentions", longnames = FALSE)

msf_bothyes <- msfit(cox_agree, newdata = newd_bothy, trans = tmat)
pt.bothyes <- probtrans(msf_bothyes, predt = 0)

# both no
newd_bothn <- data.frame(
  agreement_intentions = rep("both_no", 9),
  trans = 1:9, strata = 1:9
)
newd_bothn$agreement_intentions <- factor(newd_bothn$agreement_intentions, 
                                          levels = levels(data$agreement_intentions))
attr(newd_bothn, "trans") <- tmat
class(newd_bothn) <- c("msdata", "data.frame")
newd_bothn <- expand.covs(newd_bothn, "agreement_intentions", longnames = FALSE)
msf_bothno <- msfit(cox_agree, newdata = newd_bothn, trans = tmat)
pt.bothno <- probtrans(msf_bothno, predt = 0)

# He yes
newd_hey <- data.frame(
  agreement_intentions = rep("he_yes_she_no", 9),
  trans = 1:9, strata = 1:9
)
newd_hey$agreement_intentions <- factor(newd_hey$agreement_intentions,
                                        levels = levels(data$agreement_intentions)
)
attr(newd_hey, "trans") <- tmat
class(newd_hey) <- c("msdata", "data.frame")
newd_hey <- expand.covs(newd_hey, "agreement_intentions", longnames = FALSE)
msf_heyes <- msfit(cox_agree, newdata = newd_hey, trans = tmat)
pt.heyes <- probtrans(msf_heyes, predt = 0)

# She yes
newd_shey <- data.frame(
  agreement_intentions = rep("he_no_she_yes", 9),
  trans = 1:9, strata = 1:9
)
newd_shey$agreement_intentions <- factor(newd_shey$agreement_intentions, levels = levels(data$agreement_intentions))
attr(newd_shey, "trans") <- tmat
class(newd_shey) <- c("msdata", "data.frame")
newd_shey <- expand.covs(newd_shey, "agreement_intentions", longnames = FALSE)
msf_sheyes <- msfit(cox_agree, newdata = newd_shey, trans = tmat)
pt.sheyes <- probtrans(msf_sheyes, predt = 0)

# ---- Plots of prob trans - Figure 1 ----
col_trans <- (colorBlindness::Blue2DarkOrange12Steps)[c(1, 4, 7, 10, 12)]

mat_lay <- matrix(c(
  1, 1, 2, 2,
  3, 3, 4, 4
), nrow = 2, byrow = TRUE)

svg("Results/Figure1.svg", width = 10, height = 8, pointsize = 13)
{
  layout(mat_lay)
  par(
    mar = c(5.5, 5, 2, 3),
    oma = c(1, 2, 2, 2),
    cex.main = 1.9,
    cex.lab = 1.8,
    font.main = 1,
    mgp = c(2.5, 1, 0)
  )
  
  
  ggp_bothno <- plot(pt.bothno,
                     from = 1, label = "annotate",
                     xlab = "Time since fertility agreement report (months)",
                     main = "Both do not want a kid in next 2 years",
                     cex = 1.8,
                     cols = col_trans
  )
  
  ggp_bothyes <- plot(pt.bothyes,
                      from = 1, label = "annotate",
                      xlab = "Time since fertility agreement report (months)",
                      main = "Both want a kid in next 2 years",
                      cex = 1.8,
                      cols = col_trans
  )
  
  ggp_heyes <- plot(pt.heyes,
                    from = 1, label = "annotate",
                    xlab = "Time since fertility agreement report (months)",
                    main = "He wants a kid in next 2 years, she does not",
                    cex = 1.8,
                    cols = col_trans
  )
  ggp_sheyes <- plot(pt.sheyes,
                     from = 1, label = "annotate",
                     xlab = "Time since fertility agreement report (months)",
                     main = "She wants a kid in next 2 years, he does not",
                     cex = 1.8,
                     cols = col_trans
  )
}
dev.off()
