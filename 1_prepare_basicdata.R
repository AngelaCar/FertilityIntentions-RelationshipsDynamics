# ----- Fertility intentions agreement and partnership's dynamics ----- 
## Data preparation
## Two datasets have been prepared externally in Stata
##      1. - all waves all anchors
##      2. - all waves all partners
## Here: merge the two data, retain only covariates that we are going to use 
##       for models, identify wave when both partners answer fertility intentions 
##       questions, code intentions and agreement, restrict on age, restrict on
##       presence of children at the time of fertility questions

# ---- libraries ----
library(readstata13)
library(data.table)
library(Hmisc)

# ---- data ----
anchors <- read.dta13(file = "Data/allwaves_anchor.dta",
                      nonint.factors = TRUE,
                      generate.factors = TRUE)
partners <- read.dta13(file = "Data/allwaves_partner.dta", 
                       nonint.factors = TRUE,
                       generate.factors = TRUE)

levels(anchors$wave)
levels(partners$wave)
anchors$wave <- factor(anchors$wave, 
                        levels = c("1 2008/09", "2 2009/10", "3 2010/11" , 
                                   "4 2011/12", "5 2012/13", "6 2013/14",
                                   "7 2014/15", "8 2015/16", "9 2016/17",
                                   "10 2017/18", "11 2018/19", "12 2019/20",
                                   "13 2020/21", "14 2021/22"), 
                        labels = c("2008/09", "2009/10", "2010/11" , 
                                   "2011/12", "2012/13", "2013/14",
                                   "2014/15", "2015/16", "2016/17",
                                   "2017/18", "2018/19", "2019/20",
                                   "2020/21", "2021/22"))

partners$wave <- factor(partners$wave, 
                           levels = c("1 2008/09", "2", "3", "4", "5", "6", "7",
                                      "8", "9", "10", "11", "12", "13", "14"), 
                           labels = c("2008/09", "2009/10", "2010/11" , 
                                      "2011/12", "2012/13", "2013/14",
                                      "2014/15", "2015/16", "2016/17",
                                      "2017/18", "2018/19", "2019/20",
                                      "2020/21", "2021/22"))

# ---- merge anchor with partners ----
bothpartners <- merge(anchors, partners,
                      by = c("id", "pid", "wave"))
# couple id
bothpartners$coupleid <- paste(bothpartners$id, bothpartners$pid, sep = "-")

bothpartners <- as.data.table(bothpartners)
setorder(bothpartners, coupleid, wave)

bothpartners[, couplesharedint := 1:.N, by = coupleid]

bothpartners[couplesharedint == 1, cid := 1:.N]
bothpartners[, cid := zoo::na.locf(cid), by = coupleid]

# ---- Fertility Intentions Questions ----
# ---- Anchor's fertility intentions ----
# frt6 (W 1-2): When you think realistically about having [additional] children
#               how many [more] children do you think you will have?
bothpartners[, table(frt6)]
bothpartners[, frt6 := factor(frt6, 
                              levels = c("-5 Inconsistent value",
                                         "-4 Filter error / Incorrect entry",
                                         "-3 Does not apply", 
                                         "-2 No answer",
                                         "1 One (additional) child",
                                         "2 Two (additional) children",
                                         "3 Three (additional) children", 
                                         "4 Four or more (additional) children", 
                                         "5 I'm not sure",
                                         "6 I haven't thought about that",
                                         "7 No (additional) children"),
                              labels = c(-5:-2, 1:7))]
bothpartners[, table(frt6)]

# frt26 (W 3-13): When you think realistically about having children: how many 
#                 biological or adoptive children do you think you will have?
# Here -3 Does not apply indicate that anchor has already children or they are 
# expecting. 
bothpartners[, table(frt26)]
bothpartners[, frt26 := factor(frt26, 
                              levels = c("-5 Inconsistent value",
                                         "-4 Filter error / Incorrect entry",
                                         "-3 Does not apply", 
                                         "-2 No answer",
                                         "1 One child",
                                         "2 Two children",
                                         "3 Three children", 
                                         "4 Four or more children", 
                                         "5 I'm not sure",
                                         "6 I haven't thought about that",
                                         "7 No children"),
                              labels = c(-5:-2, 1:7))]
bothpartners[, table(frt26)]

# Combine to get the fertility intentions of the anchor
bothpartners[is.na(frt26), frt6_combined := frt6]
bothpartners[!is.na(frt26), frt6_combined := frt26]
View(bothpartners[, c("id", "pid", "cid", "wave", "frt6", "frt26", "frt7",
                      "pfrt6", "pfrt7", "frt6_combined")])
bothpartners[, table(frt6_combined)]

# frt7 (W 1-13): Do you intend to become a mother/father in the next two years (W1)
#                Do you intend to have [another] child within the next two years?
# Here -3 is because this question is asked only if people indicate positive intentions
# in previous question. This question is also not asked to individuals who are infertily or
# currently pregnant, nor to those whose partner is infertily or currently pregnant
# From W 11 only asked to individuals who are at least 17.

bothpartners[, table(frt7)]
bothpartners[, frt7 := factor(frt7, 
                               levels = c("-5 Inconsistent value",
                                          "-4 Filter error / Incorrect entry",
                                          "-3 Does not apply", 
                                          "-2 No answer",
                                          "-1 Don't know", "1 Yes, definitely",
                                          "2 Yes, perhaps", "3 No, probably not",
                                          "4 No, definitely not",
                                          "7 I haven't thought about that yet"),
                               labels = c(-5:-1, 1:4, 7))]
bothpartners[, table(frt7)]

bothpartners[, table(frt6_combined, frt7)]

# answer they want no children: frt6_combined == 7
bothpartners[, intentions_anchor := ifelse(frt6_combined == 7, "def_no", NA)]

# answer they want one or more children: frt6_combined 1-4
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(1:4) & frt7 == 1,
                                           "def_yes", intentions_anchor)]
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(1:4) & frt7 == 2,
                                           "prob_yes", intentions_anchor)]
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(1:4) & frt7 == 3,
                                           "prob_no", intentions_anchor)]
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(1:4) & frt7 == 4,
                                           "def_no", intentions_anchor)]
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(1:4) & frt7 == 7,
                                           "prob_no", intentions_anchor)]

# answer they are not sure about expected number of children: frt6_combined 5-6
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(5:6) & frt7 == 1,
                                           "def_yes", intentions_anchor)]
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(5:6) & frt7 == 2,
                                           "prob_yes", intentions_anchor)]
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(5:6) & frt7 == 3,
                                           "prob_no", intentions_anchor)]
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(5:6) & frt7 == 4,
                                           "def_no", intentions_anchor)]
bothpartners[, intentions_anchor := ifelse(frt6_combined %in% c(5:4) & frt7 == 7,
                                           "prob_no", intentions_anchor)]

bothpartners[, describe(intentions_anchor)]
View(bothpartners[is.na(intentions_anchor), c("id", "frt6", "frt26", "frt6_combined", "frt7")])
bothpartners[is.na(intentions_anchor), table(frt6_combined, frt7)]

# fertility intentions binary
bothpartners[!is.na(intentions_anchor), 
             intentions_anchor_binary := ifelse(intentions_anchor %in% c("def_no", "prob_no"),
                                                "no", "yes")]
bothpartners[, describe(intentions_anchor_binary)]

# ---- Partner's fertility intentions ----
# As in frt6 for anchor
bothpartners[, table(pfrt6)]
bothpartners[, pfrt6 := factor(pfrt6, 
                              levels = c("-9 Invalid multiple answer",
                                         "-4 Filter error / Incorrect entry",
                                         "-3 Does not apply", 
                                         "-2 No answer",
                                         "-1 I'm not sure",
                                         "1 One (additonal) child",
                                         "2 Two (additional) children",
                                         "3 Three (additional) children", 
                                         "4 Four or more (additional) children", 
                                         "5 I'm not sure",
                                         "6 I haven't thought about that",
                                         "7 No (addtional) children"),
                              labels = c(-9, -4:-1, 1:7))]
bothpartners[, table(pfrt6)]

# As in frt26, except only Yes/No/Not sure answer
bothpartners[, table(pfrt27)]
bothpartners[, levels(pfrt27)]
bothpartners[, pfrt27 := factor(pfrt27, 
                               levels = c("-9 Invalid multiple answer",
                                          "-4 Filter error / Incorrect entry",
                                          "-3 Does not apply", "-2 No Answer",
                                          "1 Yes",
                                          "2 No",
                                          "5 I´m not sure",
                                          "6 I haven't thought about that"),
                               labels = c(-9, -4:-2, 1:2, 5:6))]
bothpartners[, table(pfrt27)]

# As in frt7
bothpartners[, table(pfrt7)]
bothpartners[, pfrt7 := factor(pfrt7, 
                              levels = c("-9 Invalid multiple answer",
                                         "-4 Filter error / Incorrect entry",
                                         "-3 Does not apply", "-2 No answer",
                                         "-1 Don't know", "1 Yes, definitely",
                                         "2 Yes, perhaps", "3 No, probably not",
                                         "4 No, definitely not",
                                         "7 I haven't thought about that yet"),
                              labels = c(-9, -4:-1, 1:4, 7))]
bothpartners[, table(pfrt7)]

# We do not combine pfrt6 and pfrt27 because they have different answers' categories

# answer they want no children: pfrt6 == 7 | pfrt27 == 2
bothpartners[, intentions_partner := ifelse(pfrt6 == 7 | pfrt27 == 2, "def_no", NA)]

# answer they want one or more children: pfrt27 == 1, pfrt6 = 1-4
bothpartners[, intentions_partner := ifelse(pfrt27 == 1 & pfrt7 == 1,
                                           "def_yes", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt27 == 1 & pfrt7 == 2,
                                           "prob_yes", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt27 == 1 & pfrt7 == 3,
                                           "prob_no", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt27== 1 & pfrt7 == 4,
                                           "def_no", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt27 == 1 & pfrt7 == 7,
                                           "prob_no", intentions_partner)]

bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 1,
                                            "def_yes", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 2,
                                            "prob_yes", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 3,
                                            "prob_no", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 4,
                                            "def_no", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 7,
                                            "prob_no", intentions_partner)]

# answer they are not sure about expected number of children
bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(5:6) & pfrt7 == 1,
                                           "def_yes", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(5:6) & pfrt7 == 2,
                                           "prob_yes", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(5:6) & pfrt7 == 3,
                                           "prob_no", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(5:6) & pfrt7 == 4,
                                           "def_no", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt6 %in% c(5:4) & pfrt7 == 7,
                                           "prob_no", intentions_partner)]

bothpartners[, intentions_partner := ifelse(pfrt27 %in% c(5:6) & pfrt7 == 1,
                                            "def_yes", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt27 %in% c(5:6) & pfrt7 == 2,
                                            "prob_yes", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt27 %in% c(5:6) & pfrt7 == 3,
                                            "prob_no", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt27 %in% c(5:6) & pfrt7 == 4,
                                            "def_no", intentions_partner)]
bothpartners[, intentions_partner := ifelse(pfrt27 %in% c(5:4) & pfrt7 == 7,
                                            "prob_no", intentions_partner)]
bothpartners[, describe(intentions_partner)]
bothpartners[is.na(intentions_partner), table(pfrt6, pfrt7)]
bothpartners[is.na(intentions_partner), table(pfrt27, pfrt7)]

# fertility intentions binary
bothpartners[!is.na(intentions_partner), 
             intentions_partner_binary := ifelse(intentions_partner %in% c("def_no", "prob_no"), 
                                                 "no", "yes")]
bothpartners[, describe(intentions_partner_binary)]

table(bothpartners$intentions_anchor_binary, bothpartners$intentions_partner_binary)

# ---- Alternative variable for sensitivity analysis: Intentions timing ---- 
# answer they want no children -> do not want children now nor later
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined == 7, "never", NA)]  # do not want children now or later

# answer they want one or more children
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(1:4) & frt7 == 1,
                                                      "now", intentions_anchor_timing)] # want it now
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(1:4) & frt7 == 2,
                                                      "now", intentions_anchor_timing)] # maybe want it now
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(1:4) & frt7 == 3,
                                                      "later", intentions_anchor_timing)] # probably not now
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(1:4) & frt7 == 4,
                                                      "later", intentions_anchor_timing)] # definetely not now
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(1:4) & frt7 == 7,
                                                      "later", intentions_anchor_timing)] # not sure, so probably later

# answer they are not sure about expected number of children
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(5:6) & frt7 == 1,
                                                      "now", intentions_anchor_timing)]
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(5:6) & frt7 == 2,
                                                      "now", intentions_anchor_timing)]
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(5:6) & frt7 == 3,
                                                      "later", intentions_anchor_timing)]
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(5:6) & frt7 == 4,
                                                      "later", intentions_anchor_timing)]
bothpartners[, intentions_anchor_timing := ifelse(frt6_combined %in% c(5:4) & frt7 == 7,
                                                      "later", intentions_anchor_timing)]

bothpartners[, describe(intentions_anchor_timing)]

# code the same for partner
# answer they want no children
bothpartners[, intentions_partner_timing := ifelse(pfrt6 == 7 | pfrt27 == 2, "never", NA)]
# answer they want one or more children
bothpartners[, intentions_partner_timing := ifelse(pfrt27 == 1 & pfrt7 == 1,
                                                       "now", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt27 == 1 & pfrt7 == 2,
                                                       "now", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt27 == 1 & pfrt7 == 3,
                                                       "later", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt27== 1 & pfrt7 == 4,
                                                       "later", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt27 == 1 & pfrt7 == 7,
                                                       "later", intentions_partner_timing)]

bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 1,
                                                       "now", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 2,
                                                       "now", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 3,
                                                       "later", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 4,
                                                       "later", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(1:4) & pfrt7 == 7,
                                                       "later", intentions_partner_timing)]

# answer they are not sure about expected number of children
bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(5:6) & pfrt7 == 1,
                                                       "now", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(5:6) & pfrt7 == 2,
                                                       "now", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(5:6) & pfrt7 == 3,
                                                       "later", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(5:6) & pfrt7 == 4,
                                                       "later", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt6 %in% c(5:4) & pfrt7 == 7,
                                                       "later", intentions_partner_timing)]

bothpartners[, intentions_partner_timing := ifelse(pfrt27 %in% c(5:6) & pfrt7 == 1,
                                                       "now", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt27 %in% c(5:6) & pfrt7 == 2,
                                                       "now", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt27 %in% c(5:6) & pfrt7 == 3,
                                                       "later", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt27 %in% c(5:6) & pfrt7 == 4,
                                                       "later", intentions_partner_timing)]
bothpartners[, intentions_partner_timing := ifelse(pfrt27 %in% c(5:4) & pfrt7 == 7,
                                                       "later", intentions_partner_timing)]
bothpartners[, describe(intentions_partner_timing)]

table(bothpartners$intentions_anchor_timing, bothpartners$intentions_partner_timing)

# Now, identify the first time that both partners in the couple have a valid answer
bothpartners[!is.na(intentions_anchor_binary) & !is.na(intentions_partner_binary), flag := 1]
table(bothpartners$flag)

bothpartners[flag == 1, n_valid_waves := 1:.N, by = cid]
# View(bothpartners[, c("id", "pid", "cid", "wave", "intentions_anchor_binary", 
#                       "intentions_partner_binary", "flag", "n_valid_waves")])

# keep only first wave when both partners answer the fertility question
bothpartners_baseline <- bothpartners[n_valid_waves == 1]

# save results
#save(bothpartners_baseline, file = "./Data/bothpartners_baseline_updated.Rda")

# ---- Cleaning other variables ----
# only heterosex relationships
table(bothpartners_baseline$homosex)
bothpartners_baseline <- bothpartners_baseline[homosex == "0 Heterosexual"]

# remove one case where psex_gen == -4
bothpartners_baseline <- bothpartners_baseline[psex_gen != "-4 Filter error / Incorrect entry"]

# keep only individuals with 0 kids
bothpartners_baseline_nokids <- bothpartners_baseline[nkids == 0]

# education of the anchor
bothpartners_baseline_nokids[, table(isced)]

# remove incomplete data
bothpartners_baseline_nokids <- bothpartners_baseline_nokids[isced != "-7 Incomplete data"]

bothpartners_baseline_nokids[, education_anchor := ifelse(isced == "0 currently enrolled",
                                                          "in_school", NA)]
bothpartners_baseline_nokids[, education_anchor := ifelse(isced %in% c("1 no degree (1b)",
                                                                       "2 lower secondary education (2b)",
                                                                       "3 lower secondary education (2a)"),
                                                          "No_degree/lower", education_anchor)]
bothpartners_baseline_nokids[, education_anchor := ifelse(isced %in% c("4 upper secondary education vocational (3b)",
                                                                       "5 upper secondary education general (3a)",
                                                                       "6 post-secondary non tertiary education general (4a)"),
                                                          "Upper_secondary", education_anchor)]
bothpartners_baseline_nokids[, education_anchor := ifelse(isced %in% c("7 first stage of tertiary education (5)",
                                                                       "8 second stage of tertiary education (6)"),
                                                          "Tertiary", education_anchor)]
bothpartners_baseline_nokids[, table(education_anchor)]

# education of the anchor
bothpartners_baseline_nokids[, table(pisced)]

# remove incomplete data
bothpartners_baseline_nokids <- bothpartners_baseline_nokids[pisced != "-7 Incomplete data"]

bothpartners_baseline_nokids[, education_partner := ifelse(pisced == "0 currently enrolled", "in_school", NA)]
bothpartners_baseline_nokids[, education_partner := ifelse(pisced %in% c("1 no degree (1b)", 
                                                                         "2 lower secondary education (2b)",
                                                                         "3 lower secondary education (2a)"),
                                                          "No_degree/lower", education_partner)]
bothpartners_baseline_nokids[, education_partner := ifelse(pisced %in% c("4 upper secondary education vocational (3b)",
                                                                         "5 upper secondary education general (3a)",
                                                                         "6 post-secondary non tertiary education general (4a)"),
                                                          "Upper_secondary", education_partner)]
bothpartners_baseline_nokids[, education_partner := ifelse(pisced %in% c("7 first stage of tertiary education (5)",
                                                                         "8 second stage of tertiary education (6)"),
                                                          "Tertiary", education_partner)]
bothpartners_baseline_nokids[, table(education_partner)]

# Age at baseline should be >= 18
# function to convert from factor to numeric
as.double.factor <- function(x) {as.numeric(levels(x))[x]}

table(bothpartners_baseline_nokids$age)
levels(bothpartners_baseline_nokids$age)
bothpartners_baseline_nokids[, age := factor(age, 
                                             levels = levels(age),
                                             labels = c(-3, -5, -7, 14:50))]
bothpartners_baseline_nokids[, agea := as.double.factor(age)]
table(bothpartners_baseline_nokids$agea)

table(bothpartners_baseline_nokids$page)
levels(bothpartners_baseline_nokids$page)
bothpartners_baseline_nokids[, page := factor(page, 
                                              levels = levels(page),
                                              labels = c(-3, -5, -7, 10, 107, 11:76, 78, 80, 81, 83, 99))]
bothpartners_baseline_nokids[, agep := as.double.factor(page)]
table(bothpartners_baseline_nokids$agep)

# restrict to both partners older than 17 at baseline
bothpartners_baseline_nokids_agerest <- bothpartners_baseline_nokids[agea > 17 & agep > 17]

## relationship satisfaction
table(bothpartners_baseline_nokids_agerest$sat3)
levels(bothpartners_baseline_nokids_agerest$sat3)
bothpartners_baseline_nokids_agerest[, relsat := factor(sat3, 
                                                levels = levels(sat3),
                                                labels = c(-1, -12, -2, -3, -4, -5,
                                                           0, 1, 10, 2:9))]
bothpartners_baseline_nokids_agerest[, relsat := as.double.factor(relsat)]
table(bothpartners_baseline_nokids_agerest$relsat)

table(bothpartners_baseline_nokids_agerest$psat3)
levels(bothpartners_baseline_nokids_agerest$psat3)
bothpartners_baseline_nokids_agerest[, prelsat := factor(psat3, 
                                                         levels = levels(psat3),
                                                         labels = c(-1, -10, -2, -3, -4, -9,
                                                                    0, 1, 10, 2:9))]
bothpartners_baseline_nokids_agerest[, prelsat := as.double.factor(prelsat)]
table(bothpartners_baseline_nokids_agerest$prelsat)

# restrict to valid answers from both partners
bothpartners_baseline_nokids_agerest_relsat <- 
  bothpartners_baseline_nokids_agerest[relsat >= 0 & prelsat >= 0]

# ---- Code gendered variables ----
table(bothpartners_baseline_nokids_agerest_relsat$sex_gen)
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "2 Female", relsatf := relsat]
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "2 Female", relsatm := prelsat]
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "1 Male", relsatf := prelsat]
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "1 Male", relsatm := relsat]

table(bothpartners_baseline_nokids_agerest_relsat$education_anchor, 
      bothpartners_baseline_nokids_agerest_relsat$education_partner)


bothpartners_baseline_nokids_agerest_relsat[sex_gen == "1 Male", 
                                            educdiff := ifelse(education_anchor == "Tertiary" &
                                                                 education_partner == "Tertiary",
                                                               "both_tertiary", NA)]
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "1 Male", 
                                            educdiff := ifelse(education_anchor != "Tertiary" & 
                                                                  education_partner == "Tertiary",
                                                               "her_tertiary", educdiff)]
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "1 Male", 
                                            educdiff := ifelse(education_anchor == "Tertiary" & 
                                                                 education_partner != "Tertiary",
                                                               "him_tertiary", educdiff)]
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "1 Male", 
                                            educdiff := ifelse(education_anchor != "Tertiary" & 
                                                                 education_partner != "Tertiary",
                                                               "both_school_lower_secondary", educdiff)]

bothpartners_baseline_nokids_agerest_relsat[sex_gen == "2 Female", 
                                            educdiff := ifelse(education_anchor == "Tertiary" &
                                                                 education_partner == "Tertiary",
                                                               "both_tertiary", educdiff)]
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "2 Female", 
                                            educdiff := ifelse(education_anchor != "Tertiary" & 
                                                                 education_partner == "Tertiary",
                                                               "him_tertiary", educdiff)]
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "2 Female", 
                                            educdiff := ifelse(education_anchor == "Tertiary" & 
                                                                 education_partner != "Tertiary",
                                                               "her_tertiary", educdiff)]
bothpartners_baseline_nokids_agerest_relsat[sex_gen == "2 Female", 
                                            educdiff := ifelse(education_anchor != "Tertiary" & 
                                                                 education_partner != "Tertiary",
                                                               "both_school_lower_secondary", educdiff)]
describe(bothpartners_baseline_nokids_agerest_relsat$educdiff)


bothpartners_baseline_nokids_agerest_relsat[, sex := ifelse(sex_gen == "1 Male", "male", "female")]

bothpartners_baseline_nokids_agerest_relsat[, agediff := agea - agep]
hist(bothpartners_baseline_nokids_agerest_relsat$agediff)
summary(bothpartners_baseline_nokids_agerest_relsat$agediff)

bothpartners_baseline_nokids_agerest_relsat[, agegap := ifelse(agediff >= -2 & agediff <= 2, 
                                                               "age_homogeneous", NA)]
bothpartners_baseline_nokids_agerest_relsat[sex == "male", agegap := ifelse(agediff > 2, 
                                                                            "man_older", agegap)]
bothpartners_baseline_nokids_agerest_relsat[sex == "male", agegap := ifelse(agediff < -2, 
                                                                            "woman_older", agegap)]
bothpartners_baseline_nokids_agerest_relsat[sex == "female", agegap := ifelse(agediff > 2, 
                                                                              "woman_older", agegap)]
bothpartners_baseline_nokids_agerest_relsat[sex == "female", agegap := ifelse(agediff < -2, 
                                                                              "man_older", agegap)]

describe(bothpartners_baseline_nokids_agerest_relsat$agegap)
bothpartners_baseline_nokids_agerest_relsat[, agegap := factor(agegap)]

data_baseline <- bothpartners_baseline_nokids_agerest_relsat

# Finally, intentions gendered

data_baseline[, agreement_intentions := ifelse(intentions_anchor_binary == "yes" &
                                                 intentions_partner_binary == "yes",
                                               "both_yes", NA)]
data_baseline[, agreement_intentions := ifelse(intentions_anchor_binary == "no" &
                                                 intentions_partner_binary == "no",
                                               "both_no", agreement_intentions)]
data_baseline[sex == "male", 
              agreement_intentions := ifelse(intentions_anchor_binary == "no" &
                                               intentions_partner_binary == "yes",
                                             "he_no_she_yes", agreement_intentions)]
data_baseline[sex == "male", 
              agreement_intentions := ifelse(intentions_anchor_binary == "yes" &
                                               intentions_partner_binary == "no",
                                             "he_yes_she_no", agreement_intentions)]
data_baseline[sex == "female", 
              agreement_intentions := ifelse(intentions_anchor_binary == "yes" &
                                               intentions_partner_binary == "no",
                                             "he_no_she_yes", agreement_intentions)]
data_baseline[sex == "female", 
              agreement_intentions := ifelse(intentions_anchor_binary == "no" &
                                               intentions_partner_binary == "yes",
                                             "he_yes_she_no", agreement_intentions)]
describe(data_baseline$agreement_intentions)
data_baseline[, agreement_intentions := factor(agreement_intentions,
                                               levels = c("both_no", "both_yes",
                                                          "he_yes_she_no", "he_no_she_yes"))]
# date of interview
data_baseline[, intdat := ((inty - 1900) * 12 + intm - 1)]

save(data_baseline, file = "Data/data_baseline_ready_updated.Rda")
