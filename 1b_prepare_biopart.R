# ----- Prepara data for multistate model -----
## In this script:
## 1. load biopart data
## 2. load biochild data
## 3. load cleaned baseline couple data
## 4. Find date of birth of first child to couple
## 5. Correctly remove couples with children before entry into risk set
## 6. restrict biopart data to only couple considered
## 7. find timing of transitions
## 8. transform in multistate format


# ---- Libraries ----
library(mstate)
library(readstata13)
library(data.table)
library(survival)
library(Hmisc)

# ---- Data ----
load("Data/data_baseline_ready_updated.Rda")

biopart <- read.dta13("Data/biopart.dta",
                      convert.factors = FALSE,
                      replace.strl = TRUE)

# labels of the variables as in stata
labels <- get.label.tables(biopart)

biopart <- as.data.table(biopart)

biopart_sel <- biopart[id %in% data_baseline$id]
biopart_sel <- biopart[pid %in% data_baseline$pid]

biopart_sel <- biopart_sel[, c("id", "pid", "sexp", "relbeg", "relend", "cohbeg",
                               "cohend", "marbeg", "marend")]

# children
biochild <- read.dta13("Data/biochild.dta",
                       convert.factors = FALSE,
                       replace.strl = TRUE)
# labels of the variables as in stata
labels <- get.label.tables(biochild)

biochild <- as.data.table(biochild)

biochild_sel <- biochild[id %in% data_baseline$id]

# In data_baseline add a variable counting to 36 months after baseline interview
data_baseline[, intdat_plus36 := intdat + 36]

# transform wave variables in integer
data_baseline[, wave_int := as.integer(wave)]
table(data_baseline$wave, data_baseline$wave_int) # check

# ---- Merge biochild with data_baseline (only selected variables) ----
biochild_merge <- merge(biochild_sel, 
                        data_baseline[, c("id", "wave_int", "intdat", "intdat_plus36")],
                        by = "id")
setorder(biochild_merge, id, wave)
# cid is child id not couple id
View(biochild_merge[, c("id", "cid", "wave", "wave_int", "dobk", "intdat.y", "intdat_plus36")])

biochild_merge_sel <- biochild_merge[dobk >= intdat.y & dobk <= intdat_plus36]
nrow(biochild_merge[dobk < intdat.y])
View(biochild_merge[dobk < intdat.y, c("id", "cid", "wave", "wave_int", "dobk", 
                                       "intdat.y", "intdat_plus36")])
biochild_merge[dobk < intdat.y, describe(id)] # 59 individuals whose children were
# born before the entry into the risk set in our model - this was not captured by
# the variable nkids so we need to correct for it

flagchild <- biochild_merge[dobk < intdat.y & wave < wave_int]
describe(flagchild$id)

biochild_merge_sel[, kid_wave := 1:.N, by = c("cid")]
biochild_merge_sel <- biochild_merge_sel[kid_wave == 1]
biochild_merge_sel_dates <- biochild_merge_sel[, c("id", "pid", "cid", "dobk")]

biochild_merge_sel_dates[, nkids_couple := 1:.N, by = c("id", "pid")]

# we only need the first child born to a couple
biochild_firstonly <- biochild_merge_sel_dates[nkids_couple == 1]

# ---- Merge biopart with data_baseline ----
biopart_merge <- merge(biopart_sel, data_baseline, by = c("id", "pid"))
biopart_merge <- merge(biopart_merge, biochild_firstonly, 
                       by = c("id", "pid"), all.x = TRUE)

names(biopart_merge)

biopart_sel_covs <- biopart_merge[, c("id", "pid", "sex", "sexp", "wave", "wave_int",
                                      "demodiff.x", "sample.x", "doby_gen", "dobm_gen",
                                      "cohort", "np", "infertile", "pregnant", "east",
                                      "livewithpar", "intentions_anchor", "intentions_partner",
                                      "agea", "agep", "relsatf", "relsatm", "educdiff",
                                      "agediff", "agegap", "agreement_intentions",
                                      "intentions_anchor_timing", "intentions_partner_timing",
                                      "intdat", "intdat_plus36", 
                                      "relbeg", "relend", "cohbeg", "cohend", 
                                      "marbeg", "marend", "dobk")]
data_biopart <- biopart_sel_covs
setnames(data_biopart, c("demodiff.x", "sample.x", "intdat", "intdat_plus36", 
                         "dobk"),
         c("demodiff", "sample", "baseint", "baseint_plus3", "kiddob"))

# ---- Clean the longitudinal format and prepare for mstate ----

# (cohabitations that end before end of the follow-up)
table((data_biopart$relend > 0) & (data_biopart$relend <= data_biopart$baseint_plus3))

# among these, how many end before the relationship ends?
table((data_biopart$cohend > 0) & 
        (data_biopart$cohend > 0) & 
        (data_biopart$cohend <= data_biopart$baseint_plus3))
table((data_biopart$cohend > 0) & (data_biopart$cohend > 0) & 
        (data_biopart$cohend <= data_biopart$baseint_plus3) &
        (data_biopart$cohend < data_biopart$relend))

# start of the observation period is the maximum between the relationship's beginning
# date and the baseline interview date
data_biopart[, start := pmax(relbeg, baseint), by = "id"]

data_biopart[, cid := 1:.N] # this is now the couple id in the new dataset

nrow(data_biopart[relbeg == -7]) # error in relationship beginning date
data_biopart <- data_biopart[relbeg != -7]

data_biopart[, reldur_base := baseint - relbeg]
View(data_biopart[reldur_base < 0]) # 1 case where the baseline interview was 1 month before the
# relationship begin. It could be that this was a coding error. We assign 0 to the baseline duration
data_biopart[reldur_base < 0, reldur_base := 0]
View(data_biopart[reldur_base == 0])

data_biopart[is.na(kiddob), kiddob := -99]

## Note: The following is probably a very cumbersome and unnecessary way to code
##       the end of the observation period. However, it has helped me to be more
##       precise about it, thinking about all the possible transitions in the 
##       data.
## Additionally: many of these transitions will actually not be captured because
## they happen after the 36 months period from the baseline interview.

# Case 1: Married couple with kid splits
View(data_biopart[cohend > 0 & marend > 0 & kiddob > 0 & relend > 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend > 0 & marend > 0 & kiddob > 0 & relend > 0, 
             end := pmin(baseint_plus3, cohend, marend, kiddob, relend)]

# Case 2: Non-divorce with kid
View(data_biopart[cohend > 0 & marend < 0 & kiddob > 0 & relend > 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend > 0 & marend < 0 & kiddob > 0 & relend > 0, 
             end := pmin(baseint_plus3, cohend, kiddob, relend)]

# Case 3: Non-divorce, no kid, cohabitation ends
View(data_biopart[cohend > 0 & marend < 0 & kiddob < 0 & relend > 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend > 0 & marend < 0 & kiddob < 0 & relend > 0, 
             end := pmin(baseint_plus3, cohend, relend)]

# Case 4: Non-divorce, kid, cohabitation ends, but relationship does not
View(data_biopart[cohend > 0 & marend < 0 & kiddob > 0 & relend < 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
# none
# data_biopart[cohend > 0 & marend < 0 & kiddob > 0 & relend < 0, 
#              end := pmin(baseint_plus3, cohend, kiddob)]

# Case 4: Non-divorce, no kid, cohabitation ends, but relationship does not
View(data_biopart[cohend > 0 & marend < 0 & kiddob < 0 & relend < 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend > 0 & marend < 0 & kiddob < 0 & relend < 0, 
             end := pmin(baseint_plus3, cohend)]

# Case 6: Divorce, no kids
View(data_biopart[cohend > 0 & marend > 0 & kiddob < 0 & relend > 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend > 0 & marend > 0 & kiddob < 0 & relend > 0, 
             end := pmin(baseint_plus3, cohend, marend, relend)]

# Case 7: Divorce, no kids, relationship does not end
View(data_biopart[cohend > 0 & marend > 0 & kiddob < 0 & relend < 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
# none
# data_biopart[cohend > 0 & marend > 0 & kiddob < 0 & relend < 0, 
#              end := pmin(baseint_plus3, cohend, marend)]

## Case 8: Divorce, with kids, relationship does not end
View(data_biopart[cohend > 0 & marend > 0 & kiddob > 0 & relend < 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
# none
# data_biopart[cohend > 0 & marend > 0 & kiddob > 0 & relend < 0, 
#              end := pmin(baseint_plus3, cohend, marend, kiddob)]

# Case 9: Marriage ends, with kids, but cohabitation does not
View(data_biopart[cohend < 0 & marend > 0 & kiddob > 0 & relend > 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
# none
# data_biopart[cohend < 0 & marend > 0 & kiddob > 0 & relend > 0, 
#              end := pmin(baseint_plus3, marend, kiddob, relend)]

# Case 10: Marriage ends, no kids, but cohabitation does not
View(data_biopart[cohend < 0 & marend > 0 & kiddob < 0 & relend > 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
# none
#data_biopart[cohend < 0 & marend > 0 & kiddob < 0 & relend > 0, 
#             end := pmin(baseint_plus3, marend, relend)]

# Case 11: Marriage ends, with kids, but cohabitation and relation do not
View(data_biopart[cohend < 0 & marend > 0 & kiddob > 0 & relend < 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend < 0 & marend > 0 & kiddob > 0 & relend < 0, 
             end := pmin(baseint_plus3, marend, kiddob)]

# Case 12: Marriage ends, no kids, but cohabitation and relation do not
View(data_biopart[cohend < 0 & marend > 0 & kiddob < 0 & relend < 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend < 0 & marend > 0 & kiddob < 0 & relend < 0, 
             end := pmin(baseint_plus3, marend)]

# Case 13: cohabitation & marriage do not end, kids and relationship ends
View(data_biopart[cohend < 0 & marend < 0 & kiddob > 0 & relend > 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend < 0 & marend < 0 & kiddob > 0 & relend > 0, 
             end := pmin(baseint_plus3, kiddob, relend)]

# Case 14: cohabitation & marriage do not end, no kid and relationship ends
View(data_biopart[cohend < 0 & marend < 0 & kiddob < 0 & relend > 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend < 0 & marend < 0 & kiddob < 0 & relend > 0, 
             end := pmin(baseint_plus3, relend)]

# Case 15: cohabitation & marriage do not end, with kids and relationship does not end
View(data_biopart[cohend < 0 & marend < 0 & kiddob > 0 & relend < 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend < 0 & marend < 0 & kiddob > 0 & relend < 0, 
             end := pmin(baseint_plus3, kiddob)]

# Case 16: no (absorbing) events
View(data_biopart[cohend < 0 & marend < 0 & kiddob < 0 & relend < 0, 
                  c("cid", "start", "relbeg", "cohbeg", "marbeg", "relend",
                    "cohend", "marend", "kiddob", "baseint_plus3")])
data_biopart[cohend < 0 & marend < 0 & kiddob < 0 & relend < 0, 
             end := pmin(baseint_plus3)]

describe(data_biopart$end)

table(data_biopart$start == data_biopart$end)
View(data_biopart[start == end, 
                  c("cid", "id", "pid", "relbeg", "cohbeg",
                    "marbeg", "kiddob", "relend", "cohend",
                    "marend", "start", "baseint", "baseint_plus3", "end")])

data_biopart <- data_biopart[!(start == end)]

# Event-specific variables
data_biopart[, cohabitation := ifelse(cohbeg > 0 & cohbeg <= baseint_plus3, 1, 0)]
data_biopart[, marriage := ifelse(marbeg > 0 & marbeg <= baseint_plus3, 1, 0)]
data_biopart[, kid := ifelse(kiddob > 0 & kiddob <= baseint_plus3, 1, 0)]
data_biopart[, dissolution := (relend > 0 & relend <= baseint_plus3 |
                                 cohend > 0 & cohend <= baseint_plus3 |
                                 marend > 0 & marend <= baseint_plus3) * 1]

describe(data_biopart$cohabitation)  # 2070
describe(data_biopart$marriage)      # 770
describe(data_biopart$dissolution)   # 499
describe(data_biopart$kid)           # 286

# find in which state everyone starts
data_biopart[, start_state := ifelse(cohabitation == 1 & cohbeg <= start, 2, NA)]
data_biopart[, start_state := ifelse(marriage == 1 & marbeg <= start, 3, start_state)]
data_biopart[is.na(start_state), start_state := 1]
table(data_biopart$start_state)
# > table(data_biopart$start_state)
# 
# 1    2    3 
# 1007 1257  420 

# correction for simultaneous transitions
nrow(data_biopart[cohbeg > 0 & cohbeg == marbeg]) # 56
nrow(data_biopart[cohbeg > 0 & cohbeg == kiddob]) # 2
nrow(data_biopart[marbeg > 0 & marbeg == kiddob]) # 2
nrow(data_biopart[start == kiddob]) # 0
nrow(data_biopart[start == cohbeg]) # 76

nrow(data_biopart[start > end]) # 23 remove
data_biopart <- data_biopart[start <= end] 

# we artificially assign one more month to the transition that happens later
data_biopart[cohbeg > 0 & cohbeg == marbeg, marbeg := marbeg + 1]
data_biopart[cohbeg > 0 & cohbeg == kiddob, kiddob := kiddob + 1]
data_biopart[marbeg > 0 & marbeg == kiddob, kiddob := kiddob + 1]

# Replace last possible date if relationship transition happens later
cols_to_replace <- c("cohbeg", "marbeg", "kiddob", "relend", "cohend", "marend")
data_biopart[, (cols_to_replace) := lapply(.SD, function(x) fifelse(x > baseint_plus3, baseint_plus3, x)),
             .SDcols = cols_to_replace]

# Replace negative values in relationship dates with end of the observation period
cols_to_replace <- c("cohbeg", "marbeg", "kiddob", "relend", "cohend", "marend")
data_biopart[, (cols_to_replace) := lapply(.SD, function(x) fifelse(x < 0, baseint_plus3, x)),
             .SDcols = cols_to_replace]

# View(data_biopart[, c("cid", "id", "pid", "relbeg", "cohbeg",
#                       "marbeg", "kiddob", "relend", "cohend",
#                       "marend", "start", "start_state", "cohabitation", "marriage",
#                       "dissolution", "kid", "baseint_plus3")])
# 
table(data_biopart[start_state == 2, cohabitation]) # 1236
table(data_biopart[start_state == 3, cohabitation]) # 417
table(data_biopart[start_state == 3, marriage]) # 419
# we do not count events that happened before the entry into the risk set
data_biopart[start_state == 2, cohabitation := 0]
data_biopart[start_state == 3, cohabitation := 0]
data_biopart[start_state == 3, marriage := 0]

# code interaction between agreement and duration
data_biopart[, hist(reldur_base)]
data_biopart[, summary(reldur_base)]
data_biopart$reldur_base_cat <- cut(data_biopart$reldur_base, c(0, 24, 389),
                                    include.lowest = T)
describe(data_biopart$reldur_base_cat)
describe(data_biopart$agreement_intentions)

data_biopart[, agreement_duration := ifelse((agreement_intentions == "both_no" & 
                                               reldur_base_cat == "[0,24]"),
                                            "BothNo_le2", NA)]
data_biopart[, agreement_duration := ifelse((agreement_intentions == "both_no" & 
                                               reldur_base_cat == "(24,389]"),
                                            "BothNo_2plus", agreement_duration)]
data_biopart[, agreement_duration := ifelse((agreement_intentions == "both_yes" & 
                                               reldur_base_cat == "[0,24]"),
                                            "BothYes_le2", agreement_duration)]
data_biopart[, agreement_duration := ifelse((agreement_intentions == "both_yes" & 
                                               reldur_base_cat == "(24,389]"),
                                            "BothYes_2plus", agreement_duration)]
data_biopart[, agreement_duration := ifelse((agreement_intentions == "he_yes_she_no" & 
                                               reldur_base_cat == "[0,24]"),
                                            "HeYesSheNo_le2", agreement_duration)]
data_biopart[, agreement_duration := ifelse((agreement_intentions == "he_yes_she_no" & 
                                               reldur_base_cat == "(24,389]"),
                                            "HeYesSheNo_2plus", agreement_duration)]
data_biopart[, agreement_duration := ifelse((agreement_intentions == "he_no_she_yes" & 
                                               reldur_base_cat == "[0,24]"),
                                            "SheYesHeNo_le2", agreement_duration)]
data_biopart[, agreement_duration := ifelse((agreement_intentions == "he_no_she_yes" & 
                                               reldur_base_cat == "(24,389]"),
                                            "SheYesHeNo_2plus", agreement_duration)]

describe(data_biopart$agreement_duration)

data_biopart$agreement_duration <- factor(data_biopart$agreement_duration, 
                                          levels = c("BothNo_le2", "BothNo_2plus", "HeYesSheNo_le2", "HeYesSheNo_2plus",
                                                     "SheYesHeNo_le2", "SheYesHeNo_2plus", "BothYes_le2", "BothYes_2plus"))

# remove missing values in livewithpar
data_biopart <- data_biopart[!is.na(livewithpar)]

# Code agreement with timing
data_biopart[, table(intentions_anchor_timing, intentions_partner_timing)]
data_biopart[, agreement_timing := ifelse(intentions_anchor_timing == "never" & 
                                            intentions_partner_timing == "never",
                                          "both_never", NA)]
data_biopart[, agreement_timing := ifelse(intentions_anchor_timing == "now" &
                                            intentions_partner_timing == "now",
                                          "both_now", agreement_timing)]
data_biopart[, agreement_timing := ifelse(intentions_anchor_timing == "later" & 
                                            intentions_partner_timing == "later",
                                          "both_later", agreement_timing)]
data_biopart[, agreement_timing := ifelse(intentions_anchor_timing == "now" & 
                                            intentions_partner_timing == "later",
                                          "timing_misalignment", agreement_timing)]
data_biopart[, agreement_timing := ifelse(intentions_anchor_timing == "later" & 
                                            intentions_partner_timing == "now",
                                          "timing_misalignment", agreement_timing)]
data_biopart[, agreement_timing := ifelse(intentions_anchor_timing == "never" & 
                                            intentions_partner_timing != "never",
                                          "fundamental_misalignment", agreement_timing)]
data_biopart[, agreement_timing := ifelse(intentions_anchor_timing != "never" & 
                                            intentions_partner_timing == "never",
                                          "fundamental_misalignment", agreement_timing)]
describe(data_biopart$agreement_timing)
data_biopart[, agreement_timing := factor(agreement_timing, 
                                          levels = c("both_never", "both_later",
                                                     "both_now", "fundamental_misalignment",
                                                     "timing_misalignment"))]
data_biopart[, table(agreement_intentions, agreement_timing)]

# ---- Transform in mstate format ----
# Specify transition matrix
tmat <- transMat(x = list(c(2, 3, 4, 5), c(3, 4, 5), c(4, 5), c(), c()),
                 names = c("dating", "cohabitation", "marriage",
                           "dissolution", "kid"))
keep <- c("start", "start_state", "baseint_plus3", "end", "baseint", "sex",
          "agreement_intentions", "livewithpar", "educdiff", "agea", "agep", "agediff",
          "agegap", "reldur_base", "relsatm", "relsatf", "relbeg", "relend", "cohbeg",
          "cohend", "marbeg", "reldur_base_cat", "agreement_duration", "agreement_timing",
          "marend", "kiddob", "cohabitation", "marriage", "dissolution", "kid", "id", "pid")
data_biopart_mstate <- msprep(time = c(NA, "cohbeg", "marbeg", "relend", "kiddob"),
                              status = c(NA, "cohabitation", "marriage", "dissolution", "kid"),
                              id = "cid", 
                              keep = keep,
                              start = list(state = data_biopart$start_state, 
                                           time = data_biopart$start),
                              data = data_biopart, trans = tmat)

covs <- c("sex", "agreement_intentions", "livewithpar", "educdiff", "agea", "agegap",
          "reldur_base", "relsatm", "relsatf", "reldur_base_cat", "agreement_duration",
          "agreement_timing")
data <- expand.covs(data_biopart_mstate, covs, longnames = FALSE)

data$Tstart_from_int1 <- data_biopart_mstate$Tstart - data_biopart_mstate$start
data$Tstop_from_int1 <- data_biopart_mstate$Tstop - data_biopart_mstate$start

save(data, file = "Data/data_mstate_ready_updated.Rda")
