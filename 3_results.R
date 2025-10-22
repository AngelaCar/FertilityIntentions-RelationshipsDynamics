# ---- Results from best fitting model ----
library(mstate)
library(data.table)
library(survival)
library(Hmisc)
library(RColorBrewer)

# ---- Figures ----
load("Results/bestmodel_updated.Rda")
load("Results/datafitting_updated.Rda")

# Function to obtain HRs and their CIs from model with interaction
get_combined_HR <- function(fit, terms, conf.level = 0.95) {
  # fit   : coxph object
  # terms : character vector of coefficient names to SUM (e.g. c("a","b","a:b"))
  # returns data.frame with logHR, SE, HR, CI
  
  # Check names
  all_names <- names(coef(fit))
  if (!all(terms %in% all_names)) {
    missing <- terms[!terms %in% all_names]
    stop("These terms not found in model coefficients: ", paste(missing, collapse = ", "))
  }
  
  beta_sub <- coef(fit)[terms] # vector of selected betas
  V_sub <- vcov(fit)[terms, terms, drop = FALSE] # submatrix of vcov
  
  a <- rep(1, length(beta_sub)) # weights (1 for sum)
  logHR <- sum(beta_sub) # a^T beta
  var_logHR <- as.numeric(t(a) %*% V_sub %*% a) # a^T V a  (includes covariances)
  se_logHR <- sqrt(var_logHR)
  
  z <- qnorm((1 + conf.level) / 2)
  ci_low_log <- logHR - z * se_logHR
  ci_high_log <- logHR + z * se_logHR
  
  SE_HR<- exp(logHR) * se_logHR
  
  data.frame(
    Combination = paste(terms, collapse = " + "),
    logHR = logHR,
    SE = se_logHR,
    HR = exp(logHR),
    SE_HR = SE_HR,
    Lower = exp(ci_low_log),
    Upper = exp(ci_high_log),
    row.names = NULL
  )
}

s_bestmod <- summary(bestmod)
coeff_bestmod <- s_bestmod$coefficients
ci_bestmod <- s_bestmod$conf.int

HR <- exp(coeff_bestmod[,1])
SE_HR<- HR * coeff_bestmod[,2]

HR_SEs <- cbind(HR, SE_HR)

#* Calculate HR and their Confidence intervals for each transition, by combination
#* of exposure and duration of the relationship

# Transition 1
both_no_2plus_tr1 <- get_combined_HR(
  fit = bestmod,
  terms = c("reldur_base_cat.1")
)

both_yes_le2_tr1 <- get_combined_HR(
  fit = bestmod,
  terms = c("agreement_intentions1.1")
)

both_yes_2plus_tr1 <- get_combined_HR(
  fit = bestmod,
  terms = c(
    "agreement_intentions1.1",
    "reldur_base_cat.1",
    "agreement_intentions1.1:reldur_base_cat.1"
  )
)

him_yes_le2_tr1 <- get_combined_HR(bestmod, "agreement_intentions2.1")

him_yes_2plus_tr1 <- get_combined_HR(bestmod, c(
  "agreement_intentions2.1",
  "reldur_base_cat.1",
  "reldur_base_cat.1:agreement_intentions2.1"
))

her_yes_le2_tr1 <- get_combined_HR(bestmod, "agreement_intentions3.1")

her_yes_2plus_tr1 <- get_combined_HR(bestmod, c(
  "agreement_intentions3.1",
  "reldur_base_cat.1",
  "reldur_base_cat.1:agreement_intentions3.1"
))

tr1 <- rbind(
  both_no_2plus_tr1, 
  him_yes_le2_tr1, him_yes_2plus_tr1, her_yes_le2_tr1,
  her_yes_2plus_tr1, both_yes_le2_tr1, both_yes_2plus_tr1
)

rng1 <- range(tr1$Lower, tr1$Upper)

# Transition 3
both_no_2plus_tr3 <- get_combined_HR(
  fit = bestmod,
  terms = c("reldur_base_cat.3")
)

both_yes_le2_tr3 <- get_combined_HR(
  fit = bestmod,
  terms = c("agreement_intentions1.3")
)

both_yes_2plus_tr3 <- get_combined_HR(
  fit = bestmod,
  terms = c(
    "agreement_intentions1.3",
    "reldur_base_cat.3",
    "agreement_intentions1.3:reldur_base_cat.3"
  )
)

him_yes_le2_tr3 <- get_combined_HR(bestmod, "agreement_intentions2.3")
him_yes_2plus_tr3 <- get_combined_HR(bestmod, c(
  "agreement_intentions2.3",
  "reldur_base_cat.3",
  "reldur_base_cat.3:agreement_intentions2.3"
))

her_yes_le2_tr3 <- get_combined_HR(bestmod, "agreement_intentions3.3")
her_yes_2plus_tr3 <- get_combined_HR(bestmod, c(
  "agreement_intentions3.3",
  "reldur_base_cat.3",
  "reldur_base_cat.3:agreement_intentions3.3"
))

tr3 <- rbind(
  both_no_2plus_tr3, 
  him_yes_le2_tr3, him_yes_2plus_tr3, her_yes_le2_tr3,
  her_yes_2plus_tr3, both_yes_le2_tr3, both_yes_2plus_tr3
)
round(tr3[,2:7], 2)

rng3 <- range(tr3$Lower, tr3$Upper)

# Transition 5
both_yes_2plus_tr5 <- get_combined_HR(
  fit = bestmod,
  terms = c(
    "agreement_intentions1.5",
    "reldur_base_cat.5",
    "agreement_intentions1.5:reldur_base_cat.5"
  )
)

both_no_2plus_tr5 <- get_combined_HR(
  fit = bestmod,
  terms = c("reldur_base_cat.5")
)

both_yes_le2_tr5 <- get_combined_HR(
  fit = bestmod,
  terms = c("agreement_intentions1.5")
)

him_yes_le2_tr5 <- get_combined_HR(bestmod, "agreement_intentions2.5")
him_yes_2plus_tr5 <- get_combined_HR(bestmod, c(
  "agreement_intentions2.5",
  "reldur_base_cat.5",
  "reldur_base_cat.5:agreement_intentions2.5"
))

her_yes_le2_tr5 <- get_combined_HR(bestmod, "agreement_intentions3.5")
her_yes_2plus_tr5 <- get_combined_HR(bestmod, c(
  "agreement_intentions3.5",
  "reldur_base_cat.5",
  "reldur_base_cat.5:agreement_intentions3.5"
))

tr5 <- rbind(
  both_no_2plus_tr5,
  him_yes_le2_tr5, him_yes_2plus_tr5, her_yes_le2_tr5,
  her_yes_2plus_tr5, both_yes_le2_tr5, both_yes_2plus_tr5
)
round(tr5[,2:7], 2)

rng5 <- range(tr5$Lower, tr5$Upper)

# Transition 6 (new)
both_yes_2plus_tr6 <- get_combined_HR(
  fit = bestmod,
  terms = c(
    "agreement_intentions1.6.new",
    "reldur_base_cat.6.new",
    "agreement_intentions1.6.new:reldur_base_cat.6.new"
  )
)

both_no_2plus_tr6 <- get_combined_HR(
  fit = bestmod,
  terms = c("reldur_base_cat.6.new")
)

both_yes_le2_tr6 <- get_combined_HR(
  fit = bestmod,
  terms = c("agreement_intentions1.6.new")
)

him_yes_le2_tr6 <- get_combined_HR(bestmod, "agreement_intentions2.6.new")
him_yes_2plus_tr6 <- get_combined_HR(bestmod, c(
  "agreement_intentions2.6.new",
  "reldur_base_cat.6.new",
  "reldur_base_cat.6.new:agreement_intentions2.6.new"
))

her_yes_le2_tr6 <- get_combined_HR(bestmod, "agreement_intentions3.6.new")
her_yes_2plus_tr6 <- get_combined_HR(bestmod, c(
  "agreement_intentions3.6.new",
  "reldur_base_cat.6.new",
  "reldur_base_cat.6.new:agreement_intentions3.6.new"
))

tr6 <- rbind(
  both_no_2plus_tr6,
  him_yes_le2_tr6, him_yes_2plus_tr6, her_yes_le2_tr6,
  her_yes_2plus_tr6, both_yes_le2_tr6, both_yes_2plus_tr6
)
round(tr6[,2:7], 2)

rng6 <- range(tr6$Lower, tr6$Upper)

# Transition 7 (new)
both_yes_2plus_tr7 <- get_combined_HR(
  fit = bestmod,
  terms = c(
    "agreement_intentions1.7.new",
    "reldur_base_cat.7.new",
    "agreement_intentions1.7.new:reldur_base_cat.7.new"
  )
)

both_no_2plus_tr7 <- get_combined_HR(
  fit = bestmod,
  terms = c("reldur_base_cat.7.new")
)

both_yes_le2_tr7 <- get_combined_HR(
  fit = bestmod,
  terms = c("agreement_intentions1.7.new")
)

him_yes_le2_tr7 <- get_combined_HR(bestmod, "agreement_intentions2.7.new")
him_yes_2plus_tr7 <- get_combined_HR(bestmod, c(
  "agreement_intentions2.7.new",
  "reldur_base_cat.7.new",
  "reldur_base_cat.7.new:agreement_intentions2.7.new"
))

her_yes_le2_tr7 <- get_combined_HR(bestmod, "agreement_intentions3.7.new")
her_yes_2plus_tr7 <- get_combined_HR(bestmod, c(
  "agreement_intentions3.7.new",
  "reldur_base_cat.7.new",
  "reldur_base_cat.7.new:agreement_intentions3.7.new"
))

tr7 <- rbind(
  both_no_2plus_tr7,
  him_yes_le2_tr7, him_yes_2plus_tr7, her_yes_le2_tr7,
  her_yes_2plus_tr7, both_yes_le2_tr7, both_yes_2plus_tr7
)
round(tr7[,2:7], 2)

rng7 <- range(tr7$Lower, tr7$Upper)


rown <- c("BothNo_2Plus", "HimYes", "HimYes_2Plus",
          "HerYes", "HerYes_2Plus", "BothYes", "BothYes_2Plus")
rownames(tr1) <- rown
rownames(tr3) <- rown
rownames(tr5) <- rown
rownames(tr6) <- rown
rownames(tr7) <- rown


# ----  Forest Plots -----
col_newrel <- "#E69F00"
col_oldrel <- "#0072B2" 

# ---- Transitions to cohabitation and marriage -----
# Reproduces Figure 2
# Set margins
rng <- range(rng1, rng5)
rng[1] <- rng[1] - 0.1
rng[2] <- rng[2] + 0.1

svg("Results/Figure2.svg", width = 9.5, height = 6, pointsize = 13)
{
  par(
    mar = c(4, 5, 2, 1),
    oma = c(1, 4, 2, 2),
    cex.main = 2,
    cex.lab = 1.9,
    font.main = 1,
    xpd = TRUE
  )
  
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
         heights = c(3.5, 1.5)
  ) # 4/5 for plots, 1/5 for legend
  
  
  # Dating to cohabitation
  plot(1, 1,
       type = "n",
       xlim = rng,
       ylim = c(-0.5, 7.5),
       main = "Dating to cohabitation",
       xlab = "Hazard Ratio (HR)",
       ylab = "",
       yaxt = "n"
  )
  axis(2,
       at = c(.5, 2.5, 4.5, 6.5),
       labels = c("Neither\n Wants", "He\n Wants", "She\n Wants", "Both\n Want"),
       las = 2, cex.axis = 1.5
  )
  segments(1, -0.5, 1, 7.5, lwd = 1.5, lty = 2, col = "grey", cex.axis = 1.2)
  
  j <- -0.25
  for (i in which(rownames(tr1) %in% c("BothNo_2Plus", "HimYes_2Plus",
                                       "HerYes_2Plus", "BothYes_2Plus"))) { # more than 2 years
    
    segments(tr1$Lower[i], j + 1, tr1$Upper[i], j + 1, # horizontal bar
             col = col_oldrel, lwd = 2
    )
    segments(tr1$Lower[i], j + 1 - .15, tr1$Lower[i], j + 1 + .15, # vertical bars
             col = col_oldrel, lwd = 2
    )
    segments(tr1$Upper[i], j + 1 - .15, tr1$Upper[i], j + 1 + .15, # vertical bars
             col = col_oldrel, lwd = 2
    )
    points(tr1$HR[i], j + 1,
           pch = 15,
           col = col_oldrel
    )
    j <- j + 2
  }
  
  j <- 0.25
  for (i in which(rownames(tr1) %in% c("BothYes", "HimYes",
                                       "HerYes"))) { # up to 2 years
    
    segments(tr1$Lower[i], j + 2, tr1$Upper[i], j + 2, # horizontal bar
             col = col_newrel, lwd = 2
    )
    segments(tr1$Lower[i], j + 2 - .15, tr1$Lower[i], j + 2 + .15, # vertical bars
             col = col_newrel, lwd = 2
    )
    segments(tr1$Upper[i], j + 2 - .15, tr1$Upper[i], j + 2 + .15, # vertical bars
             col = col_newrel, lwd = 2
    )
    points(tr1$HR[i], j + 2,
           pch = 17,
           col = col_newrel
    )
    j <- j + 2
  }
  
  points(1, 0.25,
         pch = 17,
         col = col_newrel
  )
  
  # Cohabitation to marriage
  plot(1, 1,
       type = "n",
       xlim = rng,
       ylim = c(-0.5, 7.5),
       main = "Cohabitation to marriage",
       xlab = "Hazard Ratio (HR)",
       ylab = "",
       yaxt = "n"
  )
  axis(2,
       at = c(.5, 2.5, 4.5, 6.5),
       labels = c("Neither\n Wants", "He\n Wants", "She\n Wants", "Both\n Want"),
       las = 2, cex.axis = 1.5
  )
  segments(1, -0.5, 1, 7.5, lwd = 1.5, lty = 2, col = "grey", cex.axis = 1.2)
  
  j <- -0.25
  for (i in which(rownames(tr5) %in% c("BothNo_2Plus", "HimYes_2Plus",
                                       "HerYes_2Plus", "BothYes_2Plus"))) { # more than 2 years
    
    segments(tr5$Lower[i], j + 1, tr5$Upper[i], j + 1, # horizontal bar
             col = col_oldrel, lwd = 2
    )
    segments(tr5$Lower[i], j + 1 - .15, tr5$Lower[i], j + 1 + .15, # vertical bars
             col = col_oldrel, lwd = 2
    )
    segments(tr5$Upper[i], j + 1 - .15, tr5$Upper[i], j + 1 + .15, # vertical bars
             col = col_oldrel, lwd = 2
    )
    points(tr5$HR[i], j + 1,
           pch = 15,
           col = col_oldrel
    )
    j <- j + 2
  }
  
  j <- 0.25
  for (i in which(rownames(tr5) %in% c("BothYes", "HimYes",
                                       "HerYes"))) { # up to 2 years
    
    segments(tr5$Lower[i], j + 2, tr5$Upper[i], j + 2, # horizontal bar
             col = col_newrel, lwd = 2
    )
    segments(tr5$Lower[i], j + 2 - .15, tr5$Lower[i], j + 2 + .15, # vertical bars
             col = col_newrel, lwd = 2
    )
    segments(tr5$Upper[i], j + 2 - .15, tr5$Upper[i], j + 2 + .15, # vertical bars
             col = col_newrel, lwd = 2
    )
    points(tr5$HR[i], j + 2,
           pch = 17,
           col = col_newrel
    )
    j <- j + 2
  }
  
  points(1, 0.25,
         pch = 17,
         col = col_newrel
  )
  plot.new()
  legend("center",
         horiz = TRUE,
         legend = c("Up to 2 years", "More than 2 years"),
         cex = 1.8, pch = c(17, 15),
         col = c(col_newrel, col_oldrel), lwd = 2, bty = "n"
  )
}
dev.off()

# ---- Transitions to dissolution ----
# Reproduces Figure 3

rng <- range(rng3, rng6)
rng[1] <- rng[1] - 0.1
rng[2] <- rng[2] + 0.1

svg("Results/Figure3.svg", width = 9.5, height = 6, pointsize = 13)
{
  par(
    mar = c(5, 5, 4, 1),
    oma = c(1, 4, 2, 2),
    cex.main = 2,
    cex.lab = 1.9,
    font.main = 1,
    xpd = TRUE
  )
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
         heights = c(3.5, 1.5)
  ) # 4/5 for plots, 1/5 for legend
  
  # Dating to dissolution
  plot(1, 1,
       type = "n",
       xlim = rng,
       ylim = c(-0.5, 7.5),
       main = "Dating to dissolution",
       xlab = "Hazard Ratio (HR)",
       ylab = "",
       yaxt = "n"
  )
  axis(2,
       at = c(.5, 2.5, 4.5, 6.5),
       labels = c("Neither\n Wants", "He\n Wants", "She\n Wants", "Both\n Want"),
       las = 2, cex.axis = 1.5
  )
  segments(1, -0.5, 1, 7.5, lwd = 1.5, lty = 2, col = "grey", cex.axis = 1.2)
  j <- -0.25
  for (i in which(rownames(tr3) %in% c("BothNo_2Plus", "HimYes_2Plus",
                                       "HerYes_2Plus", "BothYes_2Plus"))) { # more than 2 years
    
    segments(tr3$Lower[i], j + 1, tr3$Upper[i], j + 1, # horizontal bar
             col = col_oldrel, lwd = 2
    )
    segments(tr3$Lower[i], j + 1 - .15, tr3$Lower[i], j + 1 + .15, # vertical bars
             col = col_oldrel, lwd = 2
    )
    segments(tr3$Upper[i], j + 1 - .15, tr3$Upper[i], j + 1 + .15, # vertical bars
             col = col_oldrel, lwd = 2
    )
    points(tr3$HR[i], j + 1,
           pch = 15,
           col = col_oldrel
    )
    j <- j + 2
  }
  
  j <- 0.25
  for (i in which(rownames(tr3) %in% c("BothYes", "HimYes",
                                       "HerYes"))) { # up to 2 years
    
    segments(tr3$Lower[i], j + 2, tr3$Upper[i], j + 2, # horizontal bar
             col = col_newrel, lwd = 2
    )
    segments(tr3$Lower[i], j + 2 - .15, tr3$Lower[i], j + 2 + .15, # vertical bars
             col = col_newrel, lwd = 2
    )
    segments(tr3$Upper[i], j + 2 - .15, tr3$Upper[i], j + 2 + .15, # vertical bars
             col = col_newrel, lwd = 2
    )
    points(tr3$HR[i], j + 2,
           pch = 17,
           col = col_newrel
    )
    j <- j + 2
  }
  
  points(1, 0.25,
         pch = 17,
         col = col_newrel
  )
  
  # Cohabitation/marriage to dissolution
  plot(1, 1,
       type = "n",
       xlim = rng,
       ylim = c(-0.5, 7.5),
       main = "Cohabitation/marriage\n to dissolution",
       xlab = "Hazard Ratio (HR)",
       ylab = "",
       yaxt = "n"
  )
  axis(2,
       at = c(.5, 2.5, 4.5, 6.5),
       labels = c("Neither\n Wants", "He\n Wants", "She\n Wants", "Both\n Want"),
       las = 2, cex.axis = 1.5
  )
  segments(1, -0.5, 1, 7.5, lwd = 1.5, lty = 2, col = "grey", cex.axis = 1.2)
  
  j <- -0.25
  for (i in which(rownames(tr6) %in% c("BothNo_2Plus", "HimYes_2Plus",
                                       "HerYes_2Plus", "BothYes_2Plus"))) { # more than 2 years
    
    segments(tr6$Lower[i], j + 1, tr6$Upper[i], j + 1, # horizontal bar
             col = col_oldrel, lwd = 2
    )
    segments(tr6$Lower[i], j + 1 - .15, tr6$Lower[i], j + 1 + .15, # vertical bars
             col = col_oldrel, lwd = 2
    )
    segments(tr6$Upper[i], j + 1 - .15, tr6$Upper[i], j + 1 + .15, # vertical bars
             col = col_oldrel, lwd = 2
    )
    points(tr6$HR[i], j + 1,
           pch = 15,
           col = col_oldrel
    )
    j <- j + 2
  }
  
  j <- 0.25
  for (i in which(rownames(tr6) %in% c("BothYes", "HimYes",
                                       "HerYes"))) { # up to 2 years
    
    segments(tr6$Lower[i], j + 2, tr6$Upper[i], j + 2, # horizontal bar
             col = col_newrel, lwd = 2
    )
    segments(tr6$Lower[i], j + 2 - .15, tr6$Lower[i], j + 2 + .15, # vertical bars
             col = col_newrel, lwd = 2
    )
    segments(tr6$Upper[i], j + 2 - .15, tr6$Upper[i], j + 2 + .15, # vertical bars
             col = col_newrel, lwd = 2
    )
    points(tr6$HR[i], j + 2,
           pch = 17,
           col = col_newrel
    )
    j <- j + 2
  }
  
  points(1, 0.25,
         pch = 17,
         col = col_newrel
  )
  
  plot.new()
  legend("center",
         horiz = TRUE,
         legend = c("Up to 2 years", "More than 2 years"),
         cex = 1.8, pch = c(17, 15),
         col = c(col_newrel, col_oldrel), lwd = 2, bty = "n"
  )
}
dev.off()

# ---- Cohabitation/Marriage to parenthood ----
# Reproduces Figure 2 Appendix
svg("Results/Figure2Appendix.svg", width = 9, height = 7, pointsize = 13)

{
  par(
    mar = c(7, 4, 2, 1),
    oma = c(3, 4, 1, 2),
    cex.main = 2,
    cex.lab = 1.9,
    font.main = 1,
    xpd = TRUE
  )
  
  plot(1, 1,
     type = "n",
     xlim = rng7,
     ylim = c(-0.5, 7.5),
     main = "Cohabitation/marriage to parenthood",
     xlab = "Hazard Ratio (HR)",
     ylab = "",
     yaxt = "n"
)
axis(2,
     at = c(.5, 2.5, 4.5, 6.5),
     labels = c("Neither\n Wants", "He\n Wants", "She\n Wants", "Both\n Want"),
     las = 2, cex.axis = 1.5
)
segments(1, -0.5, 1, 7.5, lwd = 1.5, lty = 2, col = "grey", cex.axis = 1.2)

j <- -0.25
for (i in which(rownames(tr7) %in% c("BothNo_2Plus", "HimYes_2Plus",
                                     "HerYes_2Plus", "BothYes_2Plus"))) { # more than 2 years
  
  segments(tr7$Lower[i], j + 1, tr7$Upper[i], j + 1, # horizontal bar
           col = col_oldrel, lwd = 2
  )
  segments(tr7$Lower[i], j + 1 - .15, tr7$Lower[i], j + 1 + .15, # vertical bars
           col = col_oldrel, lwd = 2
  )
  segments(tr7$Upper[i], j + 1 - .15, tr7$Upper[i], j + 1 + .15, # vertical bars
           col = col_oldrel, lwd = 2
  )
  points(tr7$HR[i], j + 1,
         pch = 15,
         col = col_oldrel
  )
  j <- j + 2
}

j <- 0.25
for (i in which(rownames(tr7) %in% c("BothYes", "HimYes",
                                     "HerYes"))) { # up to 2 years
  
  segments(tr7$Lower[i], j + 2, tr7$Upper[i], j + 2, # horizontal bar
           col = col_newrel, lwd = 2
  )
  segments(tr7$Lower[i], j + 2 - .15, tr7$Lower[i], j + 2 + .15, # vertical bars
           col = col_newrel, lwd = 2
  )
  segments(tr7$Upper[i], j + 2 - .15, tr7$Upper[i], j + 2 + .15, # vertical bars
           col = col_newrel, lwd = 2
  )
  points(tr7$HR[i], j + 2,
         pch = 17,
         col = col_newrel
  )
  j <- j + 2
}

points(1, 0.25,
       pch = 17,
       col = col_newrel
)
legend("bottom",
       inset = c(0, -.4),
       horiz = TRUE,
       xpd = TRUE,
       legend = c("Up to 2 years", "More than 2 years"),
       cex = 1.8, pch = c(17, 15),
       col = c(col_newrel, col_oldrel), lwd = 2, bty = "n"
)
}
dev.off()
# --- Age effects ----
secol <- "#009E73"
agecol <- "#56B4E9"

svg("Results/Figure3Appendix.svg", width = 11, height = 10, pointsize = 13)
{
  par(
    mar = c(5, 5, 4, 1),
    oma = c(1, 4, 2, 2),
    cex.main = 2,
    cex.lab = 1.9,
    font.main = 1,
    xpd = TRUE
  )
  
  layout(matrix(c(1, 1, 2, 2, 3, 3, 0, 4, 4, 5, 5, 0), nrow = 2, byrow = TRUE)
  ) # 4/5 for plots, 1/5 for legend
  
  termplot(bestmod, terms = 38, se = T,
         xlab = "Age at baseline", 
         ylab = "Partial effect",
         main = "Dating to cohabitation",
         lwd.term = 2,
         col.term = agecol, col.se = secol)
termplot(bestmod, terms = 39, se = T,
         xlab = "Age at baseline", 
         ylab = "Partial effect",
         main = "Dating to dissolution",
         lwd.term = 2,
         col.term = agecol, col.se = secol,)
termplot(bestmod, terms = 40, se = T,
         xlab = "Age at baseline", 
         ylab = "Partial effect",
         main = "Cohabitation to marriage",
         lwd.term = 2,
         col.term = agecol, col.se = secol)
termplot(bestmod, terms = 41, se = T,
         xlab = "Age at baseline", 
         ylab = "Partial effect",
         main = "Cohabitation/marriage \n to dissolution",
         lwd.term = 2,
         col.term = agecol, col.se = secol)
termplot(bestmod, terms = 42, se = T,
         xlab = "Age at baseline", 
         ylab = "Partial effect",
         main = "Cohabitation/marriage \n to parenthood",
         lwd.term = 2,
         col.term = agecol, col.se = secol)
}
dev.off()
