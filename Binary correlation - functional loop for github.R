library('tidyverse')
library("lme4")
library('gee')
library("openxlsx")

#multinomial function - define before sims
simulate_binary_pair <- function(n, p1, p2, rho) {
  P11 <- p1*p2 + rho*sqrt(p1*(1-p1)*p2*(1-p2))
  probs <- c(
    P11,
    p1 - P11,
    p2 - P11,
    1 - p1 - p2 + P11
  )
  outcomes <- sample(1:4, size=n, replace=TRUE, prob=probs)
  RE <- ifelse(outcomes %in% c(1,2), 1, 0)
  LE <- ifelse(outcomes %in% c(1,3), 1, 0)
  cbind(RE, LE)
}


loopsize <- 10000
  
set.seed(16)

for (O in 0:2) {
  
  #sets estimand
  ORsim <- 1 + 0.25*O

for(i in 0:6) {
  
  #sets correlation from 0 to 0.9 in 0.15 steps
  r <- 0.15*i
  
for(ss in 1:2) {
  
  #sets sample size
  n <- 200*ss
  
  #sets empty lists for data export
  Rightdf <- list()
  Subjdf <- list()
  Bothdf <- list()
  GLMMdf <- list()
  GEEdf <- list()
  
  #counter for loops
  comp <- 0
  
  #counter for attempted loops
  attempts <- 0
  
while (comp < loopsize) {
  
  attempts <- attempts + 1
  
  #set base parameters
  ID  <- 1:n
  Expo <- rbinom(n, 1, 0.5)
  alpha <- 0.05
  z <- qnorm(1 - alpha/2)
  
  #combine into data frame
  df <- data.frame(
    ID, Expo,
    OutcR = NA_integer_,
    OutcL = NA_integer_,
    OutcT = NA_integer_
  )
  
  #simulate Bernoulli outcome for each subject, varying outcome p by exposure
  for (id in unique(df$ID)) {
    
    dfloop <- df[df$ID == id,]
    
    Outcs <- if(df$Expo[df$ID ==id] == 1) {
      simulate_binary_pair(1, ORsim/(ORsim+1), ORsim/(ORsim+1), r)
    } else {
      simulate_binary_pair(1, 0.5, 0.5, r)
    }
    
    df$OutcR[df$ID == id] <- Outcs[1]
    df$OutcL[df$ID == id] <- Outcs[2] 
  }
  
  df$OutcT <- ifelse(df$OutcR == 0 & df$OutcL ==0, 0, 1)
  
  #pivots dataframe to long format
  dflong <- df %>%
    pivot_longer(cols = c(OutcR, OutcL),
                 names_to = "Eye",
                 names_pattern = "Outc(R|L)",
                 values_to = "Outcome")
  
  TR <- table(df$Expo, df$OutcR)
  TL <- table(df$Expo, df$OutcL)
  
  if(isTRUE(all.equal(TR, TL))) {
    next
  }
  
  comp <- comp + 1
  
  #Right eye's logistic regression
  LogOne <- glm(OutcR ~ Expo, family = binomial, data=df)
  LgMOne <- summary(LogOne)
  
  CIOne <- cbind(
    OR <- exp(coef(LogOne)["Expo"]),
    SE <- LgMOne$coefficients["Expo","Std. Error"],
    LCL <- exp(coef(LogOne)["Expo"] - z*SE),
    UCL <- exp(coef(LogOne)["Expo"] + z*SE)
  )
  
  #Per subject logistic regression
  LogSubj <- glm(OutcT ~ Expo, family = binomial, data=df)
  LgMSubj <- summary(LogSubj)
  
  CISubj <- cbind(
    OR <- exp(coef(LogSubj)["Expo"]),
    SE <- LgMSubj$coefficients["Expo","Std. Error"],
    LCL <- exp(coef(LogSubj)["Expo"] - z*SE),
    UCL <- exp(coef(LogSubj)["Expo"] + z*SE)
  )
  
  #Both eyes' logistic regression
  LogBoth <- glm(Outcome ~ Expo, family = binomial, data=dflong)
  LgMBoth <- summary(LogBoth)
  
  CIBoth <- cbind(
    OR <- exp(coef(LogBoth)["Expo"]),
    SE <- LgMBoth$coefficients["Expo","Std. Error"],
    LCL <- exp(coef(LogBoth)["Expo"] - z*SE),
    UCL <- exp(coef(LogBoth)["Expo"] + z*SE)
  )
  
  #GLMM model (logistic mixed-effects model) and extraction
  GLMM <- glmer(Outcome ~ Expo + (1 | ID),
                data    = dflong,
                family  = binomial(link = "logit")
  )
  ORGLMMs <- summary(GLMM)$coefficients
  
  CIGLMM <- cbind(
    Estimate = ORGLMMs[, "Estimate"],
    SE       = ORGLMMs[, "Std. Error"],
    Lower    = ORGLMMs[, "Estimate"] - z * ORGLMMs[, "Std. Error"],
    Upper    = ORGLMMs[, "Estimate"] + z * ORGLMMs[, "Std. Error"]
  )
  
  CIGLMM2 <- within(as.data.frame(CIGLMM), {
    OR    <- exp(Estimate)
    LCL   <- exp(Lower)
    UCL   <- exp(Upper)
  })
  CIGLMM2[2,c(5,6,7)]
  
  #GEE model and extraction
  GEE <- try(suppressMessages(
    gee(
      Outcome ~ Expo,
      id      = dflong$ID,
      data    = dflong,
      family  = binomial,
      corstr  = "exchangeable"
    )), silent=TRUE)
  ORGEEs <- summary(GEE)$coefficients
  
  CIGEE <- cbind(
    Estimate = ORGEEs[, "Estimate"],
    SE       = ORGEEs[, "Robust S.E."],
    Lower    = ORGEEs[, "Estimate"] - z * ORGEEs[, "Robust S.E."],
    Upper    = ORGEEs[, "Estimate"] + z * ORGEEs[, "Robust S.E."]
  )
  
  CIGEE2 <- within(as.data.frame(CIGEE), {
    OR  <- exp(Estimate)
    LCL <- exp(Lower)
    UCL <- exp(Upper)
  })
  
  #extract modeled values into dataframes
  
  #Extraction from right eye model
  Rightdf[[comp]] <- data.frame(
    Sim = comp,
    pvalue = LgMOne$coefficients["Expo", "Pr(>|z|)"],
    Significant = ifelse(LgMOne$coefficients["Expo", "Pr(>|z|)"] < 0.05, 1, 0),
    OR = CIOne[[1]],
    LCL = CIOne[[3]],
    UCL = CIOne[[4]],
    LogOR = coef(LogOne)["Expo"],
    StdError = CIOne[[2]]
  )

  #Extraction from per subject model
  Subjdf[[comp]] <- data.frame(
    Sim = comp,
    pvalue = LgMSubj$coefficients["Expo", "Pr(>|z|)"],
    Significant = ifelse(LgMSubj$coefficients["Expo", "Pr(>|z|)"] < 0.05, 1, 0),
    OR = CISubj[[1]],
    LCL = CISubj[[3]],
    UCL = CISubj[[4]],
    LogOR = coef(LogSubj)["Expo"],
    StdError = CISubj[[2]]
  ) 
  
  #Extraction from both eyes model
  Bothdf[[comp]] <- data.frame(
    Sim = comp,
    pvalue = LgMBoth$coefficients["Expo", "Pr(>|z|)"],
    Significant = ifelse(LgMBoth$coefficients["Expo", "Pr(>|z|)"] < 0.05, 1, 0),
    OR = CIBoth[[1]],
    LCL = CIBoth[[3]],
    UCL = CIBoth[[4]],
    LogOR = coef(LogBoth)["Expo"],
    StdError = CIBoth[[2]]
  )
  
  #Extraction from GLMM
  GLMMdf[[comp]] <- data.frame(
    Sim = comp,
    pvalue = ORGLMMs[[2,4]],
    Significant = ifelse(ORGLMMs[[2,4]] < 0.05, 1, 0),
    OR = CIGLMM2[[2,7]],
    LCL = CIGLMM2[[2,6]],
    UCL = CIGLMM2[[2,5]],
    LogOR = ORGLMMs[[2,1]],
    StdError = ORGLMMs[[2,2]],
    Attempts = attempts
  )
  
  #Extraction from GEE
  GEEdf[[comp]] <- data.frame(
    Sim = comp,
    pvalue = 2 * pnorm(-abs(ORGEEs[[2,5]])),
    Significant = ifelse(2 * pnorm(-abs(ORGEEs[[2,5]])) < 0.05, 1, 0),
    OR = CIGEE2[[2,7]],
    LCL = CIGEE2[[2,6]],
    UCL = CIGEE2[[2,5]],
    LogOR = ORGEEs[[2,1]],
    StdError = ORGEEs[[2,4]]
  )
}
  
  Right_results <- do.call(rbind, Rightdf)
  Subj_results <- do.call(rbind, Subjdf)
  Both_results <- do.call(rbind, Bothdf)
  GLMM_results <- do.call(rbind, GLMMdf)
  GEE_results <- do.call(rbind, GEEdf)
  
  Right_summary <- data.frame(Sim = "Mean",
                             pvalue = NA,
                             Significant = mean(Right_results$Significant, na.rm = TRUE),
                             OR = NA,
                             LCL = NA,
                             UCL = NA,
                             LogOR = mean(Right_results$LogOR, na.rm = TRUE),
                             StdError = mean(Right_results$StdError, na.rm = TRUE)
                             )
  
  Subj_summary <- data.frame(Sim = "Mean",
                              pvalue = NA,
                              Significant = mean(Subj_results$Significant, na.rm = TRUE),
                              OR = NA,
                              LCL = NA,
                              UCL = NA,
                              LogOR = mean(Subj_results$LogOR, na.rm = TRUE),
                              StdError = mean(Subj_results$StdError, na.rm = TRUE)
                              )
  
  Both_summary <- data.frame(Sim = "Mean",
                              pvalue = NA,
                              Significant = mean(Both_results$Significant, na.rm = TRUE),
                              OR = NA,
                              LCL = NA,
                              UCL = NA,
                              LogOR = mean(Both_results$LogOR, na.rm = TRUE),
                              StdError = mean(Both_results$StdError, na.rm = TRUE)
                              )
  
  GLMM_summary <- data.frame(Sim = "Mean",
                             pvalue = NA,
                             Significant = mean(GLMM_results$Significant, na.rm = TRUE),
                             OR = NA,
                             LCL = NA,
                             UCL = NA,
                             LogOR = mean(GLMM_results$LogOR, na.rm = TRUE),
                             StdError = mean(GLMM_results$StdError, na.rm = TRUE),
                             Attempts = NA
                            )
  
  GEE_summary <- data.frame(Sim = "Mean",
                             pvalue = NA,
                             Significant = mean(GEE_results$Significant, na.rm = TRUE),
                             OR = NA,
                             LCL = NA,
                             UCL = NA,
                             LogOR = mean(GEE_results$LogOR, na.rm = TRUE),
                             StdError = mean(GEE_results$StdError, na.rm = TRUE)
  )
  
  Right_results <- rbind(Right_results, Right_summary)
  Subj_results <- rbind(Subj_results, Subj_summary)
  Both_results <- rbind(Both_results, Both_summary)
  GLMM_results <- rbind(GLMM_results, GLMM_summary)
  GEE_results <- rbind(GEE_results, GEE_summary)
  
  filenameRight <- paste0("Right_n", n, "_r", r, "_OR", ORsim, ".xlsx")
  filenameSubj <- paste0("Subj_n", n, "_r", r, "_OR", ORsim, ".xlsx")
  filenameBoth <- paste0("Both_n", n, "_r", r, "_OR", ORsim, ".xlsx")
  filenameGLMM <- paste0("GLMM_n", n, "_r", r, "_OR", ORsim, ".xlsx")
  filenameGEE <- paste0("GEE_n", n, "_r", r, "_OR", ORsim, ".xlsx")
  
  write.xlsx(Right_results, filenameRight, rownames = FALSE)
  write.xlsx(Subj_results, filenameSubj, rownames = FALSE)
  write.xlsx(Both_results, filenameBoth, rownames = FALSE)
  write.xlsx(GLMM_results, filenameGLMM, rownames = FALSE)
  write.xlsx(GEE_results, filenameGEE, rowNames = FALSE)
  
}
  
}
  
}
