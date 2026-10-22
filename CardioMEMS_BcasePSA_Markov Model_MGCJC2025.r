#clear my environment
rm(list=ls())

#clear my R console
cat("\014")

# Libraries
library(tidyverse)
library(twig)
library(MASS)
library(matrixStats)

#Definition for key acronyms
#Imp_comp - implant-related complication
#pHosp_StableHFPostDischarge30D_cond - probability of hospitalization given survival, in the stableHFpostdischarge30D state
#pHosp_StableHF_cond - probability of hospitalization given survival, in the stableHF state



DiscRate <-	0.015 

n_samples <- 5000 
set.seed(12345)

# log-normal distribution ===========
rlnorm_disutility <- function(n_samples, util_mean, util_sd, N){ # generate random draws of disutilities using a log-normal distribution
  disutil_mean <- 1 - util_mean 
  util_se <- util_sd / sqrt(N)  
  adjustment_factor <- log((util_se/disutil_mean)^2+1) 
  disutil_mean_adj <- log(disutil_mean) - 0.5 * adjustment_factor 
  disutil_se_adj <- sqrt(adjustment_factor) 
  1 - rlnorm(n = n_samples, meanlog = disutil_mean_adj, sdlog = disutil_se_adj)
}

# correlated log-normal distributin
mvrlnorm_disutility <- function(n_samples, util_mean, util_sd, util_rho, N){
  disutil_mean <- 1 - util_mean
  util_se <- util_sd / sqrt(N)
  adjustment_factor <- log((util_se/disutil_mean)^2+1)
  disutil_mean_adj <- log(disutil_mean) - 0.5 * adjustment_factor
  disutil_se_adj <- sqrt(adjustment_factor)
  disutil_Sigma_adj <- disutil_se_adj %*% t(disutil_se_adj) * util_rho # covariance matrix, which combines the SD of the adjusted disutilities with correlation matrix util_rho
  1 - exp(mvrnorm(n = n_samples, mu = disutil_mean_adj, Sigma = disutil_Sigma_adj)) # sigma = covaraince matrix, mu= mean vector (mean of log disutilities) that describe the mutlivariate distribution we are sampling from
}
# Gamma distribution ==========
rgamma_mean_sd <- function(n_samples, mean_x, sd_x, N){
  se_x <- sd_x / sqrt(N)
  shape <- (mean_x / se_x)^2
  scale <- se_x^2 / mean_x
  rgamma(n = n_samples, shape = shape, scale = scale)
}

# beta distribution =========
rbeta_mean_sd <- function(n_samples, mean_x, N){
  rbeta(n = n_samples, shape1 = mean_x * N, shape2 = N * (1 - mean_x))
}
# Shape 1= apha, and shape 2= beta
# lognormal risks ===========
rlnorm_ci <- function(n_samples, mean_x, lb_x, ub_x){
  rlnorm(n_samples, meanlog = log(mean_x), 
         sdlog = (log(ub_x) - log(lb_x)) / (1.96 * 2))
}


mytwig <- twig() +
  decisions(names=c(UsualCare,CardioMEMS)) +
  states(names=c(StableHFPostDischarge30D,StableHF,Dead), 
         init_probs=c(1,0,0), 
         max_cycles=c(1,1,1)) + 
  event(name=die, 
        options=c(yes,none), 
        probs=c(pDie,leftover), 
        transitions=c(Dead,Hosp))+
  event(name=Hosp, 
        options=c(yes,none), 
        probs=c(pHosp,leftover), 
        transitions=c(StableHFPostDischarge30D,ProgStableHF))+
  event(name=ProgStableHF, 
        options=c(yes,none), 
        probs=c(pProgStableHF,leftover), 
        transitions=c(StableHF,stay))+
  payoffs(names=c(cost, utility),
          discount_rates=c(DiscRate,DiscRate))

# Probability of Death
## Return the probability of death based on the current health state
pDie <- function(state, pDie_StableHFPostDischarge30D, pDie_StableHF){
  ## Baseline monthly mortality probablity based on state
  base <- ifelse(state == "StableHFPostDischarge30D", pDie_StableHFPostDischarge30D,
                 ifelse(state == "StableHF", pDie_StableHF,0))
}

# Probability of Hospitalization, define how pHosp be calculated
pHosp <- function(state, decision, pHosp_StableHFPostDischarge30D_cond, pHosp_StableHF_cond, HRatio_Hosp){
  ## Baseline hospitalziation probablity based on state
  base <- ifelse(state == "StableHFPostDischarge30D", pHosp_StableHFPostDischarge30D_cond,
                 ifelse(state == "StableHF", pHosp_StableHF_cond,0))
  ## apply CardioMEMS effect
  ifelse(decision == "CardioMEMS",
         base * HRatio_Hosp, base)
}

# Probability of Progression to Stable HF
## Return the probability of progression to Stable HF based on health state and decision
pProgStableHF <- function(state){
  ## Progression to Stable HF is only possible from the StableHFPostDischarge30D
  ifelse(state =="StableHFPostDischarge30D", 1, 0) # 1 and 0 are logical indicators that determine whether 
  #progression to stableHF is possible from a given state. If the patient is in "StableHFPostDischarge30D", 
  #progression to StableHF is allowed, 
  #and the actual probability is the residual probability after accounting for death and hosp. 
  #For all other states, progression is not possible, so the function returns 0
}

# Cost Payoff
cost <- function(state, decision, cycle,
                 cOutpatientMonthly, cCardioMEMSDevice, cCardioMEMSImplant, cCardioMEMSMonitoringMonthly,
                 cHospEvent, cImp_comp, pImp_comp, pHosp_StableHFPostDischarge30D_cond,pHosp_StableHF_cond,HRatio_Hosp,cDead){
  #Monthly State Cost
  state_cost <- ifelse(state == "StableHFPostDischarge30D", cOutpatientMonthly,
                       ifelse (state == "StableHF", cOutpatientMonthly, cDead))
  
  # Monthly CardioMEMS monitoring cost, CardioMEMS arm only
  cMonitoring <- ifelse(decision =="CardioMEMS" & state != "Dead", cCardioMEMSMonitoringMonthly,
                        0)
  
  # Total CardioMEMS-related upfront cost (one-time cost only, for the first cycle), CardioMEMS arm only
  Upfront_cost <- ifelse(decision == "CardioMEMS" & cycle==1,
                         cCardioMEMSDevice + cCardioMEMSImplant,
                         0)
  #Expected hospitalziation cost, a cycle
  Hosp_cost <- pHosp(state, decision, pHosp_StableHFPostDischarge30D_cond,pHosp_StableHF_cond,HRatio_Hosp) * cHospEvent
  
  #Expected cost of complication, first cycle, CardioMEMS only
  Complication_cost <- ifelse(state == "StableHFPostDischarge30D" &
                                decision == "CardioMEMS" &
                                cycle == 1,
                                pImp_comp * cImp_comp,
                                0)
                              
  # Return the per-cycle cost for the current state and decision
  return(state_cost + cMonitoring + Upfront_cost + Hosp_cost + Complication_cost)
}

# Utility Payoff
# Utility adjustment if hospitalized. # Sample uHosp as utility, but not treating hospitalzation as a state utility
##Instead, use it to compute a temporary utility decrement when a hospitalzation occurs during the cycle.

#constants, not arguments
los_days <- 7.5
cycle_days <- 30
los_comp_days <- 4.15 

uHospAdj <- function(state, 
                     uStableHFPostDischarge30D, uStableHF, uHosp) {
  uState <- ifelse(state == "StableHFPostDischarge30D", uStableHFPostDischarge30D,
                   ifelse(state == "StableHF", uStableHF, 0))
  (uHosp-uState) * (los_days/cycle_days)
}

# Utility adjustment if implant-related complications, not tretaing it as a state utility
## Instead use it to compute temporary utility decrement when complications happen in the first cycle. This will only be applied for the first cycle

ucompAdj <- function(state, 
                     uStableHFPostDischarge30D, uStableHFPostDischarge30D_comp) {
            ifelse(state=="StableHFPostDischarge30D", (uStableHFPostDischarge30D_comp - uStableHFPostDischarge30D) * (los_comp_days /cycle_days), 0)
}  
  
utility <- function(state,cycle,decision,
                    uStableHFPostDischarge30D,uStableHFPostDischarge30D_comp, uStableHF, uDead, uHosp,
                    pImp_comp, pHosp_StableHFPostDischarge30D_cond,
                    pHosp_StableHF_cond,
                    HRatio_Hosp
                      ){
  
  ## Baseline utility weight for the state
  uState <- ifelse(state == "StableHFPostDischarge30D", uStableHFPostDischarge30D,
                   ifelse(state == "StableHF", uStableHF, uDead))
  
 
  ## utility adjustment if hospitalized, it represents the temporary utility drop casued by hosp during that cycle
  du_Hosp <- pHosp(state, decision, pHosp_StableHFPostDischarge30D_cond, pHosp_StableHF_cond, HRatio_Hosp) *
             uHospAdj(state, 
                      uStableHFPostDischarge30D,
                      uStableHF,
                      uHosp)
                      
  
  ## utility adjustment if complications, it represents the temporary utility drop casued by complications during the first cycle, only for CardioMEMS arm
  du_comp <- ifelse(decision == "CardioMEMS" & cycle == 1,
                    pImp_comp * ucompAdj(state,
                                         uStableHFPostDischarge30D,
                                         uStableHFPostDischarge30D_comp),
                         0)
  
  #Expected utility weight, per cycle, also ensure that utility weights are monthly, even though the inputs were yearly weights
  return((uState + du_Hosp + du_comp) * (cycle_days / 365))
}   

# correlated log-normal utilities ===============
# first create the covariance matrix from the correlation (rho) and the utility sd
mat_utilities <- mvrlnorm_disutility(n_samples = n_samples, 
                                     util_mean = c(StableHFPostDischarge30D = 0.790, StableHFPostDischarge30D_comp = 0.743, StableHF = 0.780, Hosp = 0.620), 
                                     util_sd = c(StableHFPostDischarge30D = 0.200, StableHFPostDischarge30D_comp = 0.203, StableHF = 0.200, Hosp = 0.250),
                                     util_rho =  matrix(c(
                                       1.0, 0.8, 0.8, 0.5,
                                       0.8, 1.0, 0.8, 0.5,
                                       0.8, 0.8, 1.0, 0.5,
                                       0.5, 0.5, 0.5, 1.0
                                     ), nrow =4, byrow = TRUE), 
                                     
                                     N = c(StableHFPostDischarge30D = 6943, StableHFPostDischarge30D_comp = 6943, StableHF = 6943, Hosp = 6943)
 
                                                                         )
# check the joint distribution of utilities 
plot(mat_utilities[,"StableHFPostDischarge30D"], mat_utilities[,"StableHF"])
plot(mat_utilities[,"StableHFPostDischarge30D"], mat_utilities[,"Hosp"])
plot(mat_utilities[,"StableHF"], mat_utilities[,"Hosp"])
plot(mat_utilities[,"StableHFPostDischarge30D"], mat_utilities[,"StableHFPostDischarge30D_comp"])
plot(mat_utilities[,"StableHF"], mat_utilities[,"StableHFPostDischarge30D_comp"])
plot(mat_utilities[,"Hosp"], mat_utilities[,"StableHFPostDischarge30D_comp"])
cor(mat_utilities)
colMeans(mat_utilities)
colSds(mat_utilities)



# Parameters List

params_psa = data.frame(uStableHFPostDischarge30D	= mat_utilities [, "StableHFPostDischarge30D"],
                        uStableHFPostDischarge30D_comp = mat_utilities [, "StableHFPostDischarge30D_comp"],
                        uStableHF	= mat_utilities [, "StableHF"],
                        uHosp = mat_utilities [, "Hosp"],
                        uDead	= 0,
                        
                        cOutpatientMonthly	= 100.22,
                        cHospEvent	= 10753.56,
                        cCardioMEMSDevice = 18860.22,
                        cCardioMEMSImplant = 15561.62,
                        cImp_comp = rgamma_mean_sd(n_samples = n_samples,
                                                   mean_x = 8808.44,
                                                   sd_x = 2235.75,
                                                   N = 1),
                        
                        cCardioMEMSMonitoringMonthly = 131.52,
                        cDead	= 0,
                        
                        #probability of complication
                        pImp_comp	= rbeta_mean_sd(
                                 n_samples = n_samples,
                                 mean_x = 0.026,
                                 N = 575
                        ),
                        
                        #Mortality Probabilities
                        pDie_StableHFPostDischarge30D	= 0.087,
                        pDie_StableHF	= rbeta_mean_sd(n_samples = n_samples, mean_x = 0.0108, N = 955),
                        
                        #Hospitalization Probabilities, among survivors
                        pHosp_StableHFPostDischarge30D_cond	= 0.115,
                        pHosp_StableHF_cond =	rbeta_mean_sd(n_samples = n_samples, mean_x = 0.0437, N = 955),
                         
                        #Intervention Effect
                        HRatio_Hosp = rlnorm_ci(n_samples = n_samples, mean_x = 0.70, lb_x = 0.58, ub_x = 0.86)
)
head(params_psa)


# Run Model (PSA)
results_PSA <- run_twig(twig_obj = mytwig, params = params_psa, n_cycles = 360, parallel=TRUE)

#ev means expected values
# results for the mean of the simulations, and results for each simulations can be found using the codes below:
results_PSA$mean_ev 
results_PSA$sim_ev
results_PSA$sim_ev[,,"5000"]


# All simulated costs for CardioMEMS
results_PSA$sim_ev["CardioMEMS", "cost",]
# All simulated utilities for CardioMEMS
results_PSA$sim_ev["CardioMEMS", "utility",]

str(results_PSA)

# Calculate ICERs
calculate_icers(results_PSA$mean_ev)

# Plot CEAC 
## WTP range is from 0-1,000,000; varying by 10000: which means 0, 10000, 20000, 30000, etc
plot_ceac(results_PSA$sim_ev, wtp_range = seq(0, 1000000, by = 10000))

# Scatterplot


## create dataframe for ModEx
df <- as.data.frame(as.table(results_PSA$sim_ev)) %>%
  setNames(c("decision", "payoff", "sim", "value")) %>%
  pivot_wider(names_from = payoff, values_from = value)

head(df)
# Compute incremental outcomes vs. Usual care
## Within each simulation, calculate df cost, differential
df_inc <- df %>%
  group_by(sim) %>%
  mutate(
    cost_UC    = cost[decision == "UsualCare"],
    util_UC    = utility[decision == "UsualCare"],
    dCost        = cost - cost_UC,
    dUtility     = utility - util_UC
  ) %>%
  ungroup() %>%
  filter(decision != "UsualCare") %>%
  mutate(comparison = "CardioMEMS vs UsualCare")

ggplot(df_inc, aes(x = dUtility, y = dCost, color = comparison)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Incremental Utility",
    y = "Incremental Cost",
    color = "Comparison",
    title = "Cost-Effectiveness Plane"
  ) +
  theme_minimal()



# Incremental Net Health Benefit Sensitivity Analysis

## Create a regression-ready dataset from the PSA outputs.
## Each row represents one PSA simulation and contains the costs and QALYs for both usual care and CardioMEMS

df_iNHB <- df %>%
  mutate(
    sim = as.integer(as.character(sim))
  )%>%
  dplyr::select(sim, decision, cost, utility) %>%
  tidyr::pivot_wider(
    names_from = decision,
    values_from = c(cost, utility)
  )%>%
  dplyr::rename(
    c_UsualCare = cost_UsualCare,
    c_CardioMEMS = cost_CardioMEMS,
    e_UsualCare = utility_UsualCare,
    e_CardioMEMS = utility_CardioMEMS
  )

  ## Add the parameter values sampled in each PSA simulation.
  ## The "p_" prefix distinguishes model paramteres from model outcomes.

params_psa_iNHB <- params_psa %>%
  mutate(sim = 1:n()) %>%
  dplyr::rename_with(~paste0("p_", .), -sim) 
  
  df_iNHB <- df_iNHB %>%
    left_join(params_psa_iNHB, by = "sim") %>%
  dplyr::select(-sim)

## calculate incremental net benefit (iNHB) for CardioMEMS
## compared with usual care at a WTP threshold of $50,000 per QALY.

## iNHB = incremental QALYs - (incremental costs/WTP)

## positive iNHB favours CardioMEMS
## negative iNHB favours usual care

WTP <- 50000
df_iNHB <- df_iNHB %>%
  mutate(
    iNHB = (e_CardioMEMS - e_UsualCare)-
      ((c_CardioMEMS - c_UsualCare)/WTP)
  )
## Check the resulting dataset and distribution of iNHB
print(nrow(df_iNHB))
print(head(df_iNHB))
print(names(df_iNHB))

mean(df_iNHB$iNHB)
summary(df_iNHB$iNHB)

## Fit a linear regression metamodel with iNHB as the outcome and the uncertain model parameters sampled in the PSA as predictors
  fit_iNHB <- lm(
    iNHB ~ 
      p_cImp_comp + 
      p_HRatio_Hosp + 
      p_pDie_StableHF + 
      p_pHosp_StableHF_cond + 
      p_pImp_comp + 
      p_uHosp +
      p_uStableHF +
      p_uStableHFPostDischarge30D + 
      p_uStableHFPostDischarge30D_comp,
    data = df_iNHB
  )
  
  ## Examine the regression coeffecients and model fit
  ## The R-squared indicates how well the metamodel reproduces variation in iNHB observed across the PSA simulations.
  
  summary(fit_iNHB)
  
 parameter_ranges <- df_iNHB %>%
  dplyr::select(
    p_cImp_comp, p_HRatio_Hosp, p_pDie_StableHF,
    p_pHosp_StableHF_cond, p_pImp_comp, p_uHosp,
    p_uStableHF, p_uStableHFPostDischarge30D,
    p_uStableHFPostDischarge30D_comp
  )%>%
   summarize(across(everything(),
                    list(low = ~quantile(.x, 0.025),
                         high = ~quantile(.x, 0.975))))
 parameter_ranges
  
 base_vals <- df_iNHB %>%
   dplyr::select(all_of(names(coef(fit_iNHB))[-1])) %>%
   dplyr::summarize(across(everything(), mean))
 
 tornado <- purrr::map_dfr(names(base_vals), \(p) {
   lo <- hi <- base_vals
   lo[[p]] <- parameter_ranges[[paste0(p, "_low")]]
   hi[[p]] <- parameter_ranges[[paste0(p, "_high")]]
   
   tibble(
     parameter = p,
     low = predict(fit_iNHB, lo),
     high = predict(fit_iNHB, hi)
   )
 }) %>%
   mutate(impact = abs(high-low),
          parameter = reorder(parameter, impact))
 
   ggplot(tornado, aes(y = parameter))+
     geom_segment(aes(x = low, xend = high, yend = parameter), linewidth = 6) + 
     geom_vline(xintercept = mean(df_iNHB$iNHB), linetype = "dashed") +
     scale_y_discrete(labels = c(
    "p_HRatio_Hosp" = "HRatio_Hosp",
    "p_pHosp_StableHF_cond" = "pHosp_StableHF_cond",
    "p_uStableHFPostDischarge30D_comp" = "uStableHFPostDischarge30D_comp",
    "p_uStableHFPostDischarge30D" = "uStableHFPostDischarge30D",
    "p_pDie_StableHF" = "pDie_StableHF",
    "p_uHosp" = "uHosp",
    "p_uStableHF" = "uStableHF",
    "p_Imp_comp" = "pImp_comp",
    "p_cImp_comp" = "cImp_comp"
     ))+
     labs(
       x = "iNHB",
       y = "Parameters",
       title = "Tornado Diagram"
     ) +
     theme_minimal()
     
   tornado
