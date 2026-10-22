#clear my environment
rm(list=ls())

#clear my R console
cat("\014")

library(tidyverse)
library(twig)

#Definition for key acronyms
#Imp_comp - implant-related complication
#pHosp_StableHFPostDischarge30D_cond - probability of hospitalization given survival, in the stableHFpostdischarge30D state
#pHosp_StableHF_cond - probability of hospitalization given survival, in the stableHF state


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
          discount_rates=c(0.015,0.015))

# none is the combination of no + NA
# in the death state everyone is dead already, so they will be in the leftover section

# Define Parameters
params = list(uStableHFPostDischarge30D	= 0.790,
              uStableHFPostDischarge30D_comp = 0.743,
              uStableHF	= 0.780,
              uHosp= 0.620,
              uDead	= 0,
              
              
              cOutpatientMonthly = 100.22,
              cHospEvent = 10753.56,
              cCardioMEMSDevice = 18860.22,
              cCardioMEMSImplant	= 15561.62,
              cImp_comp = 8808.44,
              cCardioMEMSMonitoringMonthly	= 131.52,
              cDead = 0,
              
              
              pDie_StableHFPostDischarge30D	= 0.087,
              pDie_StableHF	= 0.0108,
              
              ##Raw hospitalization probabilities (monthly, unconditional)
              pHosp_StableHFPostDischarge30D = 0.121,
              pHosp_StableHF = 0.048,
              
              ## Probability of complication
              pImp_comp = 0.026,
              ## Hospitaliziation among survivors 
              pHosp_StableHFPostDischarge30D_cond = 0.115,
              pHosp_StableHF_cond = 0.0437,                                     
              
              HRatio_Hosp =	0.7)
              
              
      
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

# Utility adjustmebt if implant-related complications, not tretaing it as a state utility
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

options(max.print = 1000000)
# Run Model(determinsitic)
## # Run the decision-analytic model over 360 cycles using specified parameters
results_deterministic <- run_twig(twig_obj = mytwig, params = params, n_cycles = 360, verbose= TRUE)

## View mean expected costs and utilities for each treatment strategy
results_deterministic$mean_ev 
## Inspect the state transition probabilities used by the Markov model (for validation and debugging)
### Is the model behaving as it should, do all the transition probabilities add up to 1
results_deterministic$markov_trans_probs

getAnywhere(run_twig)