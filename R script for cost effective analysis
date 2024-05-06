## setwd("/Users/priyanka/Desktop/Health Economics Project")
## install.packages("dampack")    # for CEA and calculate ICERs
library(dampack)


## Treatment 1 -> Face to Face CBT 
## Treatment 2 -> Online CBT

## General Setup
cycle_length <- 1 # cycle length 
n_cycles <- 60 # 5 year time horizon, total number of cycles
v_names_states <- c("D1", "C", "M", "R","D2") # the 5 states of the model: 
                                                # Depression (D1), CBT (C), CBT + Antidepressants (M), Remission (R), Death (D2)
n_states <- length(v_names_states)
d_e <- 0.03 # annual discount rate for QALYs of 3%
d_c <- 0.03 # annual discount rate for costs of 3%
v_names_str <- c("F2F CBT", "Online CBT")

#################### Defining transition probabilities for both treatment interventions ##############################

## Transition probabilities (monthly)

p_D1D2 <- 0.010438 #1. Depression to Death
p_D1C_trt1 <- 0.31 #2. Depression to Face to Face CBT 
p_D1C_trt2 <- 0.35 #3. Depression to Online CBT 

p_CM <- 0.33640 #4. CBT to Medication 
p_CD1_trt1 <- 0.0145  #5. Face to Face CBT to Depression 
p_CD1_trt2 <- 0.0201  #6. Face to Face CBT to Depression 
p_CR_trt1 <- 0.5174  #7. Face to Face CBT to Remission
p_CR_trt2 <- 0.6090 #8. Face to Face CBT to Remission
p_CM <- 0.33640 #9. CBT (Face to Face or Online) to Medication
p_CD2 <- 0.010438 #10. CBT (Face to Face or Online) to Death 

p_MD2 <- 0.010438 #11. Medication to Death
p_MD1 <- 0.0858 #12. Medication to Depression
p_MR <- 0.450 #13. Medication to Remission

p_RD1 <- 0.0064  #14. Remission to Depression
p_RD2 <- 0.010438 #15. Remission to Death
p_D2D2 <- 1 #16. Death to Death

## State rewards 
## Costs 
c_D1 <- 1142 # monthly cost of being depressed 
c_Ctrt1 <- 1407 # monthly cost of getting F2F CBT 
c_Mtrt1 <- 1483 # monthly cost of getting Medication and F2F CBT 
c_R <- 204 # monthly cost of being Healthy (in Remission State)
c_D2 <- 0 # monthly cost of being Dead
c_Ctrt2 <- 931 # monthly cost of getting Online CBT 
c_Mtrt2 <- 1007 # monthly cost of getting Medication and Online CBT 

## Utilities 
u_R <- 0.85 # monthly utility of being Healthy (in Remission State)
u_C <- 0.72 # monthly utility of getting CBT 
u_D1 <- 0.58 # monthly utility of being depressed 
u_D2 <- 0 # monthly utility of being dead 
u_M <- 0.72 # monthly utility of getting CBT + Medication 

v_m_init <- c(D1 = 1, C = 0, M = 0, R = 0, D2 = 0) # initial state vector 

## Initialize cohort trace for Treatment 1
m_M <- matrix(NA, nrow = (n_cycles + 1), ncol = n_states, dimnames = list(0:n_cycles, v_names_states))

# Storing the initial state vector in the first row of the cohort trace
m_M[1,] <- v_m_init

## Initialize cohort trace for Treatment 1 and 2
m_M_trt1 <- m_M # Treatment 1
m_M_trt2 <- m_M # Treatment 2

## Initialize transition probability matrix for Treatment 1 
m_P <- matrix(0, nrow = n_states, ncol = n_states, dimnames = list(v_names_states, v_names_states))

## Fill the matrix for Treatment 1
# From D1
m_P["D1","D1"] <- (1-p_D1D2) * (1-p_D1C_trt1)
m_P["D1","C"] <- (1-p_D1D2) * p_D1C_trt1
m_P["D1","D2"] <- p_D1D2

# From C
m_P["C","D1"] <- (1-p_CD2) * p_CD1_trt1
m_P["C","M"] <- (1-p_CD2) * p_CM
m_P["C","C"] <- (1-p_CD2) * (1 - (p_CD1_trt1 + p_CM + p_CR_trt1))
m_P["C","R"] <- (1-p_CD2) * p_CR_trt1
m_P["C","D2"] <- p_CD2

# from M
m_P["M","D1"] <- (1-p_MD2) * p_MD1
m_P["M","M"] <- (1-p_MD2) * (1 - (p_MD1 + p_MR))
m_P["M","R"] <- (1-p_MD2) * p_MR
m_P["M","D2"] <- p_MD2

# from R 
m_P["R","D1"] <- (1-p_RD2) * p_RD1
m_P["R","R"] <- (1-p_RD2) * (1-p_RD1)
m_P["R","D2"] <- p_RD2

# from D
m_P["D2","D2"] <- 1

## Initialize transition probability matrix for treatment 1 
m_P_trt1 <- m_P

## Initialize transition probability matrix for treatment 2 
m_P_trt2 <- m_P
### Update the transition probability matrix for treatment 2
m_P_trt2["D1","D1"] <- (1-p_D1D2) * (1-p_D1C_trt2)
m_P_trt2["D1","C"] <- (1-p_D1D2) * p_D1C_trt2
m_P_trt2["C","D1"] <- (1-p_CD2) * p_CD1_trt2
m_P_trt2["C","C"] <- (1-p_CD2) * (1 - (p_CD1_trt2 + p_CM + p_CR_trt2))
m_P_trt2["C","R"] <- (1-p_CD2) * p_CR_trt2

### Check if transition probability matrices are valid 
## Check that transition probabilities are between [0,1]
m_P >= 0 & m_P <= 1
m_P_trt1 >= 0 & m_P_trt1 <= 1
m_P_trt2 >= 0 & m_P_trt2 <= 1

## Check that all rows sum to 1
rowSums(m_P) == 1
rowSums(m_P_trt1) == 1
rowSums(m_P_trt2) == 1

# Iterative solution of time independent cSTM
for (t in 1:n_cycles){
  # For treatment 1
  m_M_trt1[t+1,] <- m_M_trt1[t,] %*% m_P_trt1
  # For treatment 2
  m_M_trt2[t+1,] <- m_M_trt2[t,] %*% m_P_trt2
  
}


## Cost and Effectiveness Outcomes 
# Vector of state utilities under Treatment 1
v_u_trt1 <- c(D1 = u_D1, C = u_C, M = u_M, R = u_R, D2 = u_D2) * cycle_length
# Vector of state utilities under Treatment 2
v_u_trt2 <- c(D1 = u_D1, C = u_C, M = u_M, R = u_R, D2 = u_D2) * cycle_length


# Vector of state costs under Treatment 1
v_c_trt1 <- c(D1 = c_D1, C = c_Ctrt1, M = c_Mtrt1, R = c_R, D2 = c_D2) * cycle_length
# Vector of state costs under Treatment 2
v_c_trt2 <- c(D1 = c_D1, C = c_Ctrt2, M = c_Mtrt2, R = c_R, D2 = c_D2) * cycle_length

# Vector of QALYs under Treatment 1
v_qaly_trt1 <- m_M_trt1 %*% v_u_trt1
# Vector of QALYs under Treatment 2
v_qaly_trt2 <- m_M_trt2 %*% v_u_trt2

# Vector of costs under Treatment 1
v_cost_trt1 <- m_M_trt1 %*% v_c_trt1
# Vector of costs under Treatment 2
v_cost_trt2 <- m_M_trt2 %*% v_c_trt2

## Discounting Future Rewards 

# Discount weights for effects 
v_dwe <- 1/((1+(d_e * cycle_length))^(0:n_cycles))
# Discount weights for costs
v_dwc <- 1/((1+(d_c * cycle_length))^(0:n_cycles))

## Expected discounted QALYs under treatment 1
n_tot_qaly_trt1 <- t(v_qaly_trt1)%*%(v_dwe)
## Expected discounted costs under treatment 1
n_tot_cost_trt1 <- t(v_cost_trt1)%*%(v_dwc)

## Expected discounted QALYs under treatment 2
n_tot_qaly_trt2 <- t(v_qaly_trt2)%*%(v_dwe)
## Expected discounted costs under treatment 2
n_tot_cost_trt2 <- t(v_cost_trt2)%*%(v_dwc)

### Vector of costs 
v_cost_str <- c(n_tot_cost_trt1,n_tot_cost_trt2)

### Vector of effectiveness 
v_qaly_str <- c(n_tot_qaly_trt1,n_tot_qaly_trt2 )

### Calculate cost effectiveness ratios (ICERs)
df_cea <- dampack::calculate_icers(cost = v_cost_str, effect = v_qaly_str, strategies = v_names_str)




