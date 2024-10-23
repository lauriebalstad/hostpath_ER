library(lhs) # sensitivity

# make sensitivity data

# # recall order: 
# comp_str <- parm_vect[1]; event_order <- parm_vect[2]; trans_type <- parm_vect[3] # structural
# d_0 <- parm_vect[4]; r_0 <- parm_vect[5] # disease intro and allele intro
# b_RR <- parm_vect[6]; b_WR <- parm_vect[7]; b_WW <- parm_vect[8]; b_sd <- parm_vect[9] # transmission
# m_SRR <- parm_vect[10]; m_SWR <- parm_vect[11]; m_SWW <- parm_vect[12]; m_IRR <- parm_vect[13]; m_IWR <- parm_vect[14]; m_IWW <- parm_vect[15]; m_Ssd <- parm_vect[16]; m_Isd <- parm_vect[17] # mortality
# r_RR <- parm_vect[18]; r_WR <- parm_vect[19]; r_WW <- parm_vect[20]; r_sd <- parm_vect[21] # recovery
# l_RR <- parm_vect[22]; l_WR <- parm_vect[23]; l_WW  <- parm_vect[24]; l_sd <- parm_vect[25]; mut_rate <- parm_vect[26] # reproduction
# K <- parm_vect[27]; K_sd <- parm_vect[28] # carrying capacity
# n <- parm_vect[27] # number of individuals at start of simulation, all WW

N <- 2000 # gives ~160 parameter sets per case, SIS comparision
num_parms <- 24 # only need 27 unique parameters because some repeat by defintion

mat <- randomLHS(n = N, k = num_parms) 
# need to update the matrices so that the right sorts of things are in the matrix
parm_mat <- matrix(NA, 2000, 31)
# fill in parm_mat
parm_mat[,1] <- rep(2, 2000) # hold at SIS, vect[1] = 2
parm_mat[,2] <- rep(2, 2000) # hold event order BGM
parm_mat[,3] <- rep(1, 2000) # hold transmission type at density dependence
parm_mat[,4] <- ceiling(qunif(mat[,1], 7, 12)) # "burn-in" pre-disease
parm_mat[,5] <- ceiling(qunif(mat[,2], 7, 12)) # "burn-in" post-disease (no adaptive allele)
parm_mat[,8] <- qnorm(mat[,5], 1, 0.05) # b_WW: transform into rate with mean 1 and sd 0.1
parm_mat[,6] <- parm_mat[,8]*mat[,3] # b_RR as percent of b_WW value (recall mat[7, ] is vector between 0-1)
parm_mat[,7] <- qunif(mat[,4],min=parm_mat[,6],max=parm_mat[,8]) # set b_WR between with mat[8,] dominance
parm_mat[, 9] <- qnorm(mat[,6], 0.05, 0.01) # compress b_sd 
parm_mat[,10] <- qnorm(mat[,7], 0.1, 0.025) # m_SXX: transform into rate with mean 0.1 and sd 0.02
parm_mat[,11] <- parm_mat[,10] # m_SWR = m_SRR
parm_mat[,12] <- parm_mat[,10] # m_SWW = m_SRR
# m_IXX is additive
parm_mat[,15] <- qnorm(mat[,10], 1, 0.1) # m_IWW
parm_mat[,13] <- parm_mat[,15]*mat[,8] # m_IRR as percent of m_IWW value
parm_mat[,14] <- qunif(mat[,9],min=parm_mat[,13],max=parm_mat[,15]) # set m_IWR between
parm_mat[,16] <- qnorm(mat[,11], 0.05, 0.01) # compress m_Xsd 
parm_mat[,17] <- parm_mat[,16] # m_Isd = m_Ssd
parm_mat[,18] <- qnorm(mat[,12], 0.5, 0.1) # r_RR: transform into rate with mean 0.2 and sd 0.1
parm_mat[,20] <- parm_mat[,18]*mat[,12] # r_WW as percent of r_RR value 
parm_mat[,19] <- qunif(mat[,13],min=parm_mat[,20],max=parm_mat[,18]) # set r_WR between with mat[8,] dominance
parm_mat[,21] <- qnorm(mat[,15], 0.05, 0.01) # compress r_sd 
parm_mat[,24] <- qnorm(mat[,16], 2.1, 0.1) # WW: population should grow, on average
parm_mat[,22] <- parm_mat[,24]*mat[,16] # RR: let mat[,22] be the cost
parm_mat[,23] <- qunif(mat[,17],min=parm_mat[,22],max=parm_mat[,24])
parm_mat[,25] <- qnorm(mat[,19], 0.1, 0.02) # compress l_sd 
parm_mat[,26] <- qnorm(mat[,20], 0.0075, 0.002) # mutation rate really little! 
parm_mat[,27] <- qunif(mat[,21], min = 70, max = 700) # K -- give good range
parm_mat[,28] <- qnorm(mat[,22], 4, 1) # keep K sd pretty small?
parm_mat[,29] <- ceiling(qunif(mat[,23], 0, 3)) # 1-3 disease cycles
parm_mat[,30] <- ceiling(qunif(mat[,24], 120, 140)) # ngens between 70-100

# check for zeros & replace them
parm_mat[which(parm_mat<0)] = 0

# indexing
parm_mat[,31] <- 1:N # this will be parm_num

# then will need to repeat the data frame 500 times, bc need 500 sims per parameter combination

# saved as sensitivity_data/mat_var.RDS
saveRDS(parm_mat, "sensitivity/mat_var.rds")