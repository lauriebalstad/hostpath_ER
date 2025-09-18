library(lhs) # sensitivity

# make sensitivity data
set.seed(204)

# # recall order: 
# comp_str <- parm_vect[1]; event_order <- parm_vect[2]; trans_type <- parm_vect[3] # structural
# d_0 <- parm_vect[4]; r_0 <- parm_vect[5] # disease intro and allele intro
# b_RR <- parm_vect[6]; b_WR <- parm_vect[7]; b_WW <- parm_vect[8]; b_sd <- parm_vect[9] # transmission
# m_SRR <- parm_vect[10]; m_SWR <- parm_vect[11]; m_SWW <- parm_vect[12]; m_IRR <- parm_vect[13]; m_IWR <- parm_vect[14]; m_IWW <- parm_vect[15]; m_Ssd <- parm_vect[16]; m_Isd <- parm_vect[17] # mortality
# r_RR <- parm_vect[18]; r_WR <- parm_vect[19]; r_WW <- parm_vect[20]; r_sd <- parm_vect[21] # recovery
# l_RR <- parm_vect[22]; l_WR <- parm_vect[23]; l_WW  <- parm_vect[24]; l_sd <- parm_vect[25]; mut_rate <- parm_vect[26] # reproduction
# K <- parm_vect[27]; K_sd <- parm_vect[28] # carrying capacity
# n <- parm_vect[27] # number of individuals at start of simulation, all WW

# 0131 N <- 2000 # gives ~160 parameter sets per case, SIS comparision
N <- 3000 # gives ~160 parameter sets per case, SIS comparision
num_parms <- 18 # only need 19 unique parameters 

mat <- randomLHS(n = N, k = num_parms) 
# need to update the matrices so that the right sorts of things are in the matrix
parm_mat <- matrix(NA, N, 37) # have 37 elements to fill in
# fill in parm_mat
parm_mat[,1] <- ceiling(qunif(mat[,1], min = 0, max = 3)) # compartments
parm_mat[,2] <- ceiling(qunif(mat[,2], min = 0, max = 6)) # event order

parm_mat[,37] <- ceiling(qunif(mat[,3], min = 0, max = 7)) # pathway: define this first, as it will inform others, but not part of the matrix
# 1: mortality, 2: transmission, 3: clearance
# 4: M+T, 5: M+C, 6: T+C, 7: M+T+C

parm_mat[,5] <- qunif(mat[,7], 1, 3) # f_WW: transform into rate with spanning negative (no enviro transmisison) to strong enviro transmission
parm_mat[,3] <- parm_mat[,5]*mat[,4] # f_RR as percent of f_WW value (recall mat is vector between 0-1)
parm_mat[,4] <- qunif(mat[,5],min=parm_mat[,3],max=parm_mat[,5]) # set f_WR between with mat dominance
parm_mat[, 6] <- qunif(mat[,6], 0.01, 0.1) # compress f_sd 

parm_mat[,9] <- qunif(mat[,8], 1, 3) # b_WW: transform into rate with mean 1 and sd 0.1
parm_mat[,7] <- parm_mat[,9]*mat[,4] # b_RR as percent of b_WW value (recall mat[7, ] is vector between 0-1)
parm_mat[,8] <- qunif(mat[,5],min=parm_mat[,7],max=parm_mat[,9]) # set b_WR between with mat[8,] dominance
parm_mat[, 10] <- qunif(mat[,6], 0.01, 0.1) # compress b_sd 

parm_mat[,11] <- qunif(mat[,9], 0.01, 0.1) # m_SXX: transform into rate with mean 0.1 and sd 0.02
parm_mat[,12] <- parm_mat[,11] # m_SWR = m_SRR
parm_mat[,13] <- parm_mat[,11] # m_SWW = m_SRR
# m_IXX is additive
parm_mat[,16] <- qunif(mat[,10], 1, 3) # m_IWW
parm_mat[,14] <- parm_mat[,16]*mat[,4] # m_IRR as percent of m_IWW value
parm_mat[,15] <- qunif(mat[,5],min=parm_mat[,14],max=parm_mat[,16]) # set m_IWR between
parm_mat[,17] <- qunif(mat[,6], 0.01, 0.1) # compress m_Xsd 
parm_mat[,18] <- parm_mat[,17] # m_Isd = m_Ssd

parm_mat[,19] <- qunif(mat[,11], 1, 3) # r_RR: now divding by the less than one value, since RR benefits increase gamma
parm_mat[,21] <- parm_mat[,19]*mat[, 4] # r_WW: avoids negative values
parm_mat[,20] <- qunif(mat[,5],min=parm_mat[,21],max=parm_mat[,19]) # set r_WR between with others with dominance
parm_mat[,22] <- qunif(mat[,6], 0.01, 0.1) # compress r_sd 

parm_mat[,25] <- qunif(mat[,12], 2, 2.5) # WW: population should grow, on average
parm_mat[,23] <- parm_mat[,25]*mat[,4] # RR: let mat[,21] be the cost
parm_mat[,24] <- qunif(mat[,5],min=parm_mat[,23],max=parm_mat[,25])
parm_mat[,26] <- qunif(mat[,6], 0.01, 0.1) # compress l_sd 
# parm_mat[,27] <- qnorm(mat[,24], 0.008, 0.001) # mutation rate really little! 
# parm_mat[,27] <- qunif(mat[,24], 0.0001, 0.05) # mutation rate really little! 
parm_mat[,27] <- qunif(mat[,13], 0.001, 0.01) # mutation rate really little! 

parm_mat[,28] <- qunif(mat[,14], 70, 700) # K -- give good range
parm_mat[,29] <- qunif(mat[,15], 3, 7) # keep K sd pretty small?
parm_mat[,30] <- rep(90) # d_0 as whole number
parm_mat[,31] <- ceiling(qunif(mat[,17], 1, 5)) # 2-5 disease cycles
parm_mat[,32] <- rep(0.1) # init r small 
parm_mat[,33] <- qunif(mat[,16], 0.05, 0.15) # proportion of I individuals (set to 0.1 in rep'd sims, if at 0 means 1 individual starts as infected)
parm_mat[,34] <- rep(120) # ngens between 120-140

# deal with pathways
# 1: mortality, 2: transmission, 3: clearance
# 4: M+T, 5: M+C, 6: T+C, 7: M+T+C
for (i in 1:dim(parm_mat)[1]) {
  # anything w mortality
  if (parm_mat[i,37] %in% c(2, 3, 6)) parm_mat[i,14:15] <- parm_mat[i,16] # no mortality: m_IWW == m_SWW
  if (parm_mat[i,37] %in% c(1, 3, 5)) parm_mat[i,4:5] <- parm_mat[i,3] # no transmission: density
  if (parm_mat[i,37] %in% c(1, 3, 5)) parm_mat[i,7:8] <- parm_mat[i,9] # no transmission: freq
  if (parm_mat[i,37] %in% c(1, 2, 4)) parm_mat[i,19:20] <- parm_mat[i,21] # no clearance
}

# density dependence only
parm_mat[,35] <- ceiling(qunif(mat[,18], min = 0, max = 3)) # 1 = density only, 2 = freq only, 3 = both
for (i in 1:dim(parm_mat)[1]) {
  if (parm_mat[i,35] == 1) parm_mat[i,3:5] <- -1 # density only
  if (parm_mat[i,35] == 2) parm_mat[i,7:9] <- -1 # freq only
}

# check for zeros & replace them
parm_mat[which(parm_mat[,c(1:2, 6, 10:29)]<0)] = 0 # allow environment to go negative, but nothing else

# indexing
parm_mat[,36] <- 1:N # this will be parm_num
bd_mat <- as.data.frame(mat[, c(4,5, 6)]); bd_mat$parm_number <- 1:N; colnames(bd_mat) <- c("benefit", "dominance", "number")

# then will need to repeat the data frame 500 times, bc need 500 sims per parameter combination

# saved as sensitivity_data/mat_var.RDS
# 0131 saveRDS(parm_mat, "dat/mat_var_0131.Rdata")
saveRDS(parm_mat, "dat/mat_var_0909.Rdata")
saveRDS(bd_mat, "dat/mat_bd_0909.Rdata")
