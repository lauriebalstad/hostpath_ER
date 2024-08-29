# make sensitivity data

mat <- randomLHS(n = N, k = length(params))# need to update the matrices so that the right sorts of things are in the matrix
mat[,1] <- ceiling(qunif(mat[,1], 0, 3))
mat[,2] <- ceiling(qunif(mat[,2], 0, 6))
mat[,3] <- ceiling(qunif(mat[,3], 0, 10)) # has to be positive whole numbers
mat[,4] <- ceiling(qunif(mat[,4], 0, 10)) # has to be positive whole numbers
mat[,5] <- ceiling(qunif(mat[,5], 0, 10)) # has to be positive whole numbers
# let init_I be any number 0-1 --> okay
mat[,9] <- qnorm(mat[,9], 1, 0.2) # transform into rate with mean 1 and sd 0.1
mat[,7] <- mat[,9]*mat[,7] # b_RR as percent of b_WW value
mat[,8] <- qunif(mat[,8],min=mat[,7],max=mat[,9]) # set b_WR between
# leave sd
mat[,13] <- qnorm(mat[,13], 0.1, 0.02) # transform into rate with mean 0.1 and sd 0.02
mat[,11] <- mat[,13]*mat[,11] # b_RR as percent of m_WW value
mat[,12] <- qunif(mat[,12],min=mat[,11],max=mat[,13]) # set m_WR between
# m_IXX is additive
mat[,16] <- qnorm(mat[,16], 0.9, 0.05) # transform into rate with mean 0.9 and sd 0.05
mat[,14] <- mat[,16]*mat[,14] # b_RR as percent of m_WW value
mat[,15] <- qunif(mat[,15],min=mat[,14],max=mat[,16]) # set m_WR between
# leave both sd
mat[,21] <- qnorm(mat[,21], 1, 0.1) # transform into rate with mean 1 and sd 0.1
mat[,19] <- mat[,21]*mat[,19] # r_RR as percent of b_WW value
mat[,20] <- qunif(mat[,20],min=mat[,19],max=mat[,21]) # set r_WR between
# leave sd
mat[,25] <- qnorm(mat[,25], 4, 0.1) # highest fitness for WWs
mat[,23] <- mat[,25]*mat[,23] # l_RR as percent of l_WW value
mat[,24] <- qunif(mat[,24],min=mat[,23],max=mat[,25]) # set l_WR between
# leave sd
mat[,27] <- qnorm(mat[,27], 180, 20) # keep pretty flat
mat[,28] <- qnorm(mat[,28], 10, 1) # keep wide-ish
# check for zeros
mat[which(mat<0)] = 0

# saved as sensitivity_data/mat_var.RDS