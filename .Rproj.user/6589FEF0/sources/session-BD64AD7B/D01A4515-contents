# model code
# looking at many types of structural uncertainity:
# 1. SIX, SIS, SIR, SIRS (4)
# 2. different orders of transmission/recovery/mortality (6 combos)
# 3. alternative disease costs (2)
# 4. alternative types of resistance (3)

library(dplyr)

# parameters & start conditions
n_RR = 5; n_WR = 0; n_WW = 20 # individuals by genotype
inf_status = c(rep("S", 2), rep("I", sum(n_RR, n_WR, n_WW)-2)) # generic to start
b_RR = 2; b_WR = 2; b_WW = 2 # transmission RATE
b_sd = 0.1 # enviro effect
m_SRR = 0.1; m_SWR = 0.1; m_SWW = 0.1 # sus mort rate
m_IRR = 1; m_IWR = 1; m_IWW = 1 # mortality RATE
m_Ssd = 0.1; m_Isd =0.1
r_RR = 1; r_WR = 0; r_WW = 0 # recovery RATE
r_sd = 0.1
l_RR = 2; l_WR = 4; l_WW = 4; l_sd = 0.5; K = 150; K_sd = 10 # reproduction things

# model input will be a data frame of individuals at start of generation
inds <- data.frame(ind_num = 1:sum(n_RR, n_WR, n_WW),
                   inf_stat = inf_status,
                   ind_geno = c(rep("RR", n_RR), rep("WR", n_WR), rep("WW", n_WW)),
                   b_pheno = c(rnorm(n_RR, mean=b_RR, sd=b_sd), 
                               rnorm(n_WR, mean=b_WR, sd=b_sd), 
                               rnorm(n_WW, mean=b_WW, sd=b_sd)),
                   mS_pheno = c(rnorm(n_RR, mean=m_SRR, sd=m_Ssd), 
                               rnorm(n_WR, mean=m_SWR, sd=m_Ssd), 
                               rnorm(n_WW, mean=m_SWW, sd=m_Ssd)),
                   mI_pheno = c(rnorm(n_RR, mean=m_IRR, sd=m_Isd), 
                                rnorm(n_WR, mean=m_IWR, sd=m_Isd), 
                                rnorm(n_WW, mean=m_IWW, sd=m_Isd)),
                   r_pheno = c(rnorm(n_RR, mean=r_RR, sd=r_sd), 
                               rnorm(n_WR, mean=r_WR, sd=r_sd), 
                               rnorm(n_WW, mean=r_WW, sd=r_sd))
                   ) # phenos still in rates

ngens = 40

# run_gens <- function(compartments, events, inds) {
  # note this takes in a compartment structure SIS, SIX, SIR 
  # second takes in an order: events order (e1, e2, e3) (B,M,R)
  # other uncertainity comes from parameters above: disease cost & alt resistances
  
  freq_R <- (length(which(inds$ind_geno == "WR")) + 2*length(which(inds$ind_geno == "RR")))/(2*dim(inds)[1])
  tot_pop <- dim(inds)[1] # every row individual 
  prev_I <- length(which(inds$inf_stat == "I"))/dim(inds)[1]
  
  for (i in 1:ngens) {
    # draw probabilities for all process each year --> need to truncate at 0?
    # phenotypes stay the same
    inds$transmission <- runif(dim(inds)[1]) 
    inds$mortalityS <- runif(dim(inds)[1])
    inds$mortalityI <- runif(dim(inds)[1])
    inds$recovery <- runif(dim(inds)[1])
    
    # do the three events per generation
    for (j in 1:length(events)) {
      if (events[j]=="B") {
        # get transmission phenotype
        inds$change_stat <- (1-exp(-inds$b_pheno*length(which(inds$inf_stat=="I")))) > inds$transmission
        tmpS <- inds %>% filter(inf_stat=="S" & change_stat==T)
        if (dim(tmpS)[1] > 0) tmpS$inf_stat <- "I" # goes to being an I
        # everyone else
        tmpI <- inds %>% filter(inf_stat!="S" | (inf_stat=="S" & change_stat==F))
        # recombine
        inds <- rbind(tmpS, tmpI)
        inds <- inds[, 1:11]
      }
      if (events[j]=="M") {
        inds$p_survS <- exp(-inds$mS_pheno)
        inds$p_survI <- exp(-inds$mI_pheno)
        # get mortality phenotype: infects
        tmpI <- inds %>% filter(inf_stat == "I") %>% filter(mI_pheno < mortalityI)
        tmpS <- inds %>% filter(inf_stat == "S") %>% filter(mS_pheno < mortalityS)
        inds <- rbind(tmpI, tmpS)
        inds <- inds[, 1:11]
      }
      if (events[j]=="R" & compartments !="SIX") {
        # get recovery phenotype
        inds$change_stat <- (1-exp(-inds$r_pheno)) > inds$recovery
        # get I --> not I list
        tmpI <- inds %>% filter(inf_stat=="I" & change_stat==T)
        if (dim(tmpI)[1] > 0) tmpI$inf_stat <- ifelse(compartments=="SIS", "S", "R")
        # everyone else
        tmpS <- inds %>% filter(inf_stat!="I" | (inf_stat=="I" & change_stat==F))
        # recombine
        inds <- rbind(tmpS, tmpI)
        inds <- inds[, 1:11]
      }
    }
    
    if (dim(inds)[1] == 0) break
    
    # survivors reproduce
    parent_lams <- paste0("l_",inds$ind_geno,sep="")
    rep_rate <- NULL
    for (m in 1:length(parent_lams)) {
         lam_val=ifelse(parent_lams[m]=="l_RR", l_RR, ifelse(parent_lams[m]=="l_WR", l_WR, l_WW))
         rep_rate <- c(rep_rate, floor(rnorm(1, lam_val, l_sd)))
    }
    rep_rate[which(rep_rate < 0)] = 0
    gamts <- NULL
    for (m in 1:length(parent_lams)) {
      # get gametes from parent
      gamts_temp1 <- if (inds$ind_geno[m]=="RR") {rep("R", rep_rate[m])}
      gamts_temp2 <- if (inds$ind_geno[m]=="WW") {rep("W", rep_rate[m])}
      # WRs will need to have coin flip if odd rep_rate
      gamts_temp3 <- if (inds$ind_geno[m]=="WR" & rep_rate[m]%%2==0) {c(rep("W", rep_rate[m]/2), rep("R", rep_rate[m]/2))}
      gamts_temp4 <- if (inds$ind_geno[m]=="WR" & rep_rate[m]%%2==1) {c(rep("W", (rep_rate[m]-1)/2), rep("R", (rep_rate[m]-1)/2), sample(c("W","R"), 1))}
      gamts <- c(gamts, gamts_temp1, gamts_temp2, gamts_temp3, gamts_temp4)
    }
    K_stoch = floor(rnorm(1, mean=K, sd=K_sd))
    off_dat <- NULL
    # if there are 2+ gametes, create offspring
    if (length(gamts) > 1 & K_stoch-length(inds$ind_num) > 0) {
      # draw them, limiting by stochastic K
      off_gamts <- sample(gamts, 
                          min(length(gamts), (K_stoch-length(inds$ind_num))*2), 
                          replace=FALSE)
      # only take even number: if odd length, drop the last one
      if (length(off_gamts)%%2==1) off_gamts <- off_gamts[1:length(off_gamts)-1]
      # grab pairs
      off_genos <- NULL
      for (n in 1:(length(off_gamts)/2)) {
        geno_temp1 <- if (off_gamts[2*n] == "W" & off_gamts[2*n-1] == "W") "WW"
        geno_temp2 <- if (off_gamts[2*n] == "R" & off_gamts[2*n-1] == "R") "RR"
        # heterozygous case
        geno_temp3 <- if (off_gamts[2*n] == "R" & off_gamts[2*n-1] == "W") "WR"
        geno_temp4 <- if (off_gamts[2*n] == "W" & off_gamts[2*n-1] == "R") "WR"
        off_genos <- c(off_genos, geno_temp1, geno_temp2, geno_temp3, geno_temp4)
      }
      off_dat <- data.frame(ind_num = 1:length(off_genos),
                            inf_stat = rep("S"),
                            ind_geno = off_genos)
      off_dat$b_pheno <- ifelse(off_dat$ind_geno == "WW", rnorm(1, mean=b_WW, sd=b_sd),
                                ifelse(off_dat$ind_geno == "RR", rnorm(1, mean=b_RR, sd=b_sd),
                                                                       rnorm(1, mean=b_WR, sd=b_sd)))
      off_dat$mS_pheno <- ifelse(off_dat$ind_geno == "WW", rnorm(1, mean=m_SWW, sd=m_Ssd),
                                ifelse(off_dat$ind_geno == "RR", rnorm(1, mean=m_SRR, sd=m_Ssd),
                                                                       rnorm(1, mean=m_SWR, sd=m_Ssd)))
      off_dat$mI_pheno <- ifelse(off_dat$ind_geno == "WW", rnorm(1, mean=m_IWW, sd=m_Isd),
                                 ifelse(off_dat$ind_geno == "RR", rnorm(1, mean=m_IRR, sd=m_Isd),
                                        rnorm(1, mean=m_IWR, sd=m_Isd)))
      off_dat$r_pheno <- ifelse(off_dat$ind_geno == "WW", rnorm(1, mean=r_WW, sd=r_sd),
                                ifelse(off_dat$ind_geno == "RR", rnorm(1, mean=r_RR, sd=r_sd),
                                                                       rnorm(1, mean=r_WR, sd=r_sd)))
                            
    }

    # combine parents and offspring
    inds <- rbind(inds[,1:7], off_dat)
    # reindex 
    inds$ind_num <- 1:dim(inds)[1]

    freq_R <- c(freq_R, (length(which(inds$ind_geno == "WR")) + 2*length(which(inds$ind_geno == "RR")))/(2*dim(inds)[1]))
    tot_pop <- c(tot_pop, dim(inds)[1]) # every row individual 
    prev_I <- c(prev_I, length(which(inds$inf_stat == "I"))/dim(inds)[1])
    
  }
  
  plot(freq_R); plot(tot_pop); plot(prev_I)
  
# }






