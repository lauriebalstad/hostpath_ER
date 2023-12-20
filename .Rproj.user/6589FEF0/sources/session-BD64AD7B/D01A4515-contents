# model code
# looking at many types of structural uncertainity:
# 1. SIX, SIS, SIR (3)
# 2. different orders of transmission B/recovery G/mortality M (6 combos)
# 3. alternative types of resistance (3)

library(dplyr)

run_gens <- function(comp_str, event_order, # structural
                     n_RR, n_WR, n_WW, init_I, # set ups
                     b_RR, b_WR, b_WW, b_sd, # transmission
                     m_SRR, m_SWR, m_SWW, m_IRR, m_IWR, m_IWW, m_Ssd, m_Isd, # mortality
                     r_RR, r_WR, r_WW, r_sd, # recovery
                     l_RR, l_WR, l_WW, l_sd, # reproduction
                     K, K_sd, ngens, bnum) { # pop size and gen time info
  # note this takes in a compartment structure SIS, SIX, SIR 
  # second takes in an order: events order (e1, e2, e3) (B,M,R)
  # and then all the parameters
  
  # get compartmetns
  if(comp_str == 1) compartments <- c("SIX")
  if(comp_str == 2) compartments <- c("SIS")
  if(comp_str == 3) compartments <- c("SIR")
  
  # get event order list
  if(event_order == 1) events <- c("B","M","G")
  if(event_order == 2) events <- c("B","G","M")
  if(event_order == 3) events <- c("M","G","B")
  if(event_order == 4) events <- c("M","B","G")
  if(event_order == 5) events <- c("G","B","M")
  if(event_order == 6) events <- c("G","M","B")
  
  # storage vectors for across simulations
  extinct <- NULL # did the population go extinct? T/F
  pop_size_tf_all <- NULL # final population size
  pop_size_tf_not_extinct <- NULL # conditioned on not being extinct
  at_K <- NULL # is the final population at the carrying capactity? T/F
  time_first_at_K <- NULL # and how fast did demographic rescue happen
  infect_class_I_tf_all <- NULL # final infection prevelence
  infect_class_I_not_extinct <- NULL # conditioned on not being extinct
  lost_I <- NULL # did I class go away? T/F
  time_lost_I <- NULL # and at what time? 
  freq_R_allele_tf <- NULL # final frequency of the R allele
  max_freq_R_allele <- NULL # was there a high point? trying to capture if back selection occured
  time_max_freq_R_allele <- NULL # and what time was the high point?
  
  for (b in 1:bnum) {
    # get inf_stat
    inf_status <- sample(c("S", "I"), 
                         size = sum(n_RR, n_WR, n_WW), 
                         replace = T, 
                         prob = c(1-init_I, init_I))
    # construct initial inds
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
    inds$mI_pheno = ifelse(inds$mI_pheno<0, inds$mS_pheno, inds$mS_pheno+inds$mI_pheno)
    # check for any below 0, change to positive
    if (any(inds<0)) {
      index_neg <- which(inds<0)
      for (i in 1:length(index_neg)) {
        col_num <- ceiling(index_neg[i]/dim(inds)[1])
        row_num <- index_neg[i]%%dim(inds)[1]
        inds[row_num,col_num] <- 0
      }
    }
    
    # within simulation info
    I_num <- sum(inds$inf_stat == "I")
    R_freq <- (2*(sum(inds$ind_geno == "RR")) + sum(inds$ind_geno == "WR"))/(2*dim(inds)[1])
    pop_size <- dim(inds)[1]
    carry_capacity <- floor(rnorm(1, mean=K, sd=K_sd)) # start with random K --> note K != pop for first init
    extinct_dummy <- FALSE
    
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
          tmpS <- inds %>% filter(inf_stat == "S" | inf_stat == "R") %>% filter(mS_pheno < mortalityS)
          inds <- rbind(tmpI, tmpS)
          inds <- inds[, 1:11]
        }
        if (events[j]=="G" & compartments !="SIX") {
          # get recovery phenotype
          inds$change_stat <- (1-exp(-inds$r_pheno)) > inds$recovery
          # get I --> not I list
          tmpI <- inds %>% filter(inf_stat=="I" & change_stat==T)
          if (dim(tmpI)[1] > 0) tmpI$inf_stat <- ifelse(compartments=="SIR", "R", "S")
          # everyone else
          tmpS <- inds %>% filter(inf_stat!="I" | (inf_stat=="I" & change_stat==F))
          # recombine
          inds <- rbind(tmpS, tmpI)
          inds <- inds[, 1:11]
        }
      }
      
      if (dim(inds)[1] == 0) {
        extinct_dummy <- TRUE # pop is extinct
        I_num <- c(I_num, 0) # no Is
        R_freq <- c(R_freq, 0) # no Rs
        pop_size <- c(pop_size, 0) # no one
        carry_capacity <- c(carry_capacity, floor(rnorm(1, mean=K, sd=K_sd))) # give extra K for matching
        break # stop generation loop
        }
      
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
                              ind_geno = off_genos, 
                              b_pheno = NA,
                              mS_pheno = NA,
                              mI_pheno=NA,
                              r_pheno=NA)
        for (p in 1:dim(off_dat)[1]) {
        off_dat$b_pheno[p] <- ifelse(off_dat$ind_geno[p] == "WW", rnorm(1, mean=b_WW, sd=b_sd),
                                  ifelse(off_dat$ind_geno[p] == "RR", rnorm(1, mean=b_RR, sd=b_sd),
                                                                         rnorm(1, mean=b_WR, sd=b_sd)))
        off_dat$mS_pheno[p] <- ifelse(off_dat$ind_geno[p] == "WW", rnorm(1, mean=m_SWW, sd=m_Ssd),
                                  ifelse(off_dat$ind_geno[p] == "RR", rnorm(1, mean=m_SRR, sd=m_Ssd),
                                                                         rnorm(1, mean=m_SWR, sd=m_Ssd)))
        off_dat$mI_pheno[p] <- ifelse(off_dat$ind_geno[p] == "WW", rnorm(1, mean=m_IWW, sd=m_Isd),
                                   ifelse(off_dat$ind_geno[p] == "RR", rnorm(1, mean=m_IRR, sd=m_Isd),
                                          rnorm(1, mean=m_IWR, sd=m_Isd)))
        off_dat$r_pheno[p] <- ifelse(off_dat$ind_geno[p] == "WW", rnorm(1, mean=r_WW, sd=r_sd),
                                  ifelse(off_dat$ind_geno[p] == "RR", rnorm(1, mean=r_RR, sd=r_sd),
                                                                         rnorm(1, mean=r_WR, sd=r_sd)))
        }
        off_dat$mI_pheno = ifelse(off_dat$mI_pheno<0, off_dat$mS_pheno, off_dat$mS_pheno+off_dat$mI_pheno)
        # check for any below 0, change to positive
        if (any(off_dat<0)) {
          index_neg <- which(off_dat<0)
          for (i in 1:length(index_neg)) {
            col_num <- ceiling(index_neg[i]/dim(off_dat)[1])
            row_num <- index_neg[i]%%dim(off_dat)[1]
            off_dat[row_num,col_num] <- 0
          }
        }
      }
      
      # combine parents and offspring
      inds <- rbind(inds[,1:7], off_dat)
      # reindex 
      inds$ind_num <- 1:dim(inds)[1]
      
      I_num <- c(I_num, sum(inds$inf_stat == "I"))
      R_freq <- c(R_freq, (2*(sum(inds$ind_geno == "RR")) + sum(inds$ind_geno == "WR"))/(2*dim(inds)[1]))
      pop_size <- c(pop_size, dim(inds)[1])
      carry_capacity <- c(carry_capacity, K_stoch)
      
    }
    
    extinct <- c(extinct, extinct_dummy) # did the population go extinct? T/F
    pop_size_tf_all <- c(pop_size_tf_all, last(pop_size)) # final population size
    pop_size_tf_not_extinct <- if (extinct_dummy==F) c(pop_size_tf_not_extinct, last(pop_size))
    at_K <- c(at_K, any(carry_capacity==pop_size)) # is the population ever at the carrying capactity? T/F
    time_first_at_K <- c(time_first_at_K, which(carry_capacity==pop_size)[1]) # and how fast did demographic rescue happen
    infect_class_I_tf_all <- c(infect_class_I_tf_all, last(I_num)) # final infection prevelence
    infect_class_I_not_extinct <- if (extinct_dummy==F) c(infect_class_I_not_extinct, last(I_num))
    lost_I <- c(lost_I, any(I_num == 0)) # did I class go away? T/F
    time_lost_I <- c(time_lost_I, which(I_num == 0)[1])  # and at what time? 
    freq_R_allele_tf <- c(freq_R_allele_tf, last(R_freq)) # final frequency of the R allele
    max_freq_R_allele <- c(max_freq_R_allele, max(R_freq)) # was there a high point? trying to capture if back selection occured
    time_max_freq_R_allele <- c(time_max_freq_R_allele, which(R_freq == max(R_freq))[1]) # and what time was the high point?
    
  }
  
  # return summary statistics
  p_extinct <- sum(extinct)/length(extinct)
  avg_pop_size_tf <- mean(pop_size_tf_all)
  sd_pop_size_tf <- sd(pop_size_tf_all)
  avg_pop_size_given_NE <- mean(pop_size_tf_not_extinct)
  sd_pop_size_given_NE <- sd(pop_size_tf_not_extinct)
  p_demo_recovery <- sum(at_K)/length(at_K)
  avg_time_K <- mean(time_first_at_K, na.rm = T)
  sd_time_K <- sd(time_first_at_K, na.rm = T)
  avg_infect_class_I_tf <- mean(infect_class_I_tf_all)
  sd_infect_class_I_tf <- sd(infect_class_I_tf_all)
  avg_infect_class_I_given_NE <- mean(infect_class_I_not_extinct)
  sd_infect_class_I_given_NE <- sd(infect_class_I_not_extinct)
  p_lost_I <- sum(lost_I)/length(lost_I)
  avg_time_lost_I <- mean(time_lost_I, na.rm = T)
  sd_time_lost_I <- sd(time_lost_I, na.rm = T)
  avg_freq_R_tf <- mean(freq_R_allele_tf)
  sd_freq_R_tf <- sd(freq_R_allele_tf)
  avg_max_R <- mean(max_freq_R_allele)
  avg_time_max_R <- mean(time_max_freq_R_allele) # all will have a maximum
  
  return(list(p_extinct, p_demo_recovery, avg_time_K, sd_time_K,
              avg_pop_size_tf, sd_pop_size_tf, avg_pop_size_given_NE, sd_pop_size_given_NE,
              p_lost_I, avg_time_lost_I, sd_time_lost_I,
              avg_infect_class_I_tf, sd_infect_class_I_tf, avg_infect_class_I_given_NE, sd_infect_class_I_given_NE,
              avg_freq_R_tf, sd_freq_R_tf, avg_max_R, avg_time_max_R)
         )
  
}






