# model code
# looking at many types of structural uncertainity:
# 1. SIX, SIS, SIR (3)
# 2. different orders of transmission B/recovery G/mortality M (6 combos)
# 3. alternative types of resistance (3)

library(dplyr)
library(parallel)
library(foreach)
library(doParallel)

# so the plan is that the bnum reps is divided across cores

run_gens <- function(parm_vect, ngens, bnum) { # pop size and gen time info
  
  # note this takes in a compartment structure SIS, SIX, SIR 
  # second takes in an order: events order (e1, e2, e3) (B,M,R)
  # and then all the parameters
  comp_str <- parm_vect[1]; event_order <- parm_vect[2]; trans_type <- parm_vect[3] # structural
  d_0 <- floor(parm_vect[4]/parm_vect[29]); r_0 <- parm_vect[5] # disease intro and init allele freq.
  b_RR <- parm_vect[6]; b_WR <- parm_vect[7]; b_WW <- parm_vect[8]; b_sd <- parm_vect[9] # transmission
  m_SRR <- parm_vect[10]; m_SWR <- parm_vect[11]; m_SWW <- parm_vect[12]; m_IRR <- parm_vect[13]; m_IWR <- parm_vect[14]; m_IWW <- parm_vect[15]; m_Ssd <- parm_vect[16]; m_Isd <- parm_vect[17] # mortality
  r_RR <- parm_vect[18]; r_WR <- parm_vect[19]; r_WW <- parm_vect[20]; r_sd <- parm_vect[21] # recovery
  l_RR <- parm_vect[22]; l_WR <- parm_vect[23]; l_WW  <- parm_vect[24]; l_sd <- parm_vect[25]; mut_rate <- parm_vect[26] # reproduction
  K <- parm_vect[27]; K_sd <- parm_vect[28] # carrying capacity
  N <- parm_vect[27] # number of individuals at start of simulation, all WW -- start at K
  disease_cycles <- parm_vect[29] # number of times to go through disease between reproduction cycles
  
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
  
  # get transmission type 
  if(trans_type == 1) trans_form <- c("density")
  if(trans_type == 2) trans_form <- c("freq")
  
  all_reps <- foreach(b = 1:bnum, .packages="dplyr", .combine = 'rbind') %dopar% {
    
    # storage vectors for across simulations: will turn into a list and then rbind lists together
    extinct <- NULL # did the population go extinct? T/F
    S_size <- NULL # vector of S individuals
    I_size <- NULL # vector of I individuals
    R_size <- NULL # vector of R individuals
    K_size <- NULL # vector of K, to compare if K has been reached
    r_allele <- NULL # vector of R allele freq through time
    
    # make allele pool for genotypes
    init_allele <- sample(c("W", "R"), size = 2*N, prob = c(1-r_0, r_0), replace = T)
    init_genos <- NULL
    for (n in 1:(length(init_allele)/2)) {
      geno_temp1 <- if (init_allele[2*n] == "W" & init_allele[2*n-1] == "W") "WW"
      geno_temp2 <- if (init_allele[2*n] == "R" & init_allele[2*n-1] == "R") "RR"
      # heterozygous case
      geno_temp3 <- if (init_allele[2*n] == "R" & init_allele[2*n-1] == "W") "WR"
      geno_temp4 <- if (init_allele[2*n] == "W" & init_allele[2*n-1] == "R") "WR"
      init_genos <- c(init_genos, geno_temp1, geno_temp2, geno_temp3, geno_temp4)
    }
    
    # construct initial inds: all clean
    inds <- data.frame(ind_num = 1:N,
                       inf_stat = rep("S", N), # all healthy
                       ind_geno = init_genos, # use init genos
                       b_pheno = c(rnorm(N, mean=b_WW, sd=b_sd)), # transmission placeholder
                       mS_pheno = c(rnorm(N, mean=m_SWW, sd=m_Ssd)), # S_WW type mortality
                       mI_pheno = c(rnorm(N, mean=m_IWW, sd=m_Isd)), # mortality placeholder
                       r_pheno = c(rnorm(N, mean=r_WW, sd=r_sd)) # revovery placeholder 
    ) # phenos still in rates
    
    # check that northing is negative --> change to zero
    if (any(inds<0)) {
      index_neg <- which(inds<0)
      for (i in 1:length(index_neg)) {
        col_num <- ceiling(index_neg[i]/dim(inds)[1])
        row_num <- index_neg[i]%%dim(inds)[1]
        if (row_num == 0) {row_num <- dim(inds)[1]}
        inds[row_num,col_num] <- 0
      }
    }    
    # note mI_pheno is additional to mS_pheno --> address here
    inds$mI_pheno = inds$mS_pheno+inds$mI_pheno # m_I is a bonus mortality, should only increase mortality
    
    extinct_dummy <- FALSE
    
    # pre-disease intro
    for (p in 1:d_0) {
      # draw K first
      K_stoch = floor(rnorm(1, mean=K, sd=K_sd))
      # save things -- pre reproduction
      S_size <- c(S_size, dim(inds)[1]) # only Ss
      I_size <- c(I_size, 0) # no Is
      R_size <- c(R_size, 0) # no Rs
      K_size <- c(K_size, K_stoch) # carrying capacity
      r_allele <- c(r_allele, 0) # no r allele
      
      # no disease
      # so just mortality
      # probability of survival
      inds$p_survS <- exp(-inds$mS_pheno)
      # coin flips
      inds$mortalityS <- vapply(inds$p_survS, function(x) rbinom(1, 1, x), as.integer(1L))
      # get mortality phenotype: infects
      inds <- inds %>% filter(inf_stat == "S" | inf_stat == "R") %>% filter(mortalityS==1)
      inds <- inds[, 1:7] # remove extra columns
      
      # is everyone dead??
      if (dim(inds)[1] == 0) {
        extinct_dummy <- TRUE # pop is extinct
        S_size <- c(S_size, 0) # no Ss
        I_size <- c(I_size, 0) # no Is
        R_size <- c(R_size, 0) # no Rs
        K_size <- c(K_size, K_stoch) # carrying capacity
        r_allele <- c(r_allele, 0) # no R allele currently
        break # stop generation loop
      }
      
      # and then reproduce
      # how many gametes 
      if (dim(inds)[1] == 0) {
        # print(c(i, "everyone dead after disease dynamics"))
        extinct_dummy <- TRUE # pop is extinct
        S_size <- c(S_size, 0) # no Ss
        I_size <- c(I_size, 0) # no Is
        R_size <- c(R_size, 0) # no Rs
        K_size <- c(K_size, K_stoch) # carrying capacity
        r_allele <- c(r_allele, 0) # no R allele currently
        break # stop generation loop
      }
      
      # survivors reproduce if there is someone around
      if (dim(inds)[1] > 0) {
        # print(c(i, "dim(inds) > 0"))
        parent_lams <- paste0("l_",inds$ind_geno,sep="")
        rep_rate <- NULL
        for (P in 1:length(parent_lams)) {
          lam_val=ifelse(parent_lams[P]=="l_RR", l_RR, ifelse(parent_lams[P]=="l_WR", l_WR, l_WW))
          rep_rate <- c(rep_rate, floor(rnorm(1, lam_val, l_sd)))
        }
        rep_rate[which(rep_rate < 0)] = 0
        gamts <- NULL
        for (m in 1:length(parent_lams)) {
          # get gametes from parent
          gamts_temp1 <- if (inds$ind_geno[m]=="RR") {rep("R", 2*rep_rate[m])}
          gamts_temp2 <- if (inds$ind_geno[m]=="WW") {rep("W", 2*rep_rate[m])}
          # WRs will need to have coin flip if odd rep_rate
          gamts_temp3 <- if (inds$ind_geno[m]=="WR" & rep_rate[m]%%2==0) {c(rep("W", rep_rate[m]/2), rep("R", rep_rate[m]/2))}
          gamts_temp4 <- if (inds$ind_geno[m]=="WR" & rep_rate[m]%%2==1) {c(rep("W", (rep_rate[m]-1)/2), rep("R", (rep_rate[m]-1)/2), sample(c("W","R"), 1))}
          gamts <- c(gamts, gamts_temp1, gamts_temp2, gamts_temp3, gamts_temp4)
        }
        off_dat <- NULL
        
        # if there are 2+ gametes, create offspring
        # add some mutation
        if (length(gamts) > 1 & K_stoch-length(inds$ind_num) > 0) {
          # draw them, limiting by stochastic K
          off_gamts <- sample(gamts, 
                              min(length(gamts), (K_stoch-length(inds$ind_num))*2), 
                              replace=FALSE)
          # only take even number: if odd length, drop the last one
          if (length(off_gamts)%%2==1) off_gamts <- off_gamts[1:length(off_gamts)-1]
          # mutate some
          muts <- rbinom(length(off_gamts), 1, mut_rate) # 1 = does mutate
          off_gamts[which(muts==1)] <- ifelse(off_gamts[which(muts==1)] == "W", "R", "W")
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
          
          # check for any below 0, change to positive
          if (any(off_dat<0)) {
            index_neg <- which(off_dat<0)
            for (i in 1:length(index_neg)) {
              col_num <- ceiling(index_neg[i]/dim(off_dat)[1])
              row_num <- index_neg[i]%%dim(off_dat)[1]
              if (row_num == 0) {row_num <- dim(off_dat)[1]}
              off_dat[row_num,col_num] <- 0
            }
          }
          
          off_dat$mI_pheno = ifelse(off_dat$mI_pheno<0, off_dat$mS_pheno, off_dat$mS_pheno+off_dat$mI_pheno)
          
        }
        
        if (is.null(dim(off_dat))) off_dat <- NULL
      }
      
      # dim(inds)
      # dim(off_dat)
      
      # combine parents and offspring
      inds <- rbind(inds[,1:7], off_dat)
      
      # reindex 
      inds$ind_num <- 1:dim(inds)[1]
      
      # check that northing is negative --> change to zero
      if (any(inds<0)) {
        index_neg <- which(inds<0)
        for (i in 1:length(index_neg)) {
          col_num <- ceiling(index_neg[i]/dim(inds)[1])
          row_num <- index_neg[i]%%dim(inds)[1]
          if (row_num == 0) {row_num <- dim(inds)[1]}
          inds[row_num,col_num] <- 0
        }
      }      
    }
    
    # change individual(s) to I
    inf_inds <- max(1, floor(0.1*dim(inds)[1])) # at least 1 individual, or ~10% of individuals
    inds$inf_stat[sample(1:length(inds$ind_num), inf_inds, replace = F)] <- "I"
    
    # # disease --> allele intro
    # for (r in 1:r_0){
    #   
    #   # these are probabilities of each action for an individual
    #   
    #   # draw K first so you have something to save
    #   K_stoch = floor(rnorm(1, mean=K, sd=K_sd))
    #   
    #   # now disease
    #   # do the three events per generation, times nubmer of cycles per gen
    #   for (disease_cycle in 1:disease_cycles){
    #     
    #     # save things 
    #     S_size <- c(S_size, dim(inds%>%filter(inf_stat=="S"))[1]) # some Ss
    #     I_size <- c(I_size, dim(inds%>%filter(inf_stat=="I"))[1]) # some Is
    #     R_size <- c(R_size, dim(inds%>%filter(inf_stat=="R"))[1]) # some Rs
    #     K_size <- c(K_size, K_stoch) # carrying capacity
    #     r_allele <- c(r_allele, 0) # no r allele
    #     
    #   for (j in 1:length(events)) {
    #    
    #     # make sure someone's alive!
    #     if (dim(inds)[1] == 0) {
    #       extinct_dummy <- TRUE # pop is extinct
    #       S_size <- c(S_size, 0) # no Ss
    #       I_size <- c(I_size, 0) # no Is
    #       R_size <- c(R_size, 0) # no Rs
    #       K_size <- c(K_size, K_stoch) # carrying capacity
    #       r_allele <- c(r_allele, 0) # no R allele currently
    #       break # stop generation loop
    #     }
    #     
    #     # place holders
    #     inds$p_transmit <- rep(NA)
    #     inds$p_survS <- rep(NA)
    #     inds$p_survI <- rep(NA)
    #     inds$p_recovery <- rep(NA)
    #     inds$mortalityS <- rep(NA)
    #     inds$mortalityI <- rep(NA)
    #     inds$change_stat <- rep(NA)
    #     
    #     if (events[j]=="B") {
    #       # get transmission phenotype
    #       if(trans_form == "density") {inds$p_transmit <- 1-exp(-inds$b_pheno*length(which(inds$inf_stat=="I")))}
    #       if(trans_form == "freq") {inds$p_transmit <- 1-exp(-inds$b_pheno)}
    #       # coin flips
    #       inds$change_stat <- vapply(inds$p_transmit, function(x) rbinom(1, 1, x), as.integer(1L))
    #       # note density dependence: length(which(inds$inf_stat=="I"))
    #       tmpS <- inds %>% filter(inf_stat=="S" & change_stat==1)
    #       if (dim(tmpS)[1] > 0) tmpS$inf_stat <- "I" # goes to being an I
    #       # everyone else
    #       tmpI <- inds %>% filter(inf_stat!="S" | (inf_stat=="S" & change_stat==0))
    #       # recombine
    #       inds <- rbind(tmpS, tmpI)
    #       inds <- inds[, 1:7]
    #     }
    #     if (events[j]=="M") {
    #       # probability of survival
    #       inds$p_survS <- exp(-inds$mS_pheno)
    #       inds$p_survI <- exp(-inds$mI_pheno)
    #       # coin flips
    #       inds$mortalityS <- vapply(inds$p_survS, function(x) rbinom(1, 1, x), as.integer(1L))
    #       inds$mortalityI <- vapply(inds$p_survI, function(x) rbinom(1, 1, x), as.integer(1L))
    #       # get mortality phenotype: infects
    #       tmpI <- inds %>% filter(inf_stat == "I") %>% filter(mortalityI==1)
    #       tmpS <- inds %>% filter(inf_stat == "S" | inf_stat == "R") %>% filter(mortalityS==1)
    #       inds <- rbind(tmpI, tmpS)
    #       inds <- inds[, 1:7]
    #     }
    #     if (events[j]=="G" & compartments !="SIX") {
    #       # get recovery phenotype
    #       inds$p_recovery <- (1-exp(-inds$r_pheno))
    #       # coin flip
    #       inds$change_stat <- vapply(inds$p_recovery, function(x) rbinom(1, 1, x), as.integer(1L))
    #       # get I --> not I list
    #       tmpI <- inds %>% filter(inf_stat=="I" & change_stat==1)
    #       if (dim(tmpI)[1] > 0) tmpI$inf_stat <- ifelse(compartments=="SIR", "R", "S") # if it's not SIR, it's SIS
    #       # everyone else
    #       tmpS <- inds %>% filter(inf_stat!="I" | (inf_stat=="I" & change_stat==0))
    #       # recombine
    #       inds <- rbind(tmpS, tmpI)
    #       inds <- inds[, 1:7]
    #     }
    #     
    #   }
    #     
    #   }
    #   
    #   # is everyone dead??
    #   if (dim(inds)[1] == 0) {
    #     extinct_dummy <- TRUE # pop is extinct
    #     S_size <- c(S_size, 0) # no Ss
    #     I_size <- c(I_size, 0) # no Is
    #     R_size <- c(R_size, 0) # no Rs
    #     K_size <- c(K_size, K_stoch) # carrying capacity
    #     r_allele <- c(r_allele, 0) # no R allele currently
    #     break # stop generation loop
    #   }
    #   
    #   # and then reproduce -- still all W's
    #   # how many gametes 
    #   rep_rate <- floor(rnorm(length(inds$ind_num), l_WW, l_sd))
    #   rep_rate[which(rep_rate < 0)] = 0
    #   # convert to gametes
    #   gamts <- rep("W", 2*sum(rep_rate))
    #   # combine gamts into offspring: need 2 gamts and open space
    #   off_dat <- NULL
    #   if (length(gamts) > 1 & K_stoch-length(inds$ind_num) > 0) {
    #     # since all WW type, take the gamts/2 or the remaining space, whichever limits
    #     n_offs <- min(floor(length(gamts)/2), K_stoch-length(inds$ind_num))
    #     # and fill in the data table
    #     off_dat <- data.frame(ind_num = 1:n_offs,
    #                           inf_stat = rep("S", n_offs), # all healthy
    #                           ind_geno = rep("WW", n_offs), # all W type
    #                           b_pheno = c(rnorm(n_offs, mean=b_WW, sd=b_sd)), # transmission placeholder
    #                           mS_pheno = c(rnorm(n_offs, mean=m_SWW, sd=m_Ssd)), # S_WW type mortality
    #                           mI_pheno = c(rnorm(n_offs, mean=m_IWW, sd=m_Isd)), # mortality placeholder
    #                           r_pheno = c(rnorm(n_offs, mean=r_WW, sd=r_sd)) # revovery placeholder 
    #     ) # phenos still in rates
    #   }
    #   
    #   # remove zeros
    #   # check that northing is negative --> change to zero
    #   if (any(off_dat<0)) {
    #     index_neg <- which(off_dat<0)
    #     for (i in 1:length(index_neg)) {
    #       col_num <- ceiling(index_neg[i]/dim(off_dat)[1])
    #       row_num <- index_neg[i]%%dim(off_dat)[1]
    #       if (row_num == 0) {row_num <- dim(off_dat[1])}
    #       off_dat[row_num,col_num] <- 0
    #     }
    #   }
    #   # note mI_pheno is additional to mS_pheno --> address here
    #   off_dat$mI_pheno = off_dat$mS_pheno+off_dat$mI_pheno # m_I is a bonus mortality, should only increase mortality
    #   
    #   # combine parents and offspring
    #   inds <- rbind(inds[,1:7], off_dat)
    #   # reindex 
    #   inds$ind_num <- 1:dim(inds)[1]
    #   
    # }
    
    for (i in 1:ngens) {
      
      # draw probabilities for all process each year --> need to truncate at 0?
      # phenotypes stay the same
      
      K_stoch = floor(rnorm(1, mean=K, sd=K_sd))
      
      # is everyone dead??
      if (dim(inds)[1] == 0) {
        # print(c(i, "all dead should break"))
        extinct_dummy <- TRUE # pop is extinct
        S_size <- c(S_size, 0) # no Ss
        I_size <- c(I_size, 0) # no Is
        R_size <- c(R_size, 0) # no Rs
        K_size <- c(K_size, K_stoch) # carrying capacity
        r_allele <- c(r_allele, 0) # no R allele currently
        break # stop generation loop
      }
      
      # do the three events per generation
      for (disease_cycle in 1:disease_cycles){
        
        r_freq <- (2*length(which(inds$ind_geno == "RR"))+length(which(inds$ind_geno == "WR")))/(2*length(inds$ind_geno == "WW"))
        # save things after each disease cycle -- note this is censusing BEFORE reproduction!
        S_size <- c(S_size, dim(inds%>%filter(inf_stat=="S"))[1]) # some Ss
        I_size <- c(I_size, dim(inds%>%filter(inf_stat=="I"))[1]) # some Is
        R_size <- c(R_size, dim(inds%>%filter(inf_stat=="R"))[1]) # some Rs
        K_size <- c(K_size, K_stoch) # carrying capacity
        r_allele <- c(r_allele, r_freq) # r allele
        
        
      for (j in 1:length(events)) {
          
          # make sure someone's alive!
          if (dim(inds)[1] == 0) {
            extinct_dummy <- TRUE # pop is extinct
            S_size <- c(S_size, 0) # no Ss
            I_size <- c(I_size, 0) # no Is
            R_size <- c(R_size, 0) # no Rs
            K_size <- c(K_size, K_stoch) # carrying capacity
            r_allele <- c(r_allele, 0) # no R allele currently
            break # stop generation loop
          }
          
          # place holders
          inds$p_transmit <- rep(NA)
          inds$p_survS <- rep(NA)
          inds$p_survI <- rep(NA)
          inds$p_recovery <- rep(NA)
          inds$mortalityS <- rep(NA)
          inds$mortalityI <- rep(NA)
          inds$change_stat <- rep(NA)
          
          if (events[j]=="B") {
            # get transmission phenotype
            if(trans_form == "density") {inds$p_transmit <- 1-exp(-inds$b_pheno*length(which(inds$inf_stat=="I")))}
            if(trans_form == "freq") {inds$p_transmit <- 1-exp(-inds$b_pheno)}
            # coin flips
            inds$change_stat <- vapply(inds$p_transmit, function(x) rbinom(1, 1, x), as.integer(1L))
            # note density dependence: length(which(inds$inf_stat=="I"))
            tmpS <- inds %>% filter(inf_stat=="S" & change_stat==1)
            if (dim(tmpS)[1] > 0) tmpS$inf_stat <- "I" # goes to being an I
            # everyone else
            tmpI <- inds %>% filter(inf_stat!="S" | (inf_stat=="S" & change_stat==0))
            # recombine
            inds <- rbind(tmpS, tmpI)
            inds <- inds[, 1:7]
          }
          if (events[j]=="M") {
            # probability of survival
            inds$p_survS <- exp(-inds$mS_pheno)
            inds$p_survI <- exp(-inds$mI_pheno)
            # coin flips
            inds$mortalityS <- vapply(inds$p_survS, function(x) rbinom(1, 1, x), as.integer(1L))
            inds$mortalityI <- vapply(inds$p_survI, function(x) rbinom(1, 1, x), as.integer(1L))
            # get mortality phenotype: infects
            tmpI <- inds %>% filter(inf_stat == "I") %>% filter(mortalityI==1)
            tmpS <- inds %>% filter(inf_stat == "S" | inf_stat == "R") %>% filter(mortalityS==1)
            inds <- rbind(tmpI, tmpS)
            inds <- inds[, 1:7]
          }
          if (events[j]=="G" & compartments !="SIX") {
            # get recovery phenotype
            inds$p_recovery <- (1-exp(-inds$r_pheno))
            # coin flip
            inds$change_stat <- vapply(inds$p_recovery, function(x) rbinom(1, 1, x), as.integer(1L))
            # get I --> not I list
            tmpI <- inds %>% filter(inf_stat=="I" & change_stat==1)
            if (dim(tmpI)[1] > 0) tmpI$inf_stat <- ifelse(compartments=="SIR", "R", "S") # if it's not SIR, it's SIS
            # everyone else
            tmpS <- inds %>% filter(inf_stat!="I" | (inf_stat=="I" & change_stat==0))
            # recombine
            inds <- rbind(tmpS, tmpI)
            inds <- inds[, 1:7]
          }
      }
        
        
      }
      
      if (dim(inds)[1] == 0) {
        # print(c(i, "everyone dead after disease dynamics"))
        extinct_dummy <- TRUE # pop is extinct
        S_size <- c(S_size, 0) # no Ss
        I_size <- c(I_size, 0) # no Is
        R_size <- c(R_size, 0) # no Rs
        K_size <- c(K_size, K_stoch) # carrying capacity
        r_allele <- c(r_allele, 0) # no R allele currently
        break # stop generation loop
      }
      
      # survivors reproduce if there is someone around
      if (dim(inds)[1] > 0) {
      # print(c(i, "dim(inds) > 0"))
      parent_lams <- paste0("l_",inds$ind_geno,sep="")
      rep_rate <- NULL
      for (P in 1:length(parent_lams)) {
           lam_val=ifelse(parent_lams[P]=="l_RR", l_RR, ifelse(parent_lams[P]=="l_WR", l_WR, l_WW))
           rep_rate <- c(rep_rate, floor(rnorm(1, lam_val, l_sd)))
      }
      rep_rate[which(rep_rate < 0)] = 0
      gamts <- NULL
      for (m in 1:length(parent_lams)) {
        # get gametes from parent
        gamts_temp1 <- if (inds$ind_geno[m]=="RR") {rep("R", 2*rep_rate[m])}
        gamts_temp2 <- if (inds$ind_geno[m]=="WW") {rep("W", 2*rep_rate[m])}
        # WRs will need to have coin flip if odd rep_rate
        gamts_temp3 <- if (inds$ind_geno[m]=="WR" & rep_rate[m]%%2==0) {c(rep("W", rep_rate[m]/2), rep("R", rep_rate[m]/2))}
        gamts_temp4 <- if (inds$ind_geno[m]=="WR" & rep_rate[m]%%2==1) {c(rep("W", (rep_rate[m]-1)/2), rep("R", (rep_rate[m]-1)/2), sample(c("W","R"), 1))}
        gamts <- c(gamts, gamts_temp1, gamts_temp2, gamts_temp3, gamts_temp4)
      }
      off_dat <- NULL
      
      # if there are 2+ gametes, create offspring
      # add some mutation
      if (length(gamts) > 1 & K_stoch-length(inds$ind_num) > 0) {
        # draw them, limiting by stochastic K
        off_gamts <- sample(gamts, 
                            min(length(gamts), (K_stoch-length(inds$ind_num))*2), 
                            replace=FALSE)
        # only take even number: if odd length, drop the last one
        if (length(off_gamts)%%2==1) off_gamts <- off_gamts[1:length(off_gamts)-1]
        # mutate some
        muts <- rbinom(length(off_gamts), 1, mut_rate) # 1 = does mutate
        off_gamts[which(muts==1)] <- ifelse(off_gamts[which(muts==1)] == "W", "R", "W")
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
        
        # check for any below 0, change to positive
        if (any(off_dat<0)) {
          index_neg <- which(off_dat<0)
          for (i in 1:length(index_neg)) {
            col_num <- ceiling(index_neg[i]/dim(off_dat)[1])
            row_num <- index_neg[i]%%dim(off_dat)[1]
            if (row_num == 0) {row_num <- dim(off_dat)[1]}
            off_dat[row_num,col_num] <- 0
          }
        }
        
        off_dat$mI_pheno = ifelse(off_dat$mI_pheno<0, off_dat$mS_pheno, off_dat$mS_pheno+off_dat$mI_pheno)
        
      }
      
      if (is.null(dim(off_dat))) off_dat <- NULL
      }
      
      # dim(inds)
      # dim(off_dat)
      
      # combine parents and offspring
      inds <- rbind(inds[,1:7], off_dat)
      # reindex 
      inds$ind_num <- 1:dim(inds)[1]
      
      if (dim(inds)[1] == 0) {
        # print(c(i, "all dead should break"))
        extinct_dummy <- TRUE # pop is extinct
        S_size <- c(S_size, 0) # no Ss
        I_size <- c(I_size, 0) # no Is
        R_size <- c(R_size, 0) # no Rs
        K_size <- c(K_size, K_stoch) # carrying capacity
        r_allele <- c(r_allele, 0) # no R allele currently
        break # stop generation loop
      }
      
    }
    
    output_dat <- c(# things to decide if ER occured
                    extinct <- extinct_dummy, # did the population go extinct? T/F
                    pop_drop20 <- any(S_size[(d_0):length(S_size)]+I_size[(d_0):length(I_size)]+R_size[(d_0):length(I_size)] < K_size[(d_0):length(I_size)]*0.20), # did the population drop? T/F
                    pop_drop50 <- any(S_size[(d_0):length(S_size)]+I_size[(d_0):length(I_size)]+R_size[(d_0):length(I_size)] < K_size[(d_0):length(I_size)]*0.50), # did the population drop? T/F
                    pop_drop80 <- any(S_size[(d_0):length(S_size)]+I_size[(d_0):length(I_size)]+R_size[(d_0):length(I_size)] < K_size[(d_0):length(I_size)]*0.80), # did the population drop? T/F
                    r_allele_peak15 <- any(r_allele[(d_0):length(S_size)] > 0.15), # did the allele spread at any point?
                    r_allele_peak45 <- any(r_allele[(d_0):length(S_size)] > 0.45), # did the allele spread at any point?
                    r_allele_peak75 <- any(r_allele[(d_0):length(S_size)] > 0.75), # did the allele spread at any point?
                    # pop gen outcomes
                    final_r_allele <- mean(tail(r_allele, 15), na.rm = TRUE), # average r allele frequency by end of simulation
                    final_pop_size <- mean(tail(S_size+I_size+R_size, 15), na.rm = TRUE),
                    # disease outcomes
                    final_inf_prev <- mean(tail(I_size/(S_size+I_size+R_size), 15), na.rm = TRUE), # helps determine if I was lost
                    # get the extremes
                    max_r_allele <- max(r_allele[d_0:length(r_allele)], na.rm = TRUE), # only consider post-disease hits
                    time_max_r_allele <- which(r_allele==max(r_allele))[1],
                    max_inf_prev <- max(I_size/(S_size+I_size+R_size), na.rm = TRUE),
                    time_last_zero_inf <- sum(c(d_0,  which(I_size[d_0+1:length(I_size)]/(S_size[d_0+1:length(S_size)]+I_size[d_0+1:length(I_size)]+R_size[d_0+1:length(I_size)])==0)[1]), na.rm = T), # only care about post-disease
                    min_pop <- min(S_size[(d_0):length(S_size)]+I_size[(d_0):length(I_size)]+R_size[(d_0):length(I_size)], na.rm = TRUE), # only care about post-disease min size
                    time_min_pop <- which(S_size+I_size+R_size == min_pop)[1], # time that min pop was hit
                    # did recvoery happen AFTER min_pop?
                    at_K95 <- ifelse(extinct_dummy, FALSE, any(0.95*K_size[time_min_pop:length(K_size)]<=S_size[time_min_pop:length(S_size)]+I_size[time_min_pop:length(I_size)]+R_size[time_min_pop:length(R_size)])),  # did we get back up? note that if extinct, will use start of time series (not good!)
                    first_K95 <- time_min_pop + which(S_size[time_min_pop:length(S_size)]+I_size[time_min_pop:length(I_size)]+R_size[time_min_pop:length(R_size)]>=0.95*K_size[time_min_pop:length(K_size)])[1], 
                    # some initial metrics re: mutation-selection
                    r_ts_d0 <- r_allele[d_0]
    )
    
  }
  
  # return summary statistics
  all_reps_dat <- as.data.frame(all_reps)
  colnames(all_reps_dat) <- c(# ER determinants
                              "extinct", 
                              "pop_drop20", "pop_drop50", "pop_drop80",
                              "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                              # pop gen outcomes
                              "final_r_allele", "final_pop_size",
                              # disease outcomes
                              "final_inf_prev",
                              # get the extremes -- not currenlty outputted, but could helpful later???
                              "max_r_allele", "time_max_r_allele",
                              "max_inf_prev", "time_last_zero_inf",
                              "min_pop", "time_min_pop",
                              "at_K95", "firstK95", "init_R")
  
  # for the early tests, want to know the summary from each simulation
  return(all_reps_dat) 

}
