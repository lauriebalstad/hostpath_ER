library(parallel) #for using multiple cores (mclapply for mac and linux, parLapply for windows)
# library(randomForest)

# load function
source("funcs/00_model_cluster.R")
# source("funcs/00_model_nodensity.R")
# load data

# base_vect <- c(2, 4, # SIS, BGM -- can comp with 4 for MBG
#                -1, -1, -1, 0.05, # 3-6 no enviro dep transmission
#                -1, -1, -1, 0.05, # 7-10 always dens dep transmission
#                0.05, 0.05, 0.05, # 11-13 + 17 background mort 
#                1, 1, 1, # 14-16 + 18 disease mort
#                0.05, 0.05, # 17-18 mort sd
#                0.05, 0.05, 0.05, 0.05, # 19-22 recovery, unused bc SIX (1)
#                0.5, 0.75, 1, 0.005, # 23-27 reproduction & mutation
#                # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
#                100, 4, # 28-29 carrying capacity things
#                120, 3, # 30-32 timing things # note change from d = 1 to d = 2???
#                0.05, 0.1, 160) # init R and init disease, ngens

base_vect <- c(2, 2, # SIS, BGM -- can comp with 4 for MBG
               0.073/3, 0.073/3, 0.073/3, 0.001, # 3-6 no enviro dep transmission
               0.073, 0.073, 0.073, 0.01, # 7-10 always dens dep transmission
               0.0002, 0.0002, 0.0002, # 11-13 + 17 background mort 
               0.5, 0.5, 0.5, # 14-16 + 18 disease mort
               0.00005, 0.05, # 17-18 mort sd
               0.005, 0.005, 0.005, 0.0005, # 19-22 recovery, unused bc SIX (1)
               exp(0.0866)*0.6, exp(0.0866)*0.8, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
               # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
               200, 4, # 28-29 carrying capacity things
               120, 5, # 30-32 timing things # note change from d = 1 to d = 2???
               0.01, 0.05, 600) # init R and init disease, ngens

nden_vect <- c(2, 2, # SIS, BGM -- can comp with 4 for MBG
               Inf, Inf, Inf, 0.001, # 3-6 no enviro dep transmission
               Inf, Inf, Inf, 0.01, # 7-10 always dens dep transmission
               0.0002, 0.0002, 0.0002, # 11-13 + 17 background mort 
               0.5, 0.5, 0.5, # 14-16 + 18 disease mort
               0.00005, 0.05, # 17-18 mort sd
               -2, -2, -2, 0.0005, # 19-22 recovery, unused bc SIX (1)
               exp(0.0866)*0.6, exp(0.0866)*0.8, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
               # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
               200, 4, # 28-29 carrying capacity things
               120, 5, # 30-32 timing things # note change from d = 1 to d = 2???
               0.01, 1, 600) # init R and init disease, ngens

#-----FIG2: COST/BENEFIT-----
# comparing different cost/benefits
# sis_cb_A <- expand.grid("benefit" = c(0, 0.25, 0.5, 1)*0.073, # for percent change in benefit: 0, 50, 75, 100
#                       "cost" = c(0.1, 0.5, 1), # 0.1 = some cost, 1 = no cost
#                       "case" = c("\u03b2"), # tranmission
#                       "compartments" = c(2))
# sis_cb_B <- expand.grid("benefit" = c(30, 10, 3, 1)*0.005, # for percent change in benefit: 0, 100, 500, 1500
#                         "cost" = c(0.1, 0.5, 1), # 0.1 = some cost, 1 = no cost
#                         "case" = c("\u03d2"), # recovery
#                         "compartments" = c(2))
sis_cb <- expand.grid("benefit" = c(0, 0.25, 0.5, 1)*0.5, # for percent change in benefit: 0, 50, 75, 100
                      "cost" = c(0.2, 0.6, 1), # 0.1 = some cost, 1 = no cost
                      # "case" = c("\u03bc"), # mortality
                      # "compartments" = c(2), 
                      "transmission" = 1:3)
# sis_cb <- rbind(sis_cb_A, sis_cb_B, sis_cb_C)
N <- dim(sis_cb)[1]
sis_cb$number <- 1:N
rep_num <- 2500
rep_sis_cb <- matrix(rep(t(sis_cb), rep_num), ncol = ncol(sis_cb), byrow = T)
rep_sis_cb <- as.data.frame(rep_sis_cb)
colnames(rep_sis_cb) <- colnames(sis_cb)

sim_results_cb <- mclapply(1:(N*rep_num), function(i){

  
  if (rep_sis_cb$transmission[i] == 3) {
    
    cb_vect <- nden_vect
    cb_vect[34] <- as.numeric(rep_sis_cb$`number`[i]) # index

    # benefits update
    cb_vect[14] <- as.numeric(rep_sis_cb$`benefit`[i]) # sis_cb$benefit[cb]*cb_vect[15]
    cb_vect[15] <- (cb_vect[14]+cb_vect[16])/2
    
    # cost update
    cb_vect[23] <- as.numeric(rep_sis_cb$cost[i])*exp(0.0866)
    cb_vect[24] <- (cb_vect[23]+cb_vect[25])/2
    
  }
  
  if (rep_sis_cb$transmission[i] %in% c(1, 2)) {
    
    cb_vect <- base_vect
    cb_vect[34] <- as.numeric(rep_sis_cb$`number`[i]) # index
    
    if (rep_sis_cb$transmission[i]==1) {
      cb_vect[3:5] <- -1 # no outside infection
      cb_vect[7:9] <- 0.073 # transmission increased (dens)
    }
    if (rep_sis_cb$transmission[i]==2) {
      cb_vect[3:5] <- 0.073/3 # both types of transmission
      cb_vect[7:9] <- 0.073 # transmission increased (dens)
    }
    
    # benefits update
    cb_vect[14] <- as.numeric(rep_sis_cb$`benefit`[i]) # sis_cb$benefit[cb]*cb_vect[15]
    cb_vect[15] <- (cb_vect[14]+cb_vect[16])/2
    
    # cost update
    cb_vect[23] <- as.numeric(rep_sis_cb$cost[i])*exp(0.0866)
    cb_vect[24] <- (cb_vect[23]+cb_vect[25])/2
    
  }
  
  # cb_vect <- base_vect
  # cb_vect[34] <- as.numeric(rep_sis_cb$`number`[i])
  # cb_vect[1] <- as.numeric(rep_sis_cb$`compartment`[i]) # SIS or SIR
  # # benefit update
  # if (rep_sis_cb$case[i] == "\u03b2") {
  #   cb_vect[7] <- as.numeric(rep_sis_cb$`benefit`[i]) # sis_cb$benefit[cb]*cb_vect[8]
  #   cb_vect[8] <- (cb_vect[7]+cb_vect[9])/2
  #   cb_vect[3] <- as.numeric(rep_sis_cb$`benefit`[i])/3 # and also outside likelihood
  #   cb_vect[4] <- (cb_vect[3]+cb_vect[5])/2}
  # if (rep_sis_cb$case[i] == "\u03bc") {
  #   cb_vect[14] <- as.numeric(rep_sis_cb$`benefit`[i]) # sis_cb$benefit[cb]*cb_vect[15]
  #   cb_vect[15] <- (cb_vect[14]+cb_vect[16])/2}
  # if (rep_sis_cb$case[i] == "\u03d2") {
  #   cb_vect[19] <- as.numeric(rep_sis_cb$`benefit`[i]) # sis_cb$benefit[cb]*cb_vect[20]
  #   cb_vect[20] <- (cb_vect[19]+cb_vect[21])/2}
  # # note for G, the meaning is opposite
  # 
  # # cost update -- not dependent on the case
  # cb_vect[23] <- as.numeric(rep_sis_cb$cost[i])*exp(0.0866)
  # cb_vect[24] <- (cb_vect[23]+cb_vect[25])/2

  # run simulation
  suppressWarnings(cluster_run(cb_vect))

}, mc.cores = parallel::detectCores())

# save things!
saveRDS(sim_results_cb, file = paste0("dat/cb_fig_", format(Sys.time(), "%m%d"), ".Rdata"))

#-----FIG1: STRUCTURALS-----
# so.... fig 1 will be comparision of structure across base model
cases <- expand.grid("transmission" = 1:3,
                     "compartments" = 1:3, # note SIR doesn't seem to have a big effect
                     "robustness" = 1:4) # 1 = mortality, 2 = transmission, 3 = recovery, 4 = demographic rescue only
# remove 1/3 compartment/robustness combo
cases <- cases[c(1:18, 22:36), ]
N <- dim(cases)[1]
cases$number <- 1:N
rep_num <- 2500 # 1 = testing run time, 400 = 5 tasks
rep_cases <- matrix(rep(t(cases), rep_num), ncol = ncol(cases), byrow = T)
rep_cases <- as.data.frame(rep_cases)
colnames(rep_cases) <- c("transmission", "compartments", "robustness", "number")

sim_results_str <- mclapply(1:(N*rep_num), function(i){

  if (rep_cases$transmission[i] == 3) {
    
    case_vect <- nden_vect
    case_vect[34] <- rep_cases$`number`[i]
    
    # mortality
    if (rep_cases$robustness[i] == 1) {
      case_vect[14] <- 0
      case_vect[15] <- (case_vect[14]+case_vect[16])/2
    }
    # transmission
    if (rep_cases$robustness[i] == 2) {
      case_vect[7:9] <- c(0, -log(1-0.5), Inf) # probability of infection for WR should be 0.5ish
      case_vect[3:5] <- c(0, -log(1-0.5), Inf) # probability of infection for WR should be 0.5ish
    } 
    # recovery
    if (rep_cases$robustness[i] == 3) {
      case_vect[19] <- 0.005*20 # about doubles probability of recovery compared to regular case?
      case_vect[20] <- -log(1-mean(c(1, 1-exp(-case_vect[19])))) # probability of recovery for WR should be half of RR and WW approx.
    } # 1 = recovery, nb: recovery higher
    
    # modify base parameters -- compartment
    case_vect[1] <- rep_cases$compartments[i]
    
  }
  
  if (rep_cases$transmission[i] %in% c(1, 2)) {
    
    case_vect <- base_vect
    case_vect[34] <- rep_cases$`number`[i]
    
    # deal w transmission first
    # transmission
    if (rep_cases$transmission[i]==1) {
      case_vect[3:5] <- -1 # no outside infection
      case_vect[7:9] <- 0.073 # transmission increased (dens)
      if (rep_cases$robustness[i] == 2) {case_vect[7:9] <- c(0, 0.073/2, 0.073)}
    }
    if (rep_cases$transmission[i]==2) {
      case_vect[3:5] <- 0.073/3 # both types of transmission
      case_vect[7:9] <- 0.073 # transmission increased (dens)
      if (rep_cases$robustness[i] == 2) {
        case_vect[3:5] <- c(0, 0.073/3/2, 0.073/3)
        case_vect[7:9] <- c(0, 0.073/2, 0.073)
        }
    }
    
    # mortality
    if (rep_cases$robustness[i] == 1) {
      case_vect[14] <- 0
      case_vect[15] <- (case_vect[14]+case_vect[16])/2
    }
    # recovery
    if (rep_cases$robustness[i] == 3) {
      case_vect[19] <- 0.005*20 # about doubles probability of recovery compared to regular case?
      case_vect[20] <- (case_vect[19] + case_vect[21])/2 # probability of recovery for WR should be half of RR and WW approx.
    } # 1 = recovery, nb: recovery higher
    
    # modify base parameters -- compartment
    case_vect[1] <- rep_cases$compartments[i]
    
  }
  
  # case_vect <- base_vect
  # case_vect[34] <- rep_cases$`number`[i]
  # # set disease transmission
  # if (rep_cases$`transmission type`[i]==1) {
  #   case_vect[3:5] <- -1 # no outside infection
  #   case_vect[7:9] <- 0.073 # transmission increased (dens)
  # }
  # if (rep_cases$`transmission type`[i]==2) {
  #   case_vect[3:5] <- 0.073/3 # both types of transmission
  #   case_vect[7:9] <- 0.073 # transmission increased (dens)
  # }
  # # mortality
  # if (rep_cases$robustness[i] == 1) {
  #   case_vect[14] <- 0
  #   case_vect[15] <- (case_vect[14]+case_vect[16])/2
  # } # 1 = transmission
  # if (rep_cases$robustness[i] == 2) {
  #   if (rep_cases$`transmission type`[i]==1) {case_vect[7:9] <- c(0, 0.073/2, 0.073)}
  #   if (rep_cases$`transmission type`[i]==2) {case_vect[3:5] <- c(0, 0.073/3/2, 0.073/3); 
  #                                             case_vect[7:9] <- c(0, 0.073/2, 0.073)}
  # } # 1 = recovery
  # if (rep_cases$robustness[i] == 3) {
  #   case_vect[19] <- 0.005*3
  #   case_vect[20] <- (case_vect[19]+case_vect[21])/2
  # } # 1 = recovery, nb: recovery higher
  # 
  # # modify base parameters -- compartment
  # case_vect[1] <- rep_cases$compartments[i]

  suppressWarnings(cluster_run(case_vect)) # gets upset about some of the NaNs in rbinom, but no obvious issue with actual calcualtions

}, mc.cores = parallel::detectCores())

# save things!
saveRDS(sim_results_str, file = paste0("dat/str_fig_", format(Sys.time(), "%m%d"), ".Rdata"))

#-----FIG3: POP SIZES-----
# comparing pop size/robustness/mutation rate
sis_pop_mut <- expand.grid("population size" = c(50, 100, 500),
                           "robustness" = 1, # M, G only bc they are strongest ER from past exploration
                           # "mutation rate" = c(0.0005, 0.005, 0.01),
                           "compartments" = 2,
                           "transmission type" = 1:3)
N <- dim(sis_pop_mut)[1]
sis_pop_mut$number <- 1:N
rep_num <- 2500
rep_sis_pm <- matrix(rep(t(sis_pop_mut), rep_num), ncol = ncol(sis_pop_mut), byrow = T)
rep_sis_pm <- as.data.frame(rep_sis_pm)
colnames(rep_sis_pm) <- colnames(sis_pop_mut)

sim_results_pm <- mclapply(1:(N*rep_num), function(i){

  # case_vect <- base_vect
  # case_vect[34] <- rep_sis_pm$number[i]
  # case_vect[1] <- rep_sis_pm$compartments[i]
  # 
  # if (rep_sis_pm$`transmission type`[i]==1) {
  #   case_vect[7:9] <- 1.25 # transmission increased (dens)
  # }
  # if (rep_sis_pm$`transmission type`[i]==2) {
  #   case_vect[3:5] <- 2 # transmission increased (freq)
  # }
  # 
  # case_vect[27] <- rep_sis_pm$`population size`[i]
  # case_vect[26] <- rep_sis_pm$`mutation rate`[i]
  # 
  # if (rep_sis_pm$robustness[i] == 1) {
  #   case_vect[14] <- 0
  #   case_vect[15] <- 0.5*case_vect[16]
  # } # 1 = mortality
  # if (rep_sis_pm$robustness[i] == 2) {
  #   case_vect[7] <- 0*case_vect[9]
  #   case_vect[8] <- 0.5*case_vect[9]
  # } # 1 = transmission
  # if (rep_sis_pm$robustness[i] == 3) {
  #   case_vect[19] <- 1.5
  #   case_vect[20] <- (case_vect[19]+case_vect[21])/2
  # } # 1 = recovery, nb: recovery higher
  # 
  # # use base case cost for consistency
  # 
  # # run simulation
  # cluster_run(case_vect)
  
  case_vect[27] <- rep_sis_pm$`population size`[i] # change population size--always true
  
  if (rep_sis_pm$transmission[i] == 3) {
    
    case_vect <- nden_vect
    case_vect[34] <- rep_sis_pm$`number`[i]
    
    # mortality
    if (rep_sis_pm$robustness[i] == 1) {
      case_vect[14] <- 0
      case_vect[15] <- (case_vect[14]+case_vect[16])/2
    }
    # transmission
    if (rep_sis_pm$robustness[i] == 2) {
      case_vect[7:9] <- c(0, -log(1-0.5), Inf) # probability of infection for WR should be 0.5ish
      case_vect[3:5] <- c(0, -log(1-0.5), Inf) # probability of infection for WR should be 0.5ish
    } 
    # recovery
    if (rep_sis_pm$robustness[i] == 3) {
      case_vect[19] <- 0.005*20 # about doubles probability of recovery compared to regular case?
      case_vect[20] <- -log(1-mean(c(1, 1-exp(-case_vect[19])))) # probability of recovery for WR should be half of RR and WW approx.
    } # 1 = recovery, nb: recovery higher
    
    # modify base parameters -- compartment
    case_vect[1] <- rep_sis_pm$compartments[i]
    
  }
  
  if (rep_sis_pm$transmission[i] %in% c(1, 2)) {
    
    case_vect <- base_vect
    case_vect[34] <- rep_sis_pm$`number`[i]
    
    # deal w transmission first
    # transmission
    if (rep_sis_pm$transmission[i]==1) {
      case_vect[3:5] <- -1 # no outside infection
      case_vect[7:9] <- 0.073 # transmission increased (dens)
      if (rep_sis_pm$robustness[i] == 2) {case_vect[7:9] <- c(0, 0.073/2, 0.073)}
    }
    if (rep_sis_pm$transmission[i]==2) {
      case_vect[3:5] <- 0.073/3 # both types of transmission
      case_vect[7:9] <- 0.073 # transmission increased (dens)
      if (rep_sis_pm$robustness[i] == 2) {
        case_vect[3:5] <- c(0, 0.073/3/2, 0.073/3)
        case_vect[7:9] <- c(0, 0.073/2, 0.073)
      }
    }
    
    # mortality
    if (rep_sis_pm$robustness[i] == 1) {
      case_vect[14] <- 0
      case_vect[15] <- (case_vect[14]+case_vect[16])/2
    }
    # recovery
    if (rep_sis_pm$robustness[i] == 3) {
      case_vect[19] <- 0.005*20 # about doubles probability of recovery compared to regular case?
      case_vect[20] <- (case_vect[19] + case_vect[21])/2 # probability of recovery for WR should be half of RR and WW approx.
    } # 1 = recovery, nb: recovery higher
    
    # modify base parameters -- compartment
    case_vect[1] <- rep_sis_pm$compartments[i]
    
  }
  
  suppressWarnings(cluster_run(case_vect)) # gets upset about some of the NaNs in rbinom, but no obvious issue with actual calcualtions

}, mc.cores = parallel::detectCores())

# save things!
saveRDS(sim_results_pm, file = paste0("dat/pm_fig_", format(Sys.time(), "%m%d"), ".Rdata"))

#-----SUPP FIG: DISEASE CYCLES-----
# other variable to check in the number of disease cycles between reproductions
sis_dc <- expand.grid("disease cycles" = c(2, 4, 6), # aiming to span a bit of a range wihtout too many simulations
                      "robustness" = c(1, 2), # M, B only 
                      "event order" = c(1, 4, 6),
                      "transmission" = 1:3,
                      "compartments" = c(2)) # 1 = density, 2 = envrionemntal
N <- dim(sis_dc)[1]
sis_dc$number <- 1:N
rep_num <- 2500
rep_sis_dc <- matrix(rep(t(sis_dc), rep_num), ncol = ncol(sis_dc), byrow = T)
rep_sis_dc <- as.data.frame(rep_sis_dc)
colnames(rep_sis_dc) <- colnames(sis_dc)

sim_results_dc <- mclapply(1:(N*rep_num), function(i){

  case_vect <- base_vect
  case_vect[34] <- rep_sis_dc$number[i]

  # if (rep_sis_dc$`transmission`[i] == 2) {
  #   case_vect[3:5] <- 2
  # }
  # 
  # if (rep_sis_dc$`transmission`[i] == 1) {
  #   case_vect[7:9] <- 1.25
  # }
  # 
  # if (rep_sis_dc$robustness[i] == 1) {
  #   case_vect[14] <- 0
  #   case_vect[15] <- case_vect[16]*0.5
  # } # 1 = mortality
  # if (rep_sis_dc$robustness[i] == 3) {
  #   case_vect[19] <- 1.5
  #   case_vect[20] <- (case_vect[19]+case_vect[21])/2
  # } # 3 = recovery, nb: recovery higher
  
  if (rep_sis_dc$transmission[i] == 3) {
    
    case_vect <- nden_vect
    case_vect[34] <- rep_sis_dc$`number`[i]
    
    # mortality
    if (rep_sis_dc$robustness[i] == 1) {
      case_vect[14] <- 0
      case_vect[15] <- (case_vect[14]+case_vect[16])/2
    }
    # transmission
    if (rep_sis_dc$robustness[i] == 2) {
      case_vect[7:9] <- c(0, -log(1-0.5), Inf) # probability of infection for WR should be 0.5ish
      case_vect[3:5] <- c(0, -log(1-0.5), Inf) # probability of infection for WR should be 0.5ish
    } 
    # recovery
    if (rep_sis_dc$robustness[i] == 3) {
      case_vect[19] <- 0.005*20 # about doubles probability of recovery compared to regular case?
      case_vect[20] <- -log(1-mean(c(1, 1-exp(-case_vect[19])))) # probability of recovery for WR should be half of RR and WW approx.
    } # 1 = recovery, nb: recovery higher
    
    # modify base parameters -- compartment
    case_vect[1] <- rep_sis_dc$compartments[i]
    
  }
  
  if (rep_sis_dc$transmission[i] %in% c(1, 2)) {
    
    case_vect <- base_vect
    case_vect[34] <- rep_sis_dc$`number`[i]
    
    # deal w transmission first
    # transmission
    if (rep_sis_dc$transmission[i]==1) {
      case_vect[3:5] <- -1 # no outside infection
      case_vect[7:9] <- 0.073 # transmission increased (dens)
      if (rep_sis_dc$robustness[i] == 2) {case_vect[7:9] <- c(0, 0.073/2, 0.073)}
    }
    if (rep_sis_dc$transmission[i]==2) {
      case_vect[3:5] <- 0.073/3 # both types of transmission
      case_vect[7:9] <- 0.073 # transmission increased (dens)
      if (rep_sis_dc$robustness[i] == 2) {
        case_vect[3:5] <- c(0, 0.073/3/2, 0.073/3)
        case_vect[7:9] <- c(0, 0.073/2, 0.073)
      }
    }
    
    # mortality
    if (rep_sis_dc$robustness[i] == 1) {
      case_vect[14] <- 0
      case_vect[15] <- (case_vect[14]+case_vect[16])/2
    }
    # recovery
    if (rep_sis_dc$robustness[i] == 3) {
      case_vect[19] <- 0.005*20 # about doubles probability of recovery compared to regular case?
      case_vect[20] <- (case_vect[19] + case_vect[21])/2 # probability of recovery for WR should be half of RR and WW approx.
    } # 1 = recovery, nb: recovery higher
    
    # modify base parameters -- compartment
    case_vect[1] <- rep_sis_dc$compartments[i]
    
  }
  
  # modify base parameters
  case_vect[2] <- rep_sis_dc$`event order`[i] # diddo
  case_vect[30] <- rep_sis_dc$`disease cycles`[i]

  # run simulation
  suppressWarnings(cluster_run(case_vect))
  
}, mc.cores = parallel::detectCores())

# save things!
saveRDS(sim_results_dc, file = paste0("dat/dc_fig_", format(Sys.time(), "%m%d"), ".Rdata"))
