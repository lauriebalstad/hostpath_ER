base_vect <- c(2, 2, 1, # SIS, BGM, density
               100, 0, # timing things -- no wait between mutation and non-mutation
               1, 1, 1, 0.05, # 6-9: transmission
               0.05, 0.05, 0.05, # background mort + 16
               1, 1, 1, # 13-15 + 17: disease mort
               0.01, 0.02, # mort sd
               0.05, 0.05, 0.05, 0.01, # 18-21: recovery, unused bc SIX (1)
               1.5, 1.75, 2, 0.2, 0.005, # reproduction & mutation -- try smaller number
               100, 4, # carrying capacity things 
               1) # nubmer of disease cycles per gen -- tyring 2
# nb: in base case, no benefit of allele, but cost is fixed
# nb: disease ngens = 150

#-----FIG1: STRUCTURALS-----
# so.... fig 1 will be comparision of structure across base model

cases <- expand.grid("transmission type" = 1:2, 
                     "compartments" = 1:3, # note SIR doesn't seem to have a big effect
                     "robustness" = 1:4) # 1 = mortality, 2 = transmission, 3 = recovery, 4 = demographic rescue only
# remove 1/3 compartment/robustness combo
cases <- cases[c(1:12, 15:24), ]
sim_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores()-1)
doParallel::registerDoParallel(cl)

for (case in 19:22){
  
  print(case) # for error id and to check in
  
  case_vect <- base_vect
  # if transmission type is frequency (2), need to increase baseline transmission a bit
  if (cases$`transmission type`[case] == 2) {
    case_vect[6:8] <- 3 # transmission increased
  }
  if (cases$robustness[case] == 1) {
    case_vect[13] <- 0
    case_vect[14] <- 0.5
  } # 1 = mortality
  if (cases$robustness[case] == 2) {
    case_vect[6] <- 0*case_vect[8] # since will be scaled by WW rate
    case_vect[7] <- 0.5*case_vect[8]
  } # 1 = transmission
  if (cases$robustness[case] == 3) {
    case_vect[18] <- 1
    case_vect[19] <- 0.5
  } # 1 = recovery, nb: recovery higher
  
  # modify base parameters
  case_vect[3] <- cases$`transmission type`[case] # since numbers have meaning
  case_vect[1] <- cases$compartments[case] # diddo
  
  # run simulation
  case_dat <- run_gens(case_vect, 150, 500)
  
  # save connect to case info
  case_dat$`transmission type` <- rep(cases$`transmission type`[case])
  case_dat$compartments <- rep(cases$compartments[case])
  case_dat$robustness <- rep(cases$robustness[case])
  
  # then append to sim_dat
  sim_dat <- rbind(sim_dat, case_dat)
  
}

# end parallel cluster
stopCluster(cl)

str_dat <- sim_dat # G = 0.05

# save things!
saveRDS(str_dat, file = "figures/figure_data/str_dat.rds")

#-----FIG2: COST/BENEFIT-----
# comparing different cost/benefits
# so.... fig 1 will be comparision of structure across base model
sis_cb <- expand.grid("benefit" = c(0, 0.25, 0.5, 1), # 0 = total benefit, 1.5 = some cost for B & M, but opposite for G
                      "cost" = c(1.5, 2), # 1 = some cost, 2 = no cost
                      "case" = c("\u03b2", "\u03bc", "\u03d2")) 
cb_dat <- NULL # this is where everything will get stored
# use fig2_dat as cb_dat, it currently has cb 1:10

# doing parallel
cl <- parallel::makeCluster(detectCores()-1)
doParallel::registerDoParallel(cl)

for (cb in 21:30){
  
  print(cb) # for error id and to check in
  
  cb_vect <- base_vect
  cb_vect[1] <- 2 # SIS, with density
  # benefit update
  if (sis_cb$case[cb] == "\u03b2") {
    cb_vect[6] <- sis_cb$benefit[cb] # sis_cb$benefit[cb]*cb_vect[8]
    cb_vect[7] <- (cb_vect[6]+cb_vect[8])/2}
  if (sis_cb$case[cb] == "\u03bc") { 
    cb_vect[13] <- sis_cb$benefit[cb] # sis_cb$benefit[cb]*cb_vect[15]
    cb_vect[14] <- (cb_vect[13]+cb_vect[15])/2}
  if (sis_cb$case[cb] == "\u03d2") { 
    cb_vect[18] <- sis_cb$benefit[cb] # sis_cb$benefit[cb]*cb_vect[20]
    cb_vect[19] <- (cb_vect[18]+cb_vect[20])/2}
  # note for G, the meaning is opposite
  
  # cost update -- not dependent on the case
  cb_vect[22] <- sis_cb$cost[cb]
  cb_vect[23] <- (cb_vect[22]+cb_vect[24])/2
  
  # run simulation
  case_dat <- run_gens(cb_vect, 150, 500)
  
  # save connect to case info
  case_dat$benefit <- rep(sis_cb$benefit[cb])
  case_dat$cost <- rep(sis_cb$cost[cb])
  case_dat$case <- rep(sis_cb$case[cb])
  
  # then append to sim_dat
  cb_dat <- rbind(cb_dat, case_dat)
  
}

# end parallel cluster
stopCluster(cl)

cost_ben_dat <- cb_dat # G = 0.05

# save things!
saveRDS(cost_ben_dat, file = "figures/figure_data/cost_ben_dat.rds")

#-----FIG3: POP SIZES-----
# comparing pop size/robustness/mutation rate
sis_pop_mut <- expand.grid("population size" = c(50, 100, 500), # probably provides enough insight...
                           "robustness" = c(1, 2, 3), # M, G only bc they are strongest ER (from fig 3)
                           "mutation rate" = c(0.0005, 0.005, 0.01)) # just for variation 
pm_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)
# load pm data from fig4, need to run 13:24 still
# right now, only enviro transmission result

for (case in 20:27){
  
  print(case) # for error id and to check in
  
  case_vect <- base_vect
  case_vect[1] <- 2 # SIS
  case_vect[27] <- sis_pop_mut$`population size`[case]
  case_vect[26] <- sis_pop_mut$`mutation rate`[case]
  
  if (sis_pop_mut$robustness[case] == 1) {
    case_vect[13] <- 0
    case_vect[14] <- 0.5
  } # 1 = mortality
  if (sis_pop_mut$robustness[case] == 2) {
    case_vect[6] <- 0*case_vect[8] # since will be scaled by WW rate
    case_vect[7] <- 0.5*case_vect[8]
  } # 1 = transmission
  if (sis_pop_mut$robustness[case] == 3) {
    case_vect[18] <- 1
    case_vect[19] <- 0.5
  } # 1 = recovery, nb: recovery higher
  
  # modify base parameters
  case_vect[3] <- 1 # all density
  
  # use base case cost for consistency
  
  # run simulation
  case_dat <- run_gens(case_vect, 150, 500)
  
  # save connect to case info
  case_dat$`population size` <- rep(sis_pop_mut$`population size`[case])
  case_dat$`mutation rate` <- rep(sis_pop_mut$`mutation rate`[case])
  case_dat$`robustness` <- rep(sis_pop_mut$`robustness`[case])
  
  # then append to sim_dat
  pm_dat <- rbind(pm_dat, case_dat)
  
}

# end parallel cluster
stopCluster(cl)

ext_risk_dat <- pm_dat # G = 0.05

# save things!
saveRDS(ext_risk_dat, file = "figures/figure_data/ext_risk_dat.rds")

#-----SUPP FIG: DISEASE CYCLES-----
# other variable to check in the number of disease cycles between reproductions
sis_dc <- expand.grid("disease cycles" = 1:2, # probably provides enough insight...
                      "robustness" = c(1, 3), # M, G only bc they are strongest ER (from fig 3)
                      "event order" = c(2, 4, 6), # just for variation 
                      "transmission" = 1) # prelim suggests transmission types are similar
dc_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

for (case in 1:12){
  
  print(case) # for error id and to check in
  
  case_vect <- base_vect
  case_vect[1] <- 2 # SIS
  
  # if transmission type is frequency (2), need to increase baseline transmission a bit
  if (sis_dc$`transmission`[case] == 2) {
    case_vect[6:8] <- 2 # transmission increased
  }
  
  if (sis_dc$robustness[case] == 1) {
    case_vect[13] <- 0
    case_vect[14] <- 0.5
  } # 1 = mortality
  if (sis_dc$robustness[case] == 2) {
    case_vect[6] <- 0*case_vect[8] # since will be scaled by WW rate
    case_vect[7] <- 0.5*case_vect[8]
  } # 1 = transmission
  if (sis_dc$robustness[case] == 3) {
    case_vect[18] <- 1
    case_vect[19] <- 0.5
  } # 1 = recovery, nb: recovery higher
  
  # modify base parameters
  case_vect[3] <- sis_dc$`transmission`[case] # since numbers have meaning
  case_vect[2] <- sis_dc$`event order`[case] # diddo
  case_vect[29] <- sis_dc$`disease cycles`[case]
  
  # # cost update -- not dependent on the case and always no cost
  # case_vect[22] <- 2
  # case_vect[23] <- (case_vect[22]+case_vect[24])/2
  # use base case cost for consistency
  
  # cycles update: equal time steps
  case_vect[4] <- floor(100/sis_dc$`disease cycles`[case])
  ngen_cycle <- floor(150/sis_dc$`disease cycles`[case])
  
  # run simulation
  case_dat <- run_gens(case_vect, ngen_cycle, 500)
  
  # save connect to case info
  case_dat$`disease cycles` <- rep(sis_dc$`disease cycles`[case])
  case_dat$`event order` <- rep(sis_dc$`event order`[case])
  case_dat$`transmission` <- rep(sis_dc$`transmission`[case])
  case_dat$`robustness` <- rep(sis_dc$`robustness`[case])
  
  # then append to sim_dat
  dc_dat <- rbind(dc_dat, case_dat)
  
}

# end parallel cluster
stopCluster(cl)

cyc_ord_dat <- dc_dat # G = 0.05

# save things!
saveRDS(cyc_ord_dat, file = "figures/figure_data/cyc_ord_dat.rds")


#-----FIG4: DISEASE FORCE-----
# other variable to check in the number of disease cycles between reproductions
sis_df <- expand.grid("mortality rate" = seq(from = 0.75, to = 2.25, length.out = 7), # probably provides enough insight...
                      "transmission type" = c(1),
                      "transmission rate" = c(1, 2)) # prelim suggests transmission types are similar
df_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

for (case in 1:14){
  
  print(case) # for error id and to check in
  
  case_vect <- base_vect
  case_vect[1] <- 2 # SIS
  
  case_vect[13:15] <- sis_df$`mortality rate`[case]
  ifelse(sis_df$`transmission type`[case] == 1, 
  case_vect[6:8] <- sis_df$`transmission rate`[case], 
  case_vect[6:8] <- 2*sis_df$`transmission rate`[case] # frequency dependence uses higher forcing
  )
  
  # # if mortality
  # case_vect[13] <- 0
  # case_vect[14] <- 0.5
  
  # if recovery
  case_vect[18] <- 1
  case_vect[19] <- 0.5

  # modify base parameters
  case_vect[3] <- sis_df$`transmission type`[case] # since numbers have meaning

  # run simulation
  case_dat <- run_gens(case_vect, 150, 500)
  
  # save connect to case info
  case_dat$`transmission type` <- rep(sis_df$`transmission type`[case])
  case_dat$`base transmission rate` <- rep(sis_df$`transmission rate`[case])
  case_dat$`true transmission rate` <- rep(case_vect[6])
  case_dat$`mortality rate` <- rep(sis_df$`mortality rate`[case])
  
  # then append to sim_dat
  df_dat <- rbind(df_dat, case_dat)
  
}

# end parallel cluster
stopCluster(cl)

dis_force_dat <- df_dat # G = 0.05

# save things!
saveRDS(dis_force_dat, file = "figures/figure_data/dis_force_dat.rds")
