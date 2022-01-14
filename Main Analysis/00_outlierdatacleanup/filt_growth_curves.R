#!bin/bash
require(tidyverse)
require(growthcurver)

fit_growth_curve <- function(dat, Replicate, Treatment, OD, Time, single_outlier_cutoff = c("sd",3), replicate_outlier_cutoff = c("sd",3), subtract_start=TRUE) {
  #### FOR TROUBLESHOOTING
  # dat = cfs_as_treat_adj
  # Replicate = "Rep"
  # Treatment = "sampleID"
  # OD = "OD_adj"
  # Time = "HoursElapsed"
  # single_outlier_cutoff = c("sd",3)
  # replicate_outlier_cutoff = c("sd",3)
  # subtract_start=TRUE
  
  
  ###### Setting up thresholds #######
  stopifnot(single_outlier_cutoff[1] %in% c("sd","cutoff"))
  stopifnot(replicate_outlier_cutoff[1] %in% c("sd","cutoff"))
  if ( single_outlier_cutoff[1]=="sd" ) {
    single_outlier_func <- function(x,lwrOrupr) {
      # print(x)
        if ( lwrOrupr == "lower"){
          cutoff <- (mean(x) - sd(x)*as.numeric(single_outlier_cutoff[2]))
        } else if (lwrOrupr == "upper") {
          cutoff <- (mean(x) + sd(x)*as.numeric(single_outlier_cutoff[2]))
        }
      return(cutoff)
    }
  } else if (single_outlier_cutoff[1]=="threshold") {
    single_outlier_func <- function(x,lwrOrupr) {
      if ( lwrOrupr == "lower"){
        return(mean(x) - as.numeric(single_outlier_cutoff[2]))
      } else if (lwrOrupr == "upper") {
        return(mean(x) + as.numeric(single_outlier_cutoff[2]))
      } 
    }
  }
  if ( replicate_outlier_cutoff[1]=="sd" ) {
    replicate_outlier_func <- function(x,lwrOrupr) {
      
        if ( lwrOrupr == "lower"){
          cutoff <- (0 - sd(x)*as.numeric(replicate_outlier_cutoff[2]))
        } else if (lwrOrupr == "upper") {
          cutoff <- (0 + sd(x)*as.numeric(replicate_outlier_cutoff[2]))
        } 
      
      return(cutoff)
    }
  } else if (replicate_outlier_cutoff[1]=="threshold") {
    replicate_outlier_func <- function(x,lwrOrupr) {
      if ( lwrOrupr == "lower"){
        return(0 - as.numeric(replicate_outlier_cutoff[2]))
      } else if (lwrOrupr == "upper") {
        return(0 +as.numeric(replicate_outlier_cutoff[2]))
      } 
    }
  }
  
  
  ###### ADJUST DATA; SUBTRACT START #######
  tempDat <- dat %>%  select(paste0(Replicate), paste0(Treatment), paste0(OD), paste0(Time)) %>%
    rename(Replicate = paste0(Replicate), Treatment = paste0(Treatment), OD_raw = paste0(OD), Time = paste0(Time)) 
  if (subtract_start) {
    startingOD <- tempDat %>% filter(Time == 0) %>% mutate(startOD = OD_raw) %>% select(Treatment, Replicate,startOD)
    tempDat <- tempDat %>% left_join(startingOD) %>% mutate(OD = OD_raw-startOD) 
  } else {
    tempDat <- tempDat %>% mutate(OD = OD_raw)
  }
  
  ###### REMOVE SINGLE OUTLIERS #######
  tempDat_adj <- tempDat %>% unite(Replicate, Treatment, col = "RepSamp", remove=FALSE)%>%
    group_by(Treatment, Time) %>% mutate(meanOD_sampletime = median(OD)) %>%
    ungroup() %>% mutate(resOD = OD-meanOD_sampletime) %>% 
    mutate(lowerODres = single_outlier_func(resOD, "lower"), upperODres = single_outlier_func(resOD,"upper")) %>%
    mutate(removePoint = ifelse(resOD > upperODres | resOD < lowerODres, TRUE, FALSE))
  
  
  ###### FIT MODEL BY REP; SINGLE OUTLIERS REMOVED #######
  tempCombo <- unique(tempDat_adj$RepSamp)
  allResults_init <- data.frame()
  for ( x in tempCombo ) {
    # x=tempCombo[1]
    tempfilt <- tempDat_adj %>% filter(RepSamp == x, !removePoint)
    if (nrow(tempfilt)>2) {
      tempfit <- SummarizeGrowth(data_t = tempfilt$Time, data_n = tempfilt$OD)
    } else {
      tempfit <- list()
      tempfit$vals <- c(k=0,r=0,t_gen=0,sigma=0,n0=0)
    }
    allResults_init <- rbind(allResults_init, c(x, unlist(tempfit$vals[c("k","r","t_gen","sigma", "n0")])))
  }
  colnames(allResults_init) <- c("RepSamp", "k","r","t_gen","sigma", "n0")
  # Look at outliers by replicate
  fit_results_withoutliers <- allResults_init %>% separate(RepSamp, into = c("Replicate","Treatment"), remove=FALSE, extra="merge") %>%
    mutate(k = as.numeric(k), r = as.numeric(r), t_gen = as.numeric(t_gen), sigma = as.numeric(sigma), n0 = as.numeric(n0)) %>%
    group_by(Treatment) %>% mutate(k_centred = k-mean(k), r_centred = r-mean(r), t_gen_centred = t_gen-mean(t_gen), sigma_centred = sigma-mean(sigma)) %>%
    ungroup() %>% select(RepSamp, Treatment, k_centred, r_centred, t_gen_centred, sigma_centred) %>%
    pivot_longer(-c(RepSamp, Treatment), names_to = "Metric", values_to = "centred_value") %>%
    group_by(Metric) %>%
    mutate(lowerEst = replicate_outlier_func(centred_value, "lower"), upperEst = replicate_outlier_func(centred_value, "upper")) %>% ungroup() %>%
    mutate(removeRep = ifelse(centred_value<lowerEst | centred_value>upperEst , TRUE, FALSE)) %>%
    mutate(removeRep = ifelse(is.na(removeRep), FALSE, removeRep))
  tempDat_filt <- fit_results_withoutliers %>% filter(Metric!="k_centred") %>%
    group_by(RepSamp, Treatment) %>% summarize(removeRep = any(removeRep)) %>% ungroup() %>% select(RepSamp, removeRep)  %>%
    full_join(tempDat_adj) 
  

  ###### FIT MODEL BY REP; REPLICATE AND SINGLE OUTLIERS REMOVED #######
  tempCombo_filt <- unique(tempDat_filt %>%
                             filter(!removeRep, !removePoint) %>%
                             pull(RepSamp))
  allResults_final <- data.frame()
  for ( x in tempCombo_filt ) {
    tempfilt <- tempDat_filt %>% filter(RepSamp == x, !removePoint)
    if (nrow(tempfilt)>2) {
      tempfit <- SummarizeGrowth(data_t = tempfilt$Time, data_n = tempfilt$OD)
    } else {
      tempfit <- list()
      tempfit$vals <- c(k=0,r=0,t_gen=0,sigma=0, n0=0, t_mid=0)
    }
    allResults_final <- rbind(allResults_final, c(x, unlist(tempfit$vals[c("k","r","t_gen","sigma", "n0", "t_mid")])))
  }
  colnames(allResults_final) <- c("RepSamp", "k","r","t_gen","sigma", "n0", "t_mid")

  # Summarize growth curves
  final_fits <- allResults_final %>% separate(RepSamp, into = c("Replicate","Treatment"), remove=FALSE, extra="merge") %>%
    mutate(k = as.numeric(k), r = as.numeric(r), t_gen = as.numeric(t_gen), sigma = as.numeric(sigma), n0 = as.numeric(n0), t_mid = as.numeric(t_mid))


  ###### FIT FINAL MODEL; POOLED REPLICATES #######
  tempCombo_filt_bytype <- unique(tempDat_filt %>% 
                                    pull(Treatment))
  allResults_final_bytype <- data.frame()
  for ( x in tempCombo_filt_bytype ) {
    tempfilt <- tempDat_filt %>% filter(Treatment == x, !removePoint, !removeRep)
    tempfit <- SummarizeGrowth(data_t = tempfilt$Time, data_n = tempfilt$OD)
    allResults_final_bytype <- rbind(allResults_final_bytype, c(x, unlist(tempfit$vals[c("k","k_se","r","r_se","t_gen","sigma", "n0", "n0_se", "t_mid")])))
  }
  colnames(allResults_final_bytype) <- c("Treatment", "k","k_se","r","r_se","t_gen","sigma", "n0", "n0_se", "t_mid")
  allResults_final_bytype <- allResults_final_bytype %>%
    mutate(k = as.numeric(k),k_se = as.numeric(k_se), r = as.numeric(r),r_se = as.numeric(r_se), t_gen = as.numeric(t_gen), sigma = as.numeric(sigma), n0 = as.numeric(n0), n0_se = as.numeric(n0_se), t_mid = as.numeric(t_mid)) 
  
  ###### PLOTTING #######
  # Single point outliers, standardized by replicate sample
  gg_singleoutliers <- tempDat_adj %>%
    ggplot(aes(x=resOD, fill=Treatment)) + geom_histogram(bins=50) +geom_vline(aes(xintercept = lowerODres), col="red", lty=3) +
    geom_vline(aes(xintercept = upperODres), col="red", lty=3) 
  # Replicate outliers, standardized by sample type
  gg_repoutliers <- fit_results_withoutliers %>%
    ggplot(aes(x=centred_value, fill=Treatment)) + geom_histogram(bins=50) +geom_vline(aes(xintercept = upperEst), col="red", lty=3) +
    geom_vline(aes(xintercept = lowerEst), col="red", lty=3) +facet_wrap(.~Metric, scales="free")
  # Replicate outleirs, unstandardized (raw values)
  gg_rawoutliers <- final_fits %>% pivot_longer(-c(RepSamp, Replicate, Treatment), names_to = "Metric", values_to = "Estimate") %>%
    ggplot(aes(x=Treatment, y = Estimate)) + geom_point() + 
    facet_wrap(.~Metric, scales = "free") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  # Plot total fit curves
  gg_totalfit <- tempDat_filt %>% left_join(allResults_final_bytype) %>% group_by(Treatment) %>%
    mutate(GrowthCurve = k/(1 + ((k-n0)/n0)*exp(-r*Time))-n0) %>% 
    mutate(GrowthCurve = ifelse(is.na(GrowthCurve), 0, GrowthCurve)) %>%
    ggplot(aes(x=Time)) + geom_line(aes(y=OD, group=Replicate, col=removeRep)) + 
    # ))+
    geom_point(aes(y=OD, col=removePoint|removeRep)) +
    geom_line(aes(y=GrowthCurve), col="black") +
    facet_wrap(paste0("t_gen: ",round(t_gen,2))~Treatment)+ labs(col = "Outliers removed\nfrom fitting")
  
  ###### RETURN #######
  return(list("Pooled fit" = allResults_final_bytype # Pooled fits (reps are pooled, then fit)
              , "Fit of individual replicates" = final_fits # Individual fits (each rep gets fit separately)
              , "Plot fit of pooled model" = gg_totalfit
              , "Plot outliers of centred data" = gg_singleoutliers
              , "Plot outliers of centred fit variables" = gg_repoutliers
              , "Plot outliers of raw fit variables" =gg_rawoutliers))
}

