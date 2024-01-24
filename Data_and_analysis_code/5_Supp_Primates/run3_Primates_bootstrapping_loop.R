# 2023 03 16 I.Zliobaite
# 2024 01 20 adding hominin cleaning


data_hom <- read.csv('data_working/data_sum_hominini_man.csv', header = TRUE, sep = "\t", na.string = c('\\N','NA',''))

do_lumping <- TRUE

cik <- 100 # how many bootstrapping iterations

set.seed(1981)

results_duration_all <- c()
results_range_all <- c()

for (cik in 1:cik){
  
  print('iteration')
  print(cik)
  
  data_all <- read.csv('data_raw/now_pub_20230818.csv', header = TRUE, sep = "\t", na.string = '\\N')

  ind_hominini <- which(data_all[,'SUBFAMILY']=='Hominini')
  
  ind <- intersect(ind_hominini,which(data_all[,'SPECIES']=='rhodesiensis'))
  data_all[ind,'SPECIES'] <- 'heidelbergensis'
  
  if (do_lumping){
    for (sk in 1:dim(data_hom)[1]){
      sp_now <- data_hom[sk,'SPECIES']
      ind <- intersect(ind_hominini,which(data_all[,'SPECIES']==sp_now))
      if (is.na(data_hom[sk,'Synsp'])){
        data_all[ind,'SUBFAMILY'] <- 'Hominini0'
        print(sp_now)
      }else{
        data_all[ind,'SPECIES'] <- data_hom[sk,'Synsp']
      }
    }
  }
  
  boo <- sample(1:dim(data_all)[1],replace = TRUE)
  boo <- unique(boo)    
  
  data_all <- data_all[boo,]
  
  age_threshold <- 3 #localities with larger than this time intervals are excluded
  
  by_species <- TRUE
  # FALSE = by genera
  
  ind <- which(data_all[,'ORDER']!='Sirenia')
  data_all <- data_all[ind,]
  ind <- which(data_all[,'ORDER']!='Cetacea')
  data_all <- data_all[ind,]
  ind <- which(data_all[,'ORDER']!='incertae sedis')
  data_all <- data_all[ind,]
  
  data_all[,'LAT'] <- round(data_all[,'LAT'],digits = 1)
  data_all[,'LONG'] <- round(data_all[,'LONG'],digits = 1)
  
  if (by_species){
    genusspecies <- paste( data_all[,'GENUS'], data_all[,'SPECIES']) 
    print('by species')
  }else{
    genusspecies <- data_all[,'GENUS']  
    print('by genus')
  }
  
  midage <- (data_all[,'MAX_AGE'] + data_all[,'MIN_AGE'])/2
  midage <- round(midage,digits = 2)
  diffage <- (data_all[,'MAX_AGE'] - data_all[,'MIN_AGE'])
  data_all <- cbind(data_all,genusspecies,midage,diffage)
  
  print(dim(data_all))
  ind <- which(data_all[,'diffage'] <= age_threshold)
  data_all <- data_all[ind,]
  
  print(dim(data_all))
  indet_set <- c('Gen.', 'gen.','indet.','sp.','Indet','indet','Gen','gen','Sp','sp','incertae sedis')
  ind <- which(!(data_all[,'GENUS'] %in% indet_set))
  data_all <- data_all[ind,]
  if (by_species){
    ind <- which(!(data_all[,'SPECIES'] %in% indet_set))
    data_all <- data_all[ind,]
    print(dim(data_all))
  }
  
  # getting rid of genera or species that have alternative identificiations
  ind <- which(!grepl( '/', data_all[,'GENUS'], fixed = TRUE))
  data_all <- data_all[ind,]
  if (by_species){
    ind <- which(!grepl( '/', data_all[,'SPECIES'], fixed = TRUE))
    data_all <- data_all[ind,]
  }
  print(dim(data_all))
  
  un_sp <- unique(data_all[,'genusspecies'])
  print('No. unique taxa')
  print(length(un_sp))
  
  data_sum <- c()
  
  for (sk in 1:length(un_sp)){
    sp_now <- un_sp[sk]
    
    ind_sp_loc <- which(data_all[,'genusspecies']==sp_now)
    
    no_loc <- length(ind_sp_loc)
    no_countries <- length(unique(data_all[ind_sp_loc,'COUNTRY']))
    min_mid_age <- min(data_all[ind_sp_loc,'midage'])
    max_mid_age <- max(data_all[ind_sp_loc,'midage'])
    sp_duration <- round(max_mid_age - min_mid_age,digits = 2)
    
    min_lat <- min(data_all[ind_sp_loc,'LAT'])
    ind_geo <- which(data_all[ind_sp_loc,'LAT']==min_lat)
    max_lat <- max(data_all[ind_sp_loc,'LAT'])
    ind_geo <- c(ind_geo,which(data_all[ind_sp_loc,'LAT']==max_lat))
    min_long <- min(data_all[ind_sp_loc,'LONG'])
    ind_geo <- c(ind_geo,which(data_all[ind_sp_loc,'LONG']==min_long))
    max_long <- max(data_all[ind_sp_loc,'LONG'])
    ind_geo <- c(ind_geo,which(data_all[ind_sp_loc,'LONG']==max_long))
    age_peak <- mean(data_all[ind_sp_loc[ind_geo],'midage'])
    duration_peak <- age_peak - min_mid_age
    
    diff_lat <- round(max_lat - min_lat,digits = 2)
    diff_long <- round(max_long - min_long,digits = 2)
    sp_area_el <- 3.14*(diff_lat/2)*(diff_long/2)*111*111/1000000
    sp_area_square <- diff_lat*diff_long*111*111/1000000
    #sp_range_width<- sqrt(sp_area_el)
    sp_range_width <- sqrt(sp_area_square)
    
    data_sp <- data_all[ind_sp_loc[1],c('SIDNUM','ORDER','FAMILY','SUBFAMILY','GENUS','SPECIES','genusspecies','BODYMASS','TCRWNHT','CROWNTYP')]
    data_sum <- rbind(data_sum, c(data_sp,sp_now,no_loc,no_countries,sp_duration,max_mid_age,min_mid_age,sp_area_square,sp_area_el,sp_range_width,max_lat,min_lat,max_long,min_long,duration_peak))
  }
  
  colnames(data_sum) <- c('SIDNUM','ORDER','FAMILY','SUBFAMILY','GENUS','SPECIES','genusspecies','BODYMASS','TCRWNHT','CROWNTYP','sp_now','no_loc','no_countries','sp_duration','max_mid_age','min_mid_age','sp_area_square','sp_area_el','sp_range_width','max_lat','min_lat','max_long','min_long','duration_peak')
  
  write.table(data_sum, file = "data_working/data_sum_boo_temp.csv",col.names = TRUE,row.names = FALSE, sep = '\t')       
  
  data_sum <- read.csv('data_working/data_sum_boo_temp.csv', header = TRUE, sep = "\t")  
  
  
  ind <- which(data_sum[,'ORDER']=='Primates')
  data_sum <- data_sum[ind,]
  data_sum1 <- data_sum
  data_sum[,'ORDER'] <- 'Primates other'
  ind <- which(data_sum[,'FAMILY']== 'Hominidae')
  ind <- intersect(ind,which(data_sum[,'SUBFAMILY']!= 'Hominini'))
  data_sum[ind,'ORDER'] <- 'Hominidae not Hominini'
  ind <- which(data_sum[,'SUBFAMILY']=='Hominini')
  if (do_lump){
    data_sum[ind,'ORDER'] <- 'Hominini conservative'  
  }else{
    data_sum[ind,'ORDER'] <- 'Hominini'
  }
  
  data_sum <- rbind(data_sum,data_sum1)
  
  un_orders <- unique(data_sum[,'ORDER'])
  #un_orders <- c('Primates')
  
  do_remove_zeros <- FALSE #FALSE means replace with small values, that variant is better 
  #for the analysis such that the number of taxa stays accurate, but it makes mass plots out of scale, 
  #so for better visible mass plots change to TRUE
  
  min_points <- 9
  
  if (do_remove_zeros){
    print('removing zero durations')
    ind <- which(data_sum[,'sp_range_width']>0)
    data_sum <- data_sum[ind,] 
    ind <- which(data_sum[,'sp_area_square']>0)
    data_sum <- data_sum[ind,] 
    ind <- which(data_sum[,'sp_duration']>0)
    data_sum <- data_sum[ind,] 
    
  }else{
    #replace zeros with very small non-zero values such that logarithms do not fail
    print('replacing zeros with small values')
    ind <- which(data_sum[,'sp_range_width']==0)
    data_sum[ind,'sp_range_width'] <- 0.000000000000000000001
    ind <- which(data_sum[,'sp_area_square']==0)
    data_sum[ind,'sp_area_square'] <- 0.000000000000000000001
    ind <- which(data_sum[,'sp_duration']==0)
    data_sum[ind,'sp_duration'] <- 0.000000000000000000001
  }
  
  results_duration <- c()
  results_range <- c()
  
  for (sk in 1:length(un_orders)){
    
    order_now <- un_orders[sk]
    #print(order_now)
    
    ind <- which(data_sum[,'ORDER']==order_now)  
    data_now <- data_sum[ind,]
   
    
    # LAW DURATION
    
    law_duration <- c()
    dur <- data_now[,'sp_duration']
    # order durations
    dur <- dur[order(dur)]
    # take unique durations as datapoints
    un_dur <- unique(dur)
    n_sp <- length(dur)
    
    # if we have more than min_points unique duration points (we will discard first and last point and will be left with at least 10)
    if (length(un_dur)>=min_points){
      for (sk in 1:length(un_dur)){
        
        # how many survive past threshold
        ind <- which(dur>un_dur[sk])
        
        #statistics for the law
        law_duration <- rbind(law_duration,c(un_dur[sk],length(ind),length(ind)/n_sp))
        
        
      }
      
      law_duration <- cbind(law_duration,log10(law_duration[,3]))
      law_duration <- cbind(law_duration,log10(law_duration[,2]))
      law_duration <- cbind(law_duration,log10(law_duration[,1]))
      
      colnames(law_duration) <- c('duration','nsurvive','psurvive','log10p','log10n','log10r')
      
      n_points <- dim(law_duration)[1]
      
      #print(law_duration)
      
      if (n_points>=min_points){
        # KEY we remove step 1 and last step
        # the first is zero duration, the second is zero survival
        law_duration <- law_duration[2:(n_points-1),]
        
        
        # interim save 
        file_name <- paste('outputs/law_duration_by_orders/points_law_duration_',order_now,'.csv',sep = '')
        write.table(law_duration, file = file_name,col.names = TRUE,row.names = FALSE, sep = '\t')   
        law_duration <- read.csv(file_name, header = TRUE, sep = "\t")
        
        # regression fits
        
        fit_duration_semilog <- lm(log10p ~ duration,data = law_duration)
        cf <- coef(fit_duration_semilog)
        intercept_semilog <- round(cf[1],digits = 3)
        slope_semilog <- round(cf[2],digits = 3)
        r2_semilog <- round(summary(fit_duration_semilog)$r.squared,digits = 3)
        
        fit_duration_plain <- lm(psurvive ~ duration,data = law_duration)
        cf <- coef(fit_duration_plain)
        intercept_plain <- round(cf[1],digits = 3)
        slope_plain <- round(cf[2],digits = 3)
        r2_plain <- round(summary(fit_duration_plain)$r.squared,digits = 3)
        
        fit_duration_loglog <- lm(log10p ~ log10r,data = law_duration)
        cf <- coef(fit_duration_loglog)
        intercept_loglog <- round(cf[1],digits = 3)
        slope_loglog <- round(cf[2],digits = 3)
        r2_loglog <- round(summary(fit_duration_loglog)$r.squared,digits = 3)
        
        # plots by order 
        
        #file_name_fig_law <- paste('plots/law_duration/fig_dur_plain_',order_now,'.pdf',sep='')
        #pdf(file_name_fig_law, width = 3.5, height = 4)
        #plot(law_duration[,'duration'],law_duration[,'psurvive'],pch=16,cex = 1,xlab = 'time, Myr',ylab = 'Proportion of taxa surviving',main = paste(order_now,'plain'))
        #abline(fit_duration_plain)
        #legend('topright',legend = paste('R2 =',round(r2_plain,digits = 2)),bty = "n",cex=0.7)
        #dev.off()  
        
        #file_name_fig_law <- paste('plots/law_duration/fig_dur_semilog_',order_now,'.pdf',sep='')
        #pdf(file_name_fig_law, width = 3.5, height = 4)
        #plot(law_duration[,'duration'],log10(law_duration[,'psurvive']),pch=16,cex = 1,xlab = 'time, Myr',ylab = 'log10 (Proportion of taxa surviving)',main = paste(order_now,'semi-log'))
        #abline(fit_duration_semilog)
        #legend('topright',legend = paste('R2 =',round(r2_semilog,digits = 2)),bty = "n",cex=0.7)
        #dev.off()  
        
        #file_name_fig_law <- paste('plots/law_duration/fig_dur_loglog_',order_now,'.pdf',sep='')
        #pdf(file_name_fig_law, width = 3.5, height = 4)
        #plot(log10(law_duration[,'duration']),log10(law_duration[,'psurvive']),pch=16,cex = 1,xlab = 'log10 (time, Myr)',ylab = 'log10 (Proportion of taxa surviving)',main = paste(order_now,'log-log'))
        #legend('topright',legend = paste('R2 =',round(r2_loglog,digits = 2)),bty = "n",cex=0.7)
        #abline(fit_duration_loglog)
        #dev.off()  
        
        seq_r2 <- c(r2_semilog,r2_plain,r2_loglog)
        seq_labels <- c('exponential','linear','power')
        
        max_r2 <- max(seq_r2)
        ii <- which(seq_r2==max_r2)
        
        if (length(ii)>1){
          print('troblem, long max r2')
          bestest <- paste(seq_labels[ii],collapse = ' ')
        }else{
          bestest <- seq_labels[ii]    
          bestest_r2 <- seq_r2[ii]
        }
        
        if (length(ii)<1){
          print(seq_r2)
        }
        
        results_duration <- rbind(results_duration,c(order_now,n_sp,intercept_semilog,slope_semilog,r2_semilog,intercept_plain,slope_plain,r2_plain,intercept_loglog,slope_loglog,r2_loglog,bestest,bestest_r2)) 
        
      }
      
      
      
    }
    
    # LAW RANGE
    
    law_range <- c()
    
    ran <- data_now[,'sp_range_width']
    ran <- ran[order(ran)]
    un_ran <- unique(ran)
    
    n_sp <- length(ran)
    
    if (length(un_ran)>min_points){
      for (sk in 1:length(un_ran)){
        ind<- which(ran>un_ran[sk])
        
        law_range <- rbind(law_range,c(un_ran[sk],length(ind),length(ind)/n_sp))
      }
      
      law_range <- cbind(law_range,log10(law_range[,3]))
      law_range <- cbind(law_range,log10(law_range[,2]))
      law_range <- cbind(law_range,log10(law_range[,1]))
      
      colnames(law_range) <- c('range','nexpand','pexpand','log10p','log10n','log10r')
      
      n_points <- dim(law_range)[1]
      
      if (n_points>=min_points){
        law_range <- law_range[2:(n_points-1),]
        
        file_name <- paste('outputs/law_range_by_orders/points_law_range_',order_now,'.csv',sep='')
        
        write.table(law_range, file = file_name,col.names = TRUE,row.names = FALSE, sep = '\t')   
        law_range <- read.csv(file_name, header = TRUE, sep = "\t")
        
        fit_range_semilog <- lm(log10p ~ range,data = law_range)
        cf <- coef(fit_range_semilog)
        intercept_semilog <- round(cf[1],digits = 3)
        slope_semilog <- round(cf[2],digits = 3)
        r2_semilog <- round(summary(fit_range_semilog)$r.squared,digits = 3)
        
        fit_range_semilog2 <- lm(log10p ~ I(range^2),data = law_range)
        cf <- coef(fit_range_semilog2)
        intercept_semilog2 <- round(cf[1],digits = 3)
        slope_semilog2 <- round(cf[2],digits = 3)
        r2_semilog2 <- round(summary(fit_range_semilog2)$r.squared,digits = 3)
        
        fit_range_plain <- lm(pexpand ~ range,data = law_range)
        cf <- coef(fit_range_plain)
        intercept_plain <- round(cf[1],digits = 3)
        slope_plain <- round(cf[2],digits = 3)
        r2_plain <- round(summary(fit_range_plain)$r.squared,digits = 3)
        
        fit_range_plain2 <- lm(pexpand ~ I(range^2),data = law_range)
        cf <- coef(fit_range_plain2)
        intercept_plain2 <- round(cf[1],digits = 3)
        slope_plain2 <- round(cf[2],digits = 3)
        r2_plain2 <- round(summary(fit_range_plain2)$r.squared,digits = 3)
        
        fit_range_loglog <- lm(log10p ~ log10r,data = law_range)
        cf <- coef(fit_range_loglog)
        intercept_loglog <- round(cf[1],digits = 3)
        slope_loglog <- round(cf[2],digits = 3)
        r2_loglog <- round(summary(fit_range_loglog)$r.squared,digits = 3)
        
        #file_name_fig_law <- paste('plots/law_range/fig_ran_plain_',order_now,'.pdf',sep='')
        #pdf(file_name_fig_law, width = 3.5, height = 4)
        #plot(law_range[,'range'],law_range[,'pexpand'],pch=16,cex = 1,xlab = 'Range width, Tkm',ylab = 'Proportion of taxa expanding',main = paste(order_now,'plain'))
        #abline(fit_range_plain)
        #legend('topright',legend = paste('R2 =',round(r2_plain,digits = 2)),bty = "n",cex=0.7)
        #dev.off()
        
        #file_name_fig_law <- paste('plots/law_range/fig_ran_plain2_',order_now,'.pdf',sep='')
        #pdf(file_name_fig_law, width = 3.5, height = 4)
        #plot(law_range[,'range']^2,law_range[,'pexpand'],pch=16,cex = 1,xlab = 'Range area, Mkm2',ylab = 'Proportion of taxa expanding',main = paste(order_now,'plain2'))
        #abline(fit_range_plain2)
        #legend('topright',legend = paste('R2 =',round(r2_plain2,digits = 2)),bty = "n",cex=0.7)
        #dev.off()
        
        #file_name_fig_law <- paste('plots/law_range/fig_ran_semilog_',order_now,'.pdf',sep='')
        #pdf(file_name_fig_law, width = 3.5, height = 4)
        #plot(law_range[,'range'],log10(law_range[,'pexpand']),pch=16,cex = 1,xlab = 'Range width, Tkm',ylab = 'log10 (Proportion of taxa expanding)',main = paste(order_now,'semi-log'))
        #abline(fit_range_semilog)
        #legend('topright',legend = paste('R2 =',round(r2_semilog,digits = 2)),bty = "n",cex=0.7)
        #dev.off()
        
        #file_name_fig_law <- paste('plots/law_range/fig_ran_semilog2_',order_now,'.pdf',sep='')
        #pdf(file_name_fig_law, width = 3.5, height = 4)
        #plot(law_range[,'range']^2,log10(law_range[,'pexpand']),pch=16,cex = 1,xlab = 'Range area, Mkm2',ylab = 'log10 (Proportion of taxa expanding)',main = paste(order_now,'semi-log2'))
        #abline(fit_range_semilog2)
        #legend('topright',legend = paste('R2 =',round(r2_semilog2,digits = 2)),bty = "n",cex=0.7)
        #dev.off()
        
        #file_name_fig_law <- paste('plots/law_range/fig_ran_loglog_',order_now,'.pdf',sep='')
        #pdf(file_name_fig_law, width = 3.5, height = 4)
        #plot(log10(law_range[,'range']),log10(law_range[,'pexpand']),pch=16,cex = 1,xlab = 'log10 (Range width, Tkm2)',ylab = 'log10 (Proportion of taxa expanding)',main = paste(order_now,'log-log'))
        #abline(fit_range_loglog)
        #legend('topright',legend = paste('R2 =',round(r2_loglog,digits = 2)),bty = "n",cex=0.7)
        #dev.off()
        
        seq_r2 <- c(r2_semilog,r2_semilog2,r2_plain,r2_plain2,r2_loglog)
        seq_labels <- c('Exponential','Exponential2','Linear','Linear2','Power')
        
        max_r2 <- max(seq_r2)
        ii <- which(seq_r2==max_r2)
        
        if (length(ii)>1){
          print('More than one optimum in range')
          bestest <- paste(seq_labels[ii],collapse = ' ')
        }else{
          bestest <- seq_labels[ii]    
          bestest_r2 <- seq_r2[ii]
        }
        
        if (length(ii)<1){
          print(seq_r2)
        }
        
        results_range <- rbind(results_range,c(order_now,n_sp,intercept_semilog,slope_semilog,r2_semilog,intercept_semilog2,slope_semilog2,r2_semilog2,intercept_plain,slope_plain,r2_plain,intercept_plain2,slope_plain2,r2_plain2,intercept_loglog,slope_loglog,r2_loglog,bestest,bestest_r2))  
        
      }
      
      
    }
    
  }
  
  colnames(results_duration) <- c('Order','Nsp','InterceptSemilog','SlopeSemilog','R2Semilog','IntereptPlain','SlopePlain','R2Plain','InterceptLoglog','SlopeLoglog','R2Loglog','BestModel','R2Best')
  colnames(results_range) <- c('Order','Nsp','InterceptSemilog','SlopeSemilog','R2Semilog','InterceptSemilog2','SlopeSemilog2','R2Semilog2','IntereptPlain','SlopePlain','R2Plain','IntereptPlain2','SlopePlain2','R2Plain2','InterceptLoglog','SlopeLoglog','R2Loglog','BestModel','R2Best')
  
  results_duration_all <- rbind(results_duration_all,results_duration)
  results_range_all <- rbind(results_range_all,results_range)
  
  
}


write.table(results_duration_all, file = 'outputs/law_duration_results_boo_all.csv',col.names = TRUE,row.names = FALSE, sep = '\t')   
write.table(results_range_all, file = 'outputs/law_range_results_boo_all.csv',col.names = TRUE,row.names = FALSE, sep = '\t')   







