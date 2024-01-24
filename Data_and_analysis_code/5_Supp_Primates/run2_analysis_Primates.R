# 2023 03 16 I.Zliobaite

do_lump <- TRUE
do_boo <- FALSE

if (do_lump){
  if (do_boo){
    data_sum <- read.csv('data_working/data_sum_lump_boo.csv', header = TRUE, sep = "\t")  
  }else{
    data_sum <- read.csv('data_working/data_sum_lump.csv', header = TRUE, sep = "\t")    
  }
}else{
  if (do_boo){
    data_sum <- read.csv('data_working/data_sum_split_boo.csv', header = TRUE, sep = "\t") 
  }else{
    data_sum <- read.csv('data_working/data_sum_split.csv', header = TRUE, sep = "\t")  
  }
}


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

#ind <- which(data_sum[,'max_mid_age']<3.5)
#data_sum <- data_sum[ind,]


#ind <- which(data_sum[,'SPECIES']=='sapiens')
#data_sum[ind,'ORDER'] <- 'Primates_us'

#ind <- which(data_sum[,'SPECIES']=='erectus')
#data_sum[ind,'ORDER'] <- 'Primates_er'

#ind <- which(data_sum[,'SPECIES']=='heidelbergensis')
#data_sum[ind,'ORDER'] <- 'Primates_he'

#ind <- which(data_sum[,'SPECIES']=='afarensis')
#data_sum[ind,'ORDER'] <- 'Primates_af'

#ind <- which(data_sum[,'SPECIES']=='neanderthalensis')
#data_sum[ind,'ORDER'] <- 'Primates_ne'

#ind <- which(data_sum[,'SPECIES']=='boisei')
#data_sum[ind,'ORDER'] <- 'Primates_bo'

#ind <- which(data_sum[,'sp_duration']<1)
#data_sum <- data_sum[ind,]

write.table(data_sum, file = 'data_working/data_sum_Pri.csv',col.names = TRUE,row.names = FALSE, sep = '\t')   

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
  print(order_now)
  
  ind <- which(data_sum[,'ORDER']==order_now)  
  data_now <- data_sum[ind,]
  
  if (dim(data_now)[1]>=min_points){
    #scatterplots for the appendix, select remove zeros = TRUE
    
    file_now1 <- paste('plots/scatter/fig_duration_range_',order_now,'.pdf',sep='')
    file_now2 <- paste('plots/scatter/fig_duration_range_log_',order_now,'.pdf',sep='')
    file_now3 <- paste('plots/scatter/fig_mass_duration_',order_now,'.pdf',sep='')
    file_now4 <- paste('plots/scatter/fig_mass_range_',order_now,'.pdf',sep='')
    
    pdf(file_now1, width = 4, height = 4.5)
    plot(data_now[,'sp_duration'],data_now[,'sp_area_square'],pch=16,xlab = 'Duration of taxa, Myr',ylab = 'Range area, Mkm2',main = order_now)
    dev.off()
    
    if (sd(log10(data_now[, "sp_duration"]))>0){
      pdf(file_now2, width = 4, height = 4.5)
      plot(log10(data_now[,'sp_duration']),log10(data_now[,'sp_area_square']),pch=16,xlab = 'log10 (Duration of taxa, Myr)',ylab = 'log10 (Range area, Mkm2)',main = order_now)
      cc <- round(cor(log10(data_now[,'sp_duration']),log10(data_now[,'sp_area_square'])),digits = 2)
      legend('bottomright',legend = cc, cex = 0.8, bty = "n")
      dev.off()  
    }

    if (sum(!is.na(data_now[,'BODYMASS']))>0){
      pdf(file_now3, width = 5, height = 5)
      plot(log10(data_now[,'BODYMASS']),log10(data_now[,'sp_duration']),pch=16,xlab = 'log10 mass kg', ylab = 'log10 duration Myr',main = order_now)
      dev.off()
      
      pdf(file_now4, width = 5, height = 5)
      plot(log10(data_now[,'BODYMASS']),log10(data_now[,'sp_range_width']),pch=16,xlab = 'log10 mass kg', ylab = 'log10 max range width Tkm',main = order_now)
      dev.off()  
    }
    
    file_hist1 <- paste('plots/distributions/fig_hist_duration_',order_now,'.pdf',sep='')
    file_hist2 <- paste('plots/distributions/fig_hist_range_',order_now,'.pdf',sep='')
    
    pdf(file_hist1, width = 4.5, height = 5)
    hist(data_now[,'sp_duration'], freq = FALSE, xlab = 'duration Myr', breaks = 50,main = order_now)
    dev.off()
    
    pdf(file_hist2, width = 4.5, height = 5)
    hist(data_now[,'sp_range_width'], freq = FALSE, xlab = 'range width Tkm', breaks = 50,main = order_now)
    dev.off()
    
    
  }
  
  
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
      
      file_name_fig_law <- paste('plots/law_duration/fig_dur_plain_',order_now,'.pdf',sep='')
      pdf(file_name_fig_law, width = 3.5, height = 4)
      plot(law_duration[,'duration'],law_duration[,'psurvive'],pch=16,cex = 1,xlab = 'time, Myr',ylab = 'Proportion of taxa surviving',main = paste(order_now,'plain'))
      abline(fit_duration_plain)
      legend('topright',legend = paste('R2 =',round(r2_plain,digits = 2)),bty = "n",cex=0.7)
      dev.off()  
      
      file_name_fig_law <- paste('plots/law_duration/fig_dur_semilog_',order_now,'.pdf',sep='')
      pdf(file_name_fig_law, width = 3.5, height = 4)
      plot(law_duration[,'duration'],log10(law_duration[,'psurvive']),pch=16,cex = 1,xlab = 'time, Myr',ylab = 'log10 (Proportion of taxa surviving)',main = paste(order_now,'semi-log'))
      abline(fit_duration_semilog)
      legend('topright',legend = paste('R2 =',round(r2_semilog,digits = 2)),bty = "n",cex=0.7)
      dev.off()  
      
      file_name_fig_law <- paste('plots/law_duration/fig_dur_loglog_',order_now,'.pdf',sep='')
      pdf(file_name_fig_law, width = 3.5, height = 4)
      plot(log10(law_duration[,'duration']),log10(law_duration[,'psurvive']),pch=16,cex = 1,xlab = 'log10 (time, Myr)',ylab = 'log10 (Proportion of taxa surviving)',main = paste(order_now,'log-log'))
      legend('topright',legend = paste('R2 =',round(r2_loglog,digits = 2)),bty = "n",cex=0.7)
      abline(fit_duration_loglog)
      dev.off()  
      
      seq_r2 <- c(r2_semilog,r2_plain,r2_loglog)
      seq_labels <- c('exponential','linear','power')
      
      max_r2 <- max(seq_r2)
      ii <- which(seq_r2==max_r2)
      
      if (length(ii)>1){
        print('troblem, long max r2')
        bestest <- cat(seq_labels[ii])
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
      
      file_name_fig_law <- paste('plots/law_range/fig_ran_plain_',order_now,'.pdf',sep='')
      pdf(file_name_fig_law, width = 3.5, height = 4)
      plot(law_range[,'range'],law_range[,'pexpand'],pch=16,cex = 1,xlab = 'Range width, Tkm',ylab = 'Proportion of taxa expanding',main = paste(order_now,'plain'))
      abline(fit_range_plain)
      legend('topright',legend = paste('R2 =',round(r2_plain,digits = 2)),bty = "n",cex=0.7)
      dev.off()
      
      file_name_fig_law <- paste('plots/law_range/fig_ran_plain2_',order_now,'.pdf',sep='')
      pdf(file_name_fig_law, width = 3.5, height = 4)
      plot(law_range[,'range']^2,law_range[,'pexpand'],pch=16,cex = 1,xlab = 'Range area, Mkm2',ylab = 'Proportion of taxa expanding',main = paste(order_now,'plain2'))
      abline(fit_range_plain2)
      legend('topright',legend = paste('R2 =',round(r2_plain2,digits = 2)),bty = "n",cex=0.7)
      dev.off()
      
      file_name_fig_law <- paste('plots/law_range/fig_ran_semilog_',order_now,'.pdf',sep='')
      pdf(file_name_fig_law, width = 3.5, height = 4)
      plot(law_range[,'range'],log10(law_range[,'pexpand']),pch=16,cex = 1,xlab = 'Range width, Tkm',ylab = 'log10 (Proportion of taxa expanding)',main = paste(order_now,'semi-log'))
      abline(fit_range_semilog)
      legend('topright',legend = paste('R2 =',round(r2_semilog,digits = 2)),bty = "n",cex=0.7)
      dev.off()
      
      file_name_fig_law <- paste('plots/law_range/fig_ran_semilog2_',order_now,'.pdf',sep='')
      pdf(file_name_fig_law, width = 3.5, height = 4)
      plot(law_range[,'range']^2,log10(law_range[,'pexpand']),pch=16,cex = 1,xlab = 'Range area, Mkm2',ylab = 'log10 (Proportion of taxa expanding)',main = paste(order_now,'semi-log2'))
      abline(fit_range_semilog2)
      legend('topright',legend = paste('R2 =',round(r2_semilog2,digits = 2)),bty = "n",cex=0.7)
      dev.off()
      
      file_name_fig_law <- paste('plots/law_range/fig_ran_loglog_',order_now,'.pdf',sep='')
      pdf(file_name_fig_law, width = 3.5, height = 4)
      plot(log10(law_range[,'range']),log10(law_range[,'pexpand']),pch=16,cex = 1,xlab = 'log10 (Range width, Tkm2)',ylab = 'log10 (Proportion of taxa expanding)',main = paste(order_now,'log-log'))
      abline(fit_range_loglog)
      legend('topright',legend = paste('R2 =',round(r2_loglog,digits = 2)),bty = "n",cex=0.7)
      dev.off()
      
      seq_r2 <- c(r2_semilog,r2_semilog2,r2_plain,r2_plain2,r2_loglog)
      seq_labels <- c('Exponential','Exponential2','Linear','Linear2','Power')
      
      max_r2 <- max(seq_r2)
      ii <- which(seq_r2==max_r2)
      
      if (length(ii)>1){
        print('More than one optimum in range')
        bestest <- cat(seq_labels[ii])
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

write.table(results_duration, file = 'outputs/law_duration_results.csv',col.names = TRUE,row.names = FALSE, sep = '\t')   
write.table(results_range, file = 'outputs/law_range_results.csv',col.names = TRUE,row.names = FALSE, sep = '\t')   
