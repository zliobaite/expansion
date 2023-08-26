# 2023 03 17 I.Zliobaite
# patterns of constant expansion present day

data_traits <- read.csv('data_raw/Trait_data_copy.csv', header = TRUE, sep = ",") 

sp_cent <- paste(data_traits[,'Binomial.1.2'],'Centroid',sep='_')
data_traits <- cbind(data_traits,sp_cent)
  
do_mass_plots <- TRUE
do_by_area <- TRUE 
# else by hexabons

orders <- c('Carnivora','Rodentia','Artiodactyla','Hyracoidea','Perissodactyla','Proboscidea','Primates','Notoungulata','Litopterna','Lagomorpha',"Chiroptera","Afrosoricida","Scandentia","Eulipotyphla","Paucituberculata","Didelphimorphia","Peramelemorphia","Dermoptera","Macroscelidea","Monotremata","Dasyuromorphia","Cingulata","Pilosa","Diprotodontia","Pholidota")
# no Microbiotheria
# no Tubulidentata
# no Catecata no Sirenia
# no Notoryctemorphia

perspective <- 'Present'

for (sk in 1:length(orders)){
  order_now <- orders[sk]
  orderG <- paste(order_now,'G',sep='')
  orders <- c(orders,orderG)
}

#hexagon area from one hexagon to km2
hex_area <- (50.3/2)*(50.3/2)*(sqrt(3)/2)*3
#rectangular area from one degree2 to km2
rec_area <- 111*111

results_all <-c()

for (sk in 1:length(orders)){
  order_now <- orders[sk]
  print(order_now)
    
  file_now <- paste('data_raw/sum_occ_',perspective,'_',order_now,'.csv',sep='')
    
  if (file.exists(file_now)){
      
    sum_occ <- read.csv(file_now, header = TRUE, sep = "\t")    
    if (do_by_area){
      sp_range_area <- sum_occ[,'area']*rec_area/1000000  
      data_mass_now <- cbind(sum_occ,sp_range_area)
      sum_occ <- sum_occ[,'area']
      sum_occ <- sqrt(sum_occ*rec_area/1000000)
    }else{
      sp_range_area <- sum_occ[,'hexcount']*hex_area/1000000  
      data_mass_now <- cbind(sum_occ,sp_range_area)
      sum_occ <- sum_occ[,'hexcount']
      sum_occ <- sqrt(sum_occ*hex_area/1000000)
    }
    
    if (length(sum_occ)>0){
      
      file_out_exp <- paste('plots/distributions/fig_range_width_',perspective,'_',order_now,'.pdf',sep='')
      pdf(file_out_exp, width = 5, height = 5)
      hist(sum_occ, freq = FALSE,xlab = 'range th km',breaks = 50, main = paste(order_now, perspective))
      dev.off()
        
      file_out_hex <- paste('plots/distributions/fig_range_vs_hex_',perspective,'_',order_now,'.pdf',sep='')
      pdf(file_out_hex, width = 4.5, height = 5)
      plot(log10(data_mass_now[,'hexcount']),log10(data_mass_now[,'area']))
      dev.off()
      
      range <- sum_occ[order(sum_occ)]
      un_range <- unique(range)
      n_sp <- length(range)
        
      if (n_sp>2){
          
        law_range <- c()
          
        for (sk in 1:length(un_range)){
            
          ind<- which(range>un_range[sk])
          law_range <- rbind(law_range,c(un_range[sk],length(ind),length(ind)/n_sp))
        }
          
        law_range <- cbind(law_range,log10(law_range[,3]))
        law_range <- cbind(law_range,log10(law_range[,2]))
        law_range <- cbind(law_range,log10(law_range[,1]))
          
        colnames(law_range) <- c('range','nexpand','pexpand','log10p','log10n','log10r')
          
        n_points <- dim(law_range)[1]
        #assuming no zero-inflatedness
        law_range <- law_range[1:(n_points-1),]
          
        file_out_range <- paste('data_working/data_law_range_',perspective,'_',order_now,'.csv',sep='')
          
        write.table(law_range, file = file_out_range,col.names = TRUE,row.names = FALSE, sep = '\t')   
        law_range <- read.csv(file_out_range, header = TRUE, sep = "\t")
          
        fit_range_semilog <- lm(log10p ~ range,data = law_range)
        cf <- coef(fit_range_semilog)
        intercept_semilog <- round(cf[1],digits = 3)
        slope_semilog <- round(cf[2],digits = 3)
        r2_semilog <- round(summary(fit_range_semilog)$r.squared,digits = 4)
          
        fit_range_semilog2 <- lm(log10p ~ I(range^2),data = law_range)
        cf <- coef(fit_range_semilog2)
        intercept_semilog2 <- round(cf[1],digits = 3)
        slope_semilog2 <- round(cf[2],digits = 3)
        r2_semilog2 <- round(summary(fit_range_semilog2)$r.squared,digits = 4)
          
        fit_range_plain <- lm(pexpand ~ range,data = law_range)
        cf <- coef(fit_range_plain)
        intercept_plain <- round(cf[1],digits = 3)
        slope_plain <- round(cf[2],digits = 3)
        r2_plain <- round(summary(fit_range_plain)$r.squared,digits = 4)
          
        fit_range_plain2 <- lm(pexpand ~ I(range^2),data = law_range)
        cf <- coef(fit_range_plain2)
        intercept_plain2 <- round(cf[1],digits = 3)
        slope_plain2 <- round(cf[2],digits = 3)
        r2_plain2 <- round(summary(fit_range_plain2)$r.squared,digits = 4)
          
        fit_range_loglog <- lm(log10p ~ log10r,data = law_range)
        cf <- coef(fit_range_loglog)
        intercept_loglog <- round(cf[1],digits = 3)
        slope_loglog <- round(cf[2],digits = 3)
        r2_loglog <- round(summary(fit_range_loglog)$r.squared,digits = 4)
          
          
        file_name_fig_law <- paste('plots/law/fig_law_range_',perspective,'_',order_now,'.pdf',sep='')
        pdf(file_name_fig_law, width = 4.5, height = 5)
        plot(law_range[,'range'],log10(law_range[,'pexpand']),pch=16,cex = 1,xlab = 'Range width, Tkm',ylab = 'log10 (Prop. of species expanding)',main = paste(order_now, perspective))
        abline(fit_range_semilog)
        #text(21.2,-0.2,intercept,cex = 0.8)
        #text(22,-0.5,slope,cex = 0.8)
        dev.off()
          
        file_name_fig_law <- paste('plots/law2/fig_law2_range_',perspective,'_',order_now,'.pdf',sep='')
        pdf(file_name_fig_law, width = 4.5, height = 5)
        plot(law_range[,'range']^2,log10(law_range[,'pexpand']),pch=16,cex = 1,xlab = 'Range area, Mkm2',ylab = 'log10 (Prop. of species expanding)',main = paste(order_now, perspective))
        abline(fit_range_semilog2)
        #text(21.2,-0.2,intercept,cex = 0.8)
        #text(22,-0.5,slope,cex = 0.8)
        dev.off()
          
        if (do_mass_plots){
          if (dim(data_mass_now)[1]>1){
            masses_all <- c()
            for (sk5 in 1:dim(data_mass_now)[1]){
              sp_now <- data_mass_now[sk5,1]
              ind <- which(data_traits[,'sp_cent']==sp_now)
              if (length(ind)==0){
                ind <- which(data_traits[,'Genus.1.2']==sp_now)
              }
              if (length(ind)==1){
                masses_all <- c(masses_all,log10(data_traits[ind,'Mass.g']))  
              }else{
                mm <- mean(log10(data_traits[ind,'Mass.g']))
                masses_all <- c(masses_all,mm)
              }
            }
            
            data_mass_now <- cbind(data_mass_now,masses_all)
              
            file_name_mass <- paste('plots/mass/fig_mass_',perspective,'_',order_now,'.pdf',sep='')
            pdf(file_name_mass, width = 4, height = 4.5)
            plot(data_mass_now[,'masses_all'],log10(data_mass_now[,'sp_range_area']),pch=16,cex = 1,xlab = 'log10 (Mass, g)',ylab = 'log10 (Range area, Mkm2)')
            cc <- round(cor(data_mass_now[,'masses_all'],log10(data_mass_now[,'sp_range_area'])),digits =2)
            legend('bottomright',legend = cc, bty = "n")
            dev.off()
              
            cor1 <- cor(data_mass_now[,3],data_mass_now[,2])
          }
            
        }
          
        seq_r2 <- c(r2_semilog,r2_semilog2,r2_plain,r2_plain2,r2_loglog)
        seq_labels <- c('law','law2','lin','lin2','pow')
          
        max_r2 <- max(seq_r2)
        ii <- which(seq_r2==max_r2)
          
        if (length(ii)>1){
          #print('kar double')
          #print(seq_labels[ii])
          #print(seq_r2)
          bestest <- NA
          bestest_r2 <- NA
        }else{
          bestest <- seq_labels[ii]    
          bestest_r2 <- seq_r2[ii]
        }
          
        if (length(ii)<1){
          print(seq_r2)
        }
          
        if (n_points>=10){
          if (!do_mass_plots){
            cor1 <- NA
          }
          results_all <- rbind(results_all,c(order_now,perspective,n_sp,intercept_semilog,slope_semilog,r2_semilog,intercept_semilog2,slope_semilog2,r2_semilog2,intercept_plain,slope_plain,r2_plain,intercept_plain2,slope_plain2,r2_plain2,intercept_loglog,slope_loglog,r2_loglog,bestest,bestest_r2,cor1))    
        }
          
      }
    }
  }
}

colnames(results_all) <- c('Order','Perspective','Nsp','interceptSemilog','slopeSemilog','r2Semilog','interceptSemilog2','slopeSemilog2','r2Semilog2','interceptPlain','slopePlain','r2Plain','interceptPlain2','slopePlain2','r2Plain2','interceptLoglog','slopeLoglog','r2Loglog','BestModel','R2Best','corr')

write.table(results_all, file = 'outputs/law_range_results.csv',col.names = TRUE,row.names = FALSE, sep = '\t')   

    













