# 2024 01 15 I.Zliobaite

p_decline <- 0.05

tsteps <- 1200

taxa_alive <- matrix(0,tsteps+1,1200+1)

taxa_alive[1,1] <- 1

set.seed(1981) 

do_snapshot <- FALSE

for (sk in 1:tsteps){
  ind_rising <- which(taxa_alive[,sk]>0)
  n_rising <- length(ind_rising)
  probs <- runif(n_rising)
  ind <- which(probs<=p_decline)
  taxa_alive[ind_rising,sk+1] <- taxa_alive[ind_rising,sk] + 1
  taxa_alive[ind_rising[ind],sk+1] <- (-taxa_alive[ind_rising[ind],sk])
  
  ind_declining <- which(taxa_alive[,sk]<0)
  sum_declining <- apply(taxa_alive,1,sum)
  ind_more_decline <- which(sum_declining[ind_declining]>0)
  taxa_alive[ind_declining[ind_more_decline],sk+1] <- taxa_alive[ind_declining[ind_more_decline],sk] +1 
  
  taxa_alive[sk+1,sk+1] <- 1
}

sm <- apply(abs(taxa_alive),2,sum)

pdf('fig_apply.pdf')
plot(sm)
dev.off()

taxa_alive <- taxa_alive[201:(dim(taxa_alive)[1]-200),201:(dim(taxa_alive)[2]-200)]

durations <- c()

for (sk in 1:dim(taxa_alive)[1]){
  ii <- which(taxa_alive[sk,]!=0)
  durations <- c(durations,length(ii))
}

range_wd_by_area <- c()
range_wd_by_steps <- c()
range_wd_by_areasqrt <- c()
range_wd_by_steps2 <- c()

taxa_alive_by_areasqrt <- sqrt(abs(taxa_alive))
taxa_alive_by_area <- abs(taxa_alive)
taxa_alive_by_steps <- abs(taxa_alive)
taxa_alive_by_steps2 <- abs(taxa_alive)*abs(taxa_alive)

if (do_snapshot){
    snapshot_now <- 456 #why not
    range_wd_by_area <- taxa_alive_by_area[,snapshot_now]
    range_wd_by_areasqrt <- taxa_alive_by_areasqrt[,snapshot_now]
    range_wd_by_steps <- taxa_alive_by_steps[,snapshot_now]
    range_wd_by_steps2 <- taxa_alive_by_steps2[,snapshot_now]
    
}else{
  for (sk in 1:dim(taxa_alive)[1]){
    
    range_wd_by_area <- c(range_wd_by_area,max(taxa_alive_by_area[sk,]))
    range_wd_by_areasqrt <- c(range_wd_by_areasqrt,max(taxa_alive_by_areasqrt[sk,]))
    range_wd_by_steps <- c(range_wd_by_steps,max(taxa_alive_by_steps[sk,]))
    range_wd_by_steps2 <- c(range_wd_by_steps2,max(taxa_alive_by_steps2[sk,]))
  }
  
}




un_dur <- unique(durations)
un_dur <- un_dur[order(un_dur)]
law_duration <- c()
n_sp <- sum(durations>0)

for (sk in 1:length(un_dur)){
  
  # how many survive past threshold
  ind <- which(durations>un_dur[sk])
  
  #statistics for the law
  law_duration <- rbind(law_duration,c(un_dur[sk],length(ind),length(ind)/n_sp))
}
  
colnames(law_duration) <- c('duration','nsuv','psurvive')

pdf('fig_duration.pdf', width = 3.5, height = 4)
plot(law_duration[,'duration'],log10(law_duration[,'psurvive']),pch=16,cex = 1,xlab = 'time',ylab = 'log10 (Proportion of taxa surviving)',main = 'Duration')
dev.off()  
  

un_steps <- unique(range_wd_by_steps)
un_steps <- un_steps[order(un_steps)]
law_steps <- c()

for (sk in 1:length(un_steps)){
  
  ind <- which(range_wd_by_steps>un_steps[sk])
  law_steps <- rbind(law_steps,c(un_steps[sk],length(ind),length(ind)/n_sp))
}

law_steps <- cbind(law_steps,log10(law_steps[,3]),log10(law_steps[,1]))

colnames(law_steps) <- c('steps','nsuv','pexpand','log10p','log10r')

pdf('fig_steps.pdf', width = 3.5, height = 4)
plot(law_steps[,'steps'],log10(law_steps[,'pexpand']),pch=16,cex = 1,xlab = 'range width',ylab = 'log10 (Proportion of taxa expanding)',main = 'Expansion by steps')
dev.off()  

  
un_steps2 <- unique(range_wd_by_steps2)
un_steps2 <- un_steps2[order(un_steps2)]
law_steps2 <- c()

for (sk in 1:length(un_steps2)){
  
  ind <- which(range_wd_by_steps2>un_steps2[sk])
  law_steps2 <- rbind(law_steps2,c(un_steps2[sk],length(ind),length(ind)/n_sp))
}

law_steps2 <- cbind(law_steps2,log10(law_steps2[,3]),log10(law_steps2[,1]))

colnames(law_steps2) <- c('steps2','nsuv','pexpand','log10p','log10r')


pdf('fig_steps2.pdf', width = 3.5, height = 4)
plot(law_steps2[,'steps2'],log10(law_steps2[,'pexpand']),pch=16,cex = 1,xlab = 'range area',ylab = 'log10 (Proportion of taxa expanding)',main = 'Expansion by steps')
dev.off()  



un_area <- unique(range_wd_by_area)
un_area <- un_area[order(un_area)]
law_area <- c()

for (sk in 1:length(un_area)){
  
  ind <- which(range_wd_by_area>un_area[sk])
  law_area <- rbind(law_area,c(un_area[sk],length(ind),length(ind)/n_sp))
}

law_area <- cbind(law_area,log10(law_area[,3]),log10(law_area[,1]))

colnames(law_area) <- c('area','nsuv','pexpand','log10p','log10r')


pdf('fig_area.pdf', width = 3.5, height = 4)
plot(law_area[,'area'],log10(law_area[,'pexpand']),pch=16,cex = 1,xlab = 'range area',ylab = 'log10 (Proportion of taxa expanding)',main = 'Expansion by area')
dev.off()  



un_areasqrt <- unique(range_wd_by_areasqrt)
un_areasqrt <- un_areasqrt[order(un_areasqrt)]
law_areasqrt <- c()

for (sk in 1:length(un_areasqrt)){
  
  ind <- which(range_wd_by_areasqrt>un_areasqrt[sk])
  law_areasqrt <- rbind(law_areasqrt,c(un_areasqrt[sk],length(ind),length(ind)/n_sp))
}

law_areasqrt <- cbind(law_areasqrt,log10(law_areasqrt[,3]),log10(law_areasqrt[,1]))

colnames(law_areasqrt) <- c('areasqrt','nsuv','pexpand','log10p','log10r')


pdf('fig_areasqrt.pdf', width = 3.5, height = 4)
plot(law_areasqrt[,'areasqrt'],log10(law_areasqrt[,'pexpand']),pch=16,cex = 1,xlab = 'range width',ylab = 'log10 (Proportion of taxa expanding)',main = 'Expansion by area')
dev.off()  


results <- c()
# tests
law_steps <- as.data.frame(law_steps)
law_steps <- law_steps[2:(dim(law_steps)[1]-1),]
fit <- lm(log10p ~ steps,data = law_steps)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by steps','range width','semilog',r2)
results <- rbind(results,results_now)

fit <- lm(pexpand ~ steps,data = law_steps)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by steps','range width','plain',r2)
results <- rbind(results,results_now)

fit <- lm(log10p ~ log10r,data = law_steps)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by steps','range width','loglog',r2)
results <- rbind(results,results_now)



law_steps2 <- as.data.frame(law_steps2)
law_steps2 <- law_steps2[2:(dim(law_steps2)[1]-1),]
fit <- lm(log10p ~ steps2,data = law_steps2)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by steps','range area','semilog',r2)
results <- rbind(results,results_now)

fit <- lm(pexpand ~ steps2,data = law_steps2)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by steps','range area','plain',r2)
results <- rbind(results,results_now)

fit <- lm(log10p ~ log10r,data = law_steps2)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by steps','range area','loglog',r2)
results <- rbind(results,results_now)




law_areasqrt <- as.data.frame(law_areasqrt)
law_areasqrt <- law_areasqrt[2:(dim(law_areasqrt)[1]-1),]
fit <- lm(log10p ~ areasqrt,data = law_areasqrt)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by area','range width','semilog',r2)
results <- rbind(results,results_now)

fit <- lm(pexpand ~ areasqrt,data = law_areasqrt)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by area','range width','plain',r2)
results <- rbind(results,results_now)

fit <- lm(log10p ~ log10r,data = law_areasqrt)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by area','range width','loglog',r2)
results <- rbind(results,results_now)




law_area <- as.data.frame(law_area)
law_area <- law_area[2:(dim(law_area)[1]-1),]
fit <- lm(log10p ~ area,data = law_area)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by area','range area','semilog',r2)
results <- rbind(results,results_now)

fit <- lm(pexpand ~ area,data = law_area)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by area','range area','plain',r2)
results <- rbind(results,results_now)

fit <- lm(log10p ~ log10r,data = law_area)
r2 <- round(summary(fit)$r.squared,digits = 3)
results_now <- c('Expansion by area','range area','loglog',r2)
results <- rbind(results,results_now)

write.table(results, file = 'results_snapshot.csv',col.names = TRUE,row.names = FALSE, sep = '\t')   
