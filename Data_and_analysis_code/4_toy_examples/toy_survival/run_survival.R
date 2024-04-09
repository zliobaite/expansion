# 2023 08 11 I.Zliobaite
# toy survival

x <- seq(0,15,0.1)

lambda <- 1
a <- 0.001
c <- 3
k <- 0.0666666

y <- exp(-lambda*x)
y2 <- a*x^(-c)
y3 <- -k*x + 1

ll <- 4
wd <- 3.5
ht <- 4

pdf('fig_models_plain.pdf',width = wd,height = ht)
plot(x,y,type='l',lwd=ll,ylim = c(0,1), main = 'Plain plot')
points(x,y2,type='l',lwd=ll,lty=2)
points(x,y3,type='l',lwd=ll, lty=3)
legend('topright',c('Linear decay (I)','Exponential (II)','Power (III)'),lty = c(3,1,2),lwd=2,bty = "n",cex=0.7)
dev.off()

pdf('fig_models_plain_col.pdf',width = wd,height = ht)
plot(x,y,type='l',lwd=ll+2,ylim = c(0,1), main = 'Plain plot',col='darkred',xlab = 'Age',ylab = 'Fraction surviving')
points(x,y2,type='l',lwd=ll+2,col='darkblue')
points(x,y3,type='l',lwd=ll+2, col='darkgreen')
#legend('topright',c('Linear decay (I)','Exponential (II)','Power (III)'),lty = c(3,1,2),lwd=2,bty = "n",cex=0.7)
dev.off()

pdf('fig_models_semilog.pdf',width = wd,height = ht)
plot(x,log10(y),type='l',lwd=ll, main = 'Semi-log plot')
points(x,log10(y2),type='l',lwd=ll,lty=2)
points(x,log10(y3),type='l',lwd=ll,lty=3)
legend(9,-2,c('(I)','(II)','(III)'),lty = c(3,1,2),lwd=2,bty = "n",cex=0.7)
dev.off()

pdf('fig_models_semilog_col.pdf',width = wd,height = ht)
plot(x,log10(y),type='l',lwd=ll+2, main = 'Semi-log plot',col='darkred',xlab = 'Age',ylab = 'log Fraction surviving')
points(x,log10(y2),type='l',lwd=ll+2,col='darkblue')
points(x,log10(y3),type='l',lwd=ll+2,col='darkgreen')
#legend(9,-2,c('(I)','(II)','(III)'),lty = c(3,1,2),lwd=2,bty = "n",cex=0.7)
dev.off()


pdf('fig_models_loglog.pdf',width = wd,height = ht)
plot(log10(x),log10(y),type='l',lwd=ll, main = "Log-log plot")
points(log10(x),log10(y2),type='l',lwd=ll,lty=2)
points(log10(x),log10(y3),type='l',lwd=ll,lty=3)
legend('bottomleft',c('Linear decay (I)','Exponential (II)','Power (III)'),lty = c(3,1,2),lwd=2,bty = "n",cex=0.7)
dev.off()

pdf('fig_models_loglog_col.pdf',width = wd,height = ht)
plot(log10(x),log10(y),type='l',lwd=ll+2, main = "Log-log plot",col='darkred',xlab = 'log Age',ylab = 'log Fraction surviving')
points(log10(x),log10(y2),type='l',lwd=ll+2,col='darkblue')
points(log10(x),log10(y3),type='l',lwd=ll+2,col='darkgreen')
#legend('bottomleft',c('Linear decay (I)','Exponential (II)','Power (III)'),lty = c(3,1,2),lwd=2,bty = "n",cex=0.7)
dev.off()


pdf('fig_survivorship_individuals.pdf',width = wd,height = ht)
plot(x,log10(y),type='l',lwd=ll, main = 'Mortality of individuals',xaxt='n',yaxt='n',xlab = 'Age', ylab = 'Log Proportion of organisms surviving')
points(x,log10(y2),type='l',lwd=ll)
points(x,log10(y3),type='l',lwd=ll)
text(8.4,-1,'Type I\n (e.g. humans)',bty = "n",cex=0.8)
text(10.4,-2.9,'Type II\n (e.g. birds)',bty = "n",cex=0.8)
text(3.4,-5.6,'Type III\n (e.g. trees)',bty = "n",cex=0.8)
dev.off()

pdf('fig_survivorship_species.pdf',width = wd,height = ht)
plot(x,log10(y),type='l',lwd=ll, main = 'Survivorship over time',xaxt='n',yaxt='n',xlab = 'Age', ylab = 'Log Proportion of organisms or taxa surviving')
     #ylab = 'Log Proportion of organisms or taxa alive')
points(x,log10(y2),type='l',lwd=ll)
points(x,log10(y3),type='l',lwd=ll)
text(8,-1,'Type I \n Ageing',bty = "n",cex=0.8)
text(10.5,-3.5,"Type II\n Law",bty = "n",cex=0.8)
text(2.6,-5.4,'Type III\n Acceleration',bty = "n",cex=0.8)
dev.off()

pdf('fig_expansion_species.pdf',width = wd+0.2,height = ht)
#op <- par(mar=c(5, 5, 4, 2) + 0.1)
plot(x,log10(y),type='l',lwd=ll, main = 'Expansion in space',xaxt='n',yaxt='n',xlab = 'Range size', ylab = 'Log Proportion of taxa expanding')
     #ylab = 'Log Proportion of taxa having a larger range')
#par(op)
points(x,log10(y2),type='l',lwd=ll)
points(x,log10(y3),type='l',lwd=ll)
text(8,-1,'Type I\n Slowdown',bty = "n",cex=0.8)
text(10.9,-3,"Type II\n Constant\n expansion",bty = "n",cex=0.8)
text(4,-5.8,'Type III\n Rich get richer',bty = "n",cex=0.8)
dev.off()


#experiment log-normal
n <- 1000
yy <- rnorm(n, mean=0, sd=1)
xx <- exp(yy)
xx <- sort(xx)
dd <- c()
for (sk in 1:(n-1)){
  dd <- rbind(dd,c(xx[sk],(n-sk)/n))
}
dd <- cbind(dd,log10(sqrt(dd[,2])))
dd <- as.data.frame(dd)

fit_dd <- lm(V3 ~ V1,data = dd)


pdf('fig_lognormal.pdf',width = wd,height = ht)
plot(dd[,1],dd[,3],type='l')
abline(fit_dd)
dev.off()


