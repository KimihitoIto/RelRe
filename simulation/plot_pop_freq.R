require(dplyr)

df_frequency <- read.csv("pop_freq.csv") %>% select(-c(average_c,average_k))
df_frequency$date <- as.Date(df_frequency$date)
n <- ncol(df_frequency)
dates <- df_frequency$date
color <- rainbow(n) 

pdf("pop_freq.pdf", width=8, height=6)
plot(dates, rep(0,length(dates)), xlab="", xaxt="n", type="n", ylab="Frequency", ylim=c(0,1))
axis(side=1,at=dates,labels=format(dates,"%h %d"),tick=T,las=2)

for (i in 2:n){
    lines(dates, df_frequency[,i], lty=1, lwd=1, col=color[i])
}
legend("topright", legend=colnames(df_frequency)[2:n], col=color[2:n], lty=1, lwd=1, bty="n", cex = 0.7)
dev.off()



