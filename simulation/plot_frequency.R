library(dplyr);
n = 12

df_frequency <- read.csv("frequencies.csv");
df_frequency$date <- as.Date(df_frequency$date);


dates <- df_frequency$date;
color <- rainbow(12) 

pdf("frequencies.pdf", width=8, height=6)
plot(dates, rep(0,length(dates)), xlab="", xaxt="n", type="n", ylab="Frequency", ylim=c(0,1))

axis(side=1,at=dates,labels=format(dates,"%h %d"),tick=T,las=2);
for (i in 1:12){
    lines(dates, df_frequency[,i], lty=1, lwd=1, col=color[i]);
}


#legend("right", legend=c("Delta","Omicron_BA1","Omicron_BA2"), lty=1, pch=c(1,2,4), lwd=2, col=color, bty="n")

dev.off();



