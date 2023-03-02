library(dplyr);

df_count <- read.csv("Tokyo_BA1_BA2.csv");
df_count$date <- as.Date(df_count$date);

df_frequency <- read.csv("frequencies.csv");
df_frequency$date <- as.Date(df_frequency$date);

t_end <- max(df_count$date+6);

f_Omicron_BA1=df_count$Omicron_BA1/rowSums(df_count[c("Omicron_BA1","Omicron_BA2")],na.rm=T);
f_Omicron_BA2=df_count$Omicron_BA2/rowSums(df_count[c("Omicron_BA1","Omicron_BA2")],na.rm=T);

dates <- df_frequency$date;
color <- c("blue", "red")

pdf("frequencies.pdf", width=5, height=5)
plot(dates, rep(0,length(dates)), xlab="", xaxt="n", type="n", ylab="Frequency", ylim=c(0,1))
points(df_count$date+3, f_Omicron_BA1, pch=1, col=color[1]);
points(df_count$date+3, f_Omicron_BA2, pch=2, col=color[2]);

axis(side=1,at=dates,labels=format(dates,"%h %d"),tick=T,las=2);
dates_est <- dates<=t_end;
dates_prd <- dates>=t_end;

lines(dates[dates_est], df_frequency$Omicron_BA1[dates_est], lty=1, lwd=1, col=color[1]);
lines(dates[dates_prd], df_frequency$Omicron_BA1[dates_prd], lty=2, lwd=1, col=color[1]);
lines(dates, df_frequency$Omicron_BA1_lb, lty=3, lwd=1, col=color[1]);
lines(dates, df_frequency$Omicron_BA1_ub, lty=3, lwd=1, col=color[1]);

lines(dates[dates_est], df_frequency$Omicron_BA2[dates_est], lty=1, lwd=1, col=color[2]);
lines(dates[dates_prd], df_frequency$Omicron_BA2[dates_prd], lty=2, lwd=1, col=color[2]);
lines(dates, df_frequency$Omicron_BA2_lb, lty=3, lwd=1, col=color[2]);
lines(dates, df_frequency$Omicron_BA2_ub, lty=3, lwd=1, col=color[2]);

legend("right", legend=c("Omicron_BA1","Omicron_BA2"), lty=1, pch=c(1,2), lwd=2, col=c("blue","red"), bty="n")

dev.off();



