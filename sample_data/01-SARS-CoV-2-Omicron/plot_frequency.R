library(dplyr);

df_count <- read.csv("Denmark_Delta_Omicron.csv");
df_count$date_from <- as.Date(df_count$date_from);
df_count$date_till <- as.Date(df_count$date_till);

df_frequency <- read.csv("frequencies.csv");
df_frequency$date <- as.Date(df_frequency$date);

t_end <- max(df_count$date_till);

f_Delta=df_count$Delta/rowSums(df_count[c("Delta","Omicron")],na.rm=T);
f_Omicron=df_count$Omicron/rowSums(df_count[c("Delta","Omicron")],na.rm=T);

dates <- df_frequency$date;
mid_dates <- df_count$date_from+round((df_count$date_till-df_count$date_from)/2)
color <- c("blue", "red")

pdf("frequencies.pdf", width=7, height=5)
plot(dates, rep(0,length(dates)), xlab="", xaxt="n", type="n", ylab="Frequency", ylim=c(0,1))
points(mid_dates, f_Delta, pch=1, col=color[1]);
points(mid_dates, f_Omicron, pch=2, col=color[2]);

axis(side=1,at=dates,labels=format(dates,"%h %d"),tick=T,las=2);
dates_est <- dates<=t_end;
dates_prd <- dates>=t_end;

lines(dates[dates_est], df_frequency$Delta[dates_est], lty=1, lwd=1, col=color[1]);
lines(dates[dates_prd], df_frequency$Delta[dates_prd], lty=2, lwd=1, col=color[1]);
lines(dates, df_frequency$Delta_lb, lty=3, lwd=1, col=color[1]);
lines(dates, df_frequency$Delta_ub, lty=3, lwd=1, col=color[1]);

lines(dates[dates_est], df_frequency$Omicron[dates_est], lty=1, lwd=1, col=color[2]);
lines(dates[dates_prd], df_frequency$Omicron[dates_prd], lty=2, lwd=1, col=color[2]);
lines(dates, df_frequency$Omicron_lb, lty=3, lwd=1, col=color[2]);
lines(dates, df_frequency$Omicron_ub, lty=3, lwd=1, col=color[2]);

legend("right", legend=c("Delta","Omicron"), lty=1, pch=c(1,2), lwd=2, col=c("blue","red"), bty="n")

dev.off();



