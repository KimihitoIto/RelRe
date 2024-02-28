library(dplyr);

df_count <- read.csv("Tokyo_BA1_BA2.csv");
df_count$date_from <- as.Date(df_count$date_from);
df_count$date_till <- as.Date(df_count$date_till);

df_frequency <- read.csv("frequencies.csv");
df_frequency$date <- as.Date(df_frequency$date);

t_end <- max(df_count$date_till);

dates <- df_frequency$date;

pdf("average_k.pdf", width=5, height=5);

dates_est <- dates<=t_end;
dates_prd <- dates>=t_end;

plot(dates, rep(0,length(dates)), xlab="", xaxt="n", type="n", ylab=expression("Average " * R[t] * " relative to BA.1"), ylim=c(1.0,1.5));
axis(side=1,at=dates,labels=format(dates,"%h %d"),tick=T,las=2);

lines(dates[dates_est], df_frequency$average_k[dates_est], lty=1, lwd=1, col="black");
lines(dates[dates_prd], df_frequency$average_k[dates_prd], lty=2, lwd=1, col="black");

lines(dates,df_frequency$average_k_lb, lty=3, lwd=1, col="black");
lines(dates,df_frequency$average_k_ub, lty=3, lwd=1, col="black");

dev.off();



