# Copyright 2021 Chayada Piantham and Kimihito Ito

library(dplyr);

df_count <- read.csv("Japan_Linenage_Counts.csv");
df_count$date_from <- as.Date(df_count$date_from);
df_count$date_till <- as.Date(df_count$date_till);

df_frequency <- read.csv("frequencies.csv");
df_frequency$date <- as.Date(df_frequency$date);

t_end <- max(df_count$date_till);

first_R.1 <- df_count$date_from[min(which(df_count$R.1 > 0))];
first_Alpha <- df_count$date_from[min(which(df_count$Alpha > 0))];
first_Delta <- df_count$date_from[min(which(df_count$Delta > 0))];

t_start <- min(first_R.1,first_Alpha,first_Delta);
t_end <- max(df_count$date_till);

f_R.1=df_count$R.1/rowSums(df_count[3:6],na.rm=T);
f_Alpha=df_count$Alpha/rowSums(df_count[3:6],na.rm=T);
f_Delta=df_count$Delta/rowSums(df_count[3:6],na.rm=T);
f_other=df_count$other/rowSums(df_count[3:6],na.rm=T);

dates <- df_frequency$date;
color <- hcl.colors(4);

obs_dates <- df_count$date_from + (df_count$date_till-df_count$date_from)/2

pdf("frequencies.pdf", width=7, height=5)
plot(dates, rep(0,length(dates)), xlab="", xaxt="n", type="n", ylab="Frequency", ylim=c(0,1))
points(obs_dates, f_R.1, pch=1, col=color[1]);
points(obs_dates, f_Alpha, pch=2, col=color[2]);
points(obs_dates, f_Delta, pch=3, col=color[3]);
points(obs_dates, f_other, pch=4, col=color[4]);

axis(side=1,at=dates,labels=format(dates,"%h %d"),tick=T,las=2);
dates_est <- dates<=t_end;
dates_prd <- dates>=t_end;

lines(dates[dates_est], df_frequency$R.1[dates_est], lty=1, lwd=1, col=color[1]);
lines(dates[dates_prd], df_frequency$R.1[dates_prd], lty=2, lwd=1, col=color[1]);
lines(dates, df_frequency$R.1_lb, lty=3, lwd=1, col=color[1]);
lines(dates, df_frequency$R.1_ub, lty=3, lwd=1, col=color[1]);

lines(dates[dates_est], df_frequency$Alpha[dates_est] , lty=1, lwd=1, col=color[2]);
lines(dates[dates_prd], df_frequency$Alpha[dates_prd], lty=2, lwd=1, col=color[2]);
lines(dates, df_frequency$Alpha_lb, lty=3, lwd=1, col=color[2]);
lines(dates, df_frequency$Alpha_ub, lty=3, lwd=1, col=color[2]);

lines(dates[dates_est], df_frequency$Delta[dates_est], lty=1, lwd=1, col=color[3]);
lines(dates[dates_prd], df_frequency$Delta[dates_prd], lty=2, lwd=1, col=color[3]);
lines(dates, df_frequency$Delta_lb, lty=3, lwd=1, col=color[3]);
lines(dates, df_frequency$Delta_ub, lty=3, lwd=1, col=color[3]);

lines(dates[dates_est], df_frequency$other[dates_est], lty=1, lwd=1, col=color[4]);
lines(dates[dates_prd], df_frequency$other[dates_prd], lty=2, lwd=1, col=color[4]);
lines(dates, df_frequency$other_lb, lty=3, lwd=1, col=color[4]);
lines(dates, df_frequency$other_ub, lty=3, lwd=1, col=color[4]);

legend("topleft", legend=c("R1", "Alpha", "Delta", "other"), lty=1, pch=c(1,2,3,4), lwd=2, col=color, bty="n")

dev.off();



