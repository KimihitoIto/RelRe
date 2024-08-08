library(dplyr);

df_count <- read.csv("Denmark_all_lineages.csv");
df_count$date_from <- as.Date(df_count$date_from);
df_count$date_till <- as.Date(df_count$date_till);

df_frequency <- read.csv("frequencies.csv");
df_frequency$date <- as.Date(df_frequency$date);

t_end <- max(df_count$date_till);

f_Delta=df_count$Delta/rowSums(df_count[c("Delta","Omicron_BA1","Omicron_BA2")],na.rm=T);
f_Omicron_BA1=df_count$Omicron_BA1/rowSums(df_count[c("Delta","Omicron_BA1","Omicron_BA2")],na.rm=T);
f_Omicron_BA2=df_count$Omicron_BA2/rowSums(df_count[c("Delta","Omicron_BA1","Omicron_BA2")],na.rm=T);

lb_Omicron_BA1 <- numeric();
ub_Omicron_BA1 <- numeric();
lb_Omicron_BA2 <- numeric();
ub_Omicron_BA2 <- numeric();
lb_Delta <- numeric();
ub_Delta <- numeric();

for(i in 1:nrow(df_count)){
  n_BA1 <- df_count$Omicron_BA1[i];
  n_BA2 <- df_count$Omicron_BA2[i];
  n_Delta <- df_count$Delta[i];
  n_tot <- n_BA1 + n_BA2 + n_Delta;
  lb_Omicron_BA1 <- c(lb_Omicron_BA1, binom.test(n_BA1,n_tot)$conf.int[1]);
  ub_Omicron_BA1 <- c(ub_Omicron_BA1, binom.test(n_BA1,n_tot)$conf.int[2]);
  lb_Omicron_BA2 <- c(lb_Omicron_BA2, binom.test(n_BA2,n_tot)$conf.int[1]);
  ub_Omicron_BA2 <- c(ub_Omicron_BA2, binom.test(n_BA2,n_tot)$conf.int[2]);
  lb_Delta <- c(lb_Delta, binom.test(n_Delta,n_tot)$conf.int[1]);
  ub_Delta <- c(ub_Delta, binom.test(n_Delta,n_tot)$conf.int[2]);    
 }

dates <- df_frequency$date;
color <- c("blue", "orange", "red")

#pdf("frequencies_bar.pdf", width=8, height=6)
svg("frequencies_bar.svg", width=8, height=6)
plot(dates, rep(0,length(dates)), xlab="", xaxt="n", type="n", ylab="Frequency", ylim=c(0,1))
points(df_count$date_from, f_Delta, pch=1, col=color[1]);
points(df_count$date_from, f_Omicron_BA1, pch=2, col=color[2]);
points(df_count$date_from, f_Omicron_BA2, pch=4, col=color[3]);

arrows(df_count$date_from, lb_Delta, df_count$date_from, ub_Delta, col=color[1],length=0.02, angle=90, code=3);
arrows(df_count$date_from, lb_Omicron_BA1, df_count$date_from, ub_Omicron_BA1, col=color[2],length=0.02, angle=90, code=3);
arrows(df_count$date_from, lb_Omicron_BA2, df_count$date_from, ub_Omicron_BA2, col=color[3],length=0.02, angle=90, code=3);

axis(side=1,at=dates,labels=format(dates,"%h %d"),tick=T,las=2);
dates_est <- dates<=t_end;
dates_prd <- dates>=t_end;

lines(dates[dates_est], df_frequency$Delta[dates_est], lty=1, lwd=1, col=color[1]);
lines(dates[dates_prd], df_frequency$Delta[dates_prd], lty=2, lwd=1, col=color[1]);
lines(dates, df_frequency$Delta_lb, lty=3, lwd=1, col=color[1]);
lines(dates, df_frequency$Delta_ub, lty=3, lwd=1, col=color[1]);

lines(dates[dates_est], df_frequency$Omicron_BA1[dates_est], lty=1, lwd=1, col=color[2]);
lines(dates[dates_prd], df_frequency$Omicron_BA1[dates_prd], lty=2, lwd=1, col=color[2]);
lines(dates, df_frequency$Omicron_BA1_lb, lty=3, lwd=1, col=color[2]);
lines(dates, df_frequency$Omicron_BA1_ub, lty=3, lwd=1, col=color[2]);

lines(dates[dates_est], df_frequency$Omicron_BA2[dates_est], lty=1, lwd=1, col=color[3]);
lines(dates[dates_prd], df_frequency$Omicron_BA2[dates_prd], lty=2, lwd=1, col=color[3]);
lines(dates, df_frequency$Omicron_BA2_lb, lty=3, lwd=1, col=color[3]);
lines(dates, df_frequency$Omicron_BA2_ub, lty=3, lwd=1, col=color[3]);

legend("right", legend=c("Delta","Omicron_BA1","Omicron_BA2"), lty=1, pch=c(1,2,4), lwd=2, col=color, bty="n")

dev.off();



