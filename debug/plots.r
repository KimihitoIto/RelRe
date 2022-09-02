#Libraries
library(ggplot2)
library(viridis)
library(reshape2)
library(lubridate)
library(dplyr)

#############
#Plot counts#
#############
rreDf <- read.csv("output_6B1A5a/rre.csv")
rreDf <- melt(rreDf, id=c("date"))
colnames(rreDf) <- c("date","clade","count")
rreDf$date <- 
rreDf$date <- as_date(parse_date_time(rreDf$date, orders = c("m/d/y", "m/y", "y", "y-m-d")))

#monthly 
#Plot RE predictions
ggplot(rreDf, aes(x=date, y=count, fill=clade)) +
  geom_area(position = "identity", alpha=0.4) +
  scale_x_date(breaks = "2 years") +
  theme_minimal()

#Plot RE predictions
ggplot(rreDf, aes(x=date, y=count, fill=clade, color=clade)) +
  geom_line(size=2) +
  scale_x_date(breaks = "2 years") +
  scale_y_continuous(breaks = seq(0,10,2)) +
  scale_color_viridis(option="turbo", discrete = TRUE) +
  theme_minimal()
