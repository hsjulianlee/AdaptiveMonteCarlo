
load(url("http://www.stats.gla.ac.uk/~claire/aptsassess.RData"))

library(ggplot2) 
library(quantreg) 

#%% Figure 1 #######

ls_fit <- lm(Wmax ~ Year, MaxSpeed) 
title <- sprintf('LS: speed = %.3f%+.3f*year', ls_fit$coefficients[[1]], ls_fit$coefficients[[2]]) 
ggplot(MaxSpeed, aes(Year, Wmax)) + ggtitle(title) + geom_point() + geom_smooth(method = "lm") 

#%% Figure 2 #######

tau <- c( 0.5, 0.75, 0.95 ) 
q_fit <- rq(Wmax ~ Year, tau, MaxSpeed) 
title <- ""
for (i in 1:ncol(q_fit$coefficients)) {
  if (i > 1) title <- paste0(title, "\n") 
  title <- paste0(title, sprintf('Q(%.2f): speed = %.3f%+.3f*year', tau[i], q_fit$coefficients[1,i], q_fit$coefficients[2,i])) 
}
ggplot(MaxSpeed, aes(Year, Wmax)) + ggtitle(title) + geom_point() + geom_quantile(aes(Year, Wmax), quantiles = tau) 

#%% Figure 3 #######

plot(summary(q_fit), parm="Year", main="Slope coefficients for quantiles") 

#%% Figure 4 #######

tau <- 1:19/20 
q_fit_2 <- rq(Wmax ~ Year, tau, MaxSpeed) 
title <- ""
for (i in 1:ncol(q_fit_2$coefficients)) {
  if (i %% 3 == 1) title <- paste0(title, "\n") else title <- paste(title, "") 
  title <- paste0(title, sprintf('Q(%.2f): slope = %+.3f', tau[i], q_fit_2$coefficients[2,i])) 
}
ggplot(MaxSpeed, aes(Year, Wmax)) + ggtitle(title) + geom_point() + geom_quantile(aes(Year, Wmax), quantiles = tau) 

#%% Figure 5 #######

plot(summary(q_fit_2), parm="Year", main="Slope coefficients for quantiles") 

#%% Figure 6 #######

MaxSpeed_3 <- MaxSpeed 
for (i in 1:nrow(MaxSpeed_3)) {
  MaxSpeed_3[i, 2] <- MaxSpeed_3[i, 1] * ls_fit$coefficients[[2]] + ls_fit$coefficients[[1]] + rnorm(1) 
}
ls_fit_3 <- lm(Wmax ~ Year, MaxSpeed_3) 
title <- sprintf('LS: speed = %.3f%+.3f*year', ls_fit_3$coefficients[[1]], ls_fit_3$coefficients[[2]]) 
ggplot(MaxSpeed_3, aes(Year, Wmax)) + ggtitle(title) + geom_point() + geom_smooth(method = "lm") 

#%% Figure 7 #######

tau <- 1:19/20 
q_fit_3 <- rq(Wmax ~ Year, tau, MaxSpeed_3) 
title <- ""
for (i in 1:ncol(q_fit_3$coefficients)) {
  if (i %% 3 == 1) title <- paste0(title, "\n") else title <- paste(title, "") 
  title <- paste0(title, sprintf('Q(%.2f): slope = %+.3f', tau[i], q_fit_3$coefficients[2,i])) 
}
ggplot(MaxSpeed_3, aes(Year, Wmax)) + ggtitle(title) + geom_point() + geom_quantile(aes(Year, Wmax), quantiles = tau) 

#%% Figure 8 #######

plot(summary(q_fit_3), parm="Year", main="Slope coefficients for quantiles") 

#%% Figure 9 #######

MaxSpeed_4 <- MaxSpeed 
MaxSpeed_4[20, 2] <- 200 
ls_fit_4 <- lm(Wmax ~ Year, MaxSpeed_4) 
title <- sprintf('LS: speed = %.3f%+.3f*year', ls_fit_4$coefficients[[1]], ls_fit_4$coefficients[[2]]) 
ggplot(MaxSpeed_4, aes(Year, Wmax)) + ggtitle(title) + geom_point() + geom_smooth(method = "lm") 

#%% Figure 10 #######

tau <- 1:19/20 
q_fit_4 <- rq(Wmax ~ Year, tau, MaxSpeed_4) 
title <- ""
for (i in 1:ncol(q_fit_4$coefficients)) {
  if (i %% 3 == 1) title <- paste0(title, "\n") else title <- paste(title, "") 
  title <- paste0(title, sprintf('Q(%.2f): slope = %+.3f', tau[i], q_fit_4$coefficients[2,i])) 
}
ggplot(MaxSpeed_4, aes(Year, Wmax)) + ggtitle(title) + geom_point() + geom_quantile(aes(Year, Wmax), quantiles = tau) 

#%% Figure 11 #######

plot(summary(q_fit_4), parm="Year", main="Slope coefficients for quantiles") 
