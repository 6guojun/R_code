library("survival")
library("rms")

data("smart")

####perfome nomogram
rt = smart
ddist <- datadist(rt)
time = smart$TEVENT
event = smart$EVENT
y = Surv(time, event)

f2 <- psm(Surv(time, event) ~ DIASTBP + AGE + AAA, data = rt, dist = 'lognormal')
med  <- Quantile(f2)
surv <- Survival(f2)  # This would also work if f was from cph
nom <- nomogram(f2, fun = list(function(x) surv(3, x),
                            function(x) surv(6, x)),
                funlabel=c("3-Month Survival Probability",
                           "6-month Survival Probability"))
plot(nom, xfrac = .45)

#### Predictions form a parametric survival model
# Plots for a parametric survival model
units(t) <- "day"
ddist <- datadist(rt$SEX, rt$DIABETES)
Srv <- Surv(time, event)
f <- cph(Srv ~  AAA + IMT + HDL, data = rt, x = TRUE, y = TRUE, surv=TRUE, u = 365*5)
cal <- calibrate(f, u = 365*5, m = 3, cuts = 3, conf.int = TRUE) # usually B=200 or 300
plot(cal)

