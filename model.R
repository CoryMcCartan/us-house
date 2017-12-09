# HOUSE MODEL

library(loo)
library(shinystan)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ELECTION.DAY = as.Date("2018-11-06")


polls = read.csv("data/current_polls.csv", colClasses=c(date="Date"))
# drop missing data
polls = polls[complete.cases(polls),]
# week IDs
start.day = min(polls$date)
end.day = max(polls$date)
n.weeks = ceiling(as.numeric(end.day - start.day) / 7)
polls$week = round(as.numeric(polls$date - start.day) / 7) + 1
# pollster IDs
polls$firm.id = as.numeric(as.factor(polls$pollster))
# round
polls$n_dem = round(polls$n_dem)
polls$n_side = round(polls$n_side)
polls$n_resp = round(polls$n_resp)

model.data = list(W = n.weeks,
                  P = max(polls$firm.id),
                  N = nrow(polls),
                  w = polls$week,
                  p = polls$firm.id,
                  n_resp = polls$n_resp,
                  n_side = polls$n_side,
                  n_dem = polls$n_dem)

intent.model = stan(file="intent-model.stan", model_name="intent",
                    data=model.data, iter=10000, warmup=3000, chains=1,
                    control=list(adapt_delta = 0.99))
log.lik.int = extract_log_lik(intent.model, "log_lik")

und.model = stan(file="intent-model-und.stan", model_name="undecided",
                    data=model.data, iter=10000, warmup=3000, chains=1,
                    control=list(adapt_delta = 0.99))
log.lik.und = extract_log_lik(und.model, "log_lik")

launch_shinystan(und.model)

print(intent.model, pars=c("poll_error", "logit_poll", "log_lik", "logit_dem",
                       "RP_pollster", "RP_poll", "delta_dem"), include=F)
print(und.model, pars=c("poll_error", "logit_poll", "log_lik", "logit_dem", "logit_undecided",
                       "RP_pollster", "RP_poll", "delta_dem", "delta_undecided"), include=F)

# plot dem. support estimates
estimates = rstan::extract(intent.model, pars="dem_margin")$dem_margin
mrg = as.data.frame(t(apply(estimates, 2, quantile, probs=c(0.05, 0.5, 0.95))))
names(mrg) = c("low", "median", "high")
mrg$week = 1:nrow(mrg)

ggplot(mrg, aes(x=week)) + geom_line(aes(y=median)) + 
    geom_ribbon(aes(ymin=low, ymax=high), alpha=0.5) +
    ylim(-0.05, 0.15) + geom_hline(yintercept=0, linetype="dashed")
