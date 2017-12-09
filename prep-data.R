# DATA PREPARATION 

library(pollstR)
library(tidyverse)
library(lubridate)

# Past Presidential approval
president = read.csv("data/past_pres_approval_raw.csv", sep="\t") %>%
    transmute(approval=Approving / 100,
              date=mdy(End.Date),
              president=President,
              party=Party) %>%
    filter(year(date) >= 1980, year(date) <= 2030) %>%
    arrange(date)
write.csv(president, "data/past_pres_approval.csv", row.names=F)

# Current polling
polls = pollster_charts_polls("2018-national-house-race")[["content"]] %>%
    filter(partisanship == "Nonpartisan") %>%
    transmute(dem = Democrat / 100,
              gop = Republican / 100,
              other = Other / 100,
              undecided = Undecided / 100,
              n_resp = observations,
              n_side = (dem + gop) * n_resp,
              n_dem = dem * n_resp,
              pollster = survey_house,
              date = ymd(end_date),
              type = recode(sample_subpopulation, `Registered Voters`="RV",
                            `Likely Votesrs`="LV")) %>%
    as.data.frame
write.csv(polls, "data/current_polls.csv", row.names=F)

# past polling
past.polls = read.csv("data/past_generic_polling.tsv", sep="\t")

# Merge
economy = read.csv("data/economy.csv", na.strings=".") %>%
    transmute(date = ymd(DATE),
              unemp = as.numeric(UNRATE) / 100,
              gdp = as.numeric(A939RX0Q048SBEA_PC1) / 100,
              infl = as.numeric(CPIAUCSL_PC1) / 100,
              earn = as.numeric(AHETPI_PC1) / 100) %>%
    filter(year(date) > 1980)
congress = read.csv("data/congress_approval.csv")
control = read.csv("data/party_control.csv")
