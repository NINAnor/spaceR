# code to generate dummy dataset 'mypop'

library(tidyverse)

s1 <- s2 <- 1:60 - 0.5

# set the per strata sample size
n <- 5

# expand to get every combination
# These will act as coordinates for our dummy data
mypop <- expand.grid(s1, s2)

# name columns
names(mypop) <- c("Easting", "Northing")

# add column 'stratum'
# Strata '0' for when s1 is between -Inf and 4;
# Strata '1' for when s1 is between 4 and 6; etc.
# mypop$stratum <- as.factor(findInterval(mypop$Easting, c(4, 6, 14)))
mypop$stratum <- as.factor(findInterval(mypop$Easting, c(12, 18, 42)))

# give the strata names as capital letters
levels(mypop$stratum) <- LETTERS[1:4]

# add a column for area that can be used for pps
mypop <- mypop |>
  mutate(
    area = case_when(
      Northing > 30 ~ 1,
      .default = 0.3
    ),
    area2 = case_when(
      Northing > 10 ~ 1,
      .default = 0.5
    )
  )

mypop <- mypop |>
    # get population size N,
    # and the population area A_h,
    # for each strata
  group_by(stratum) |>
  summarise(
    N_h = n(),
    A_h = sum(area),
    A_h2 = sum(area2),
    .groups = "drop"
  ) |>
    # calculate inclusion probability (pi) per strata, assuming we will sample 5 units
  mutate(pi_h = rep(n, 4) / N_h) |>
    # join with the original data frame
  right_join(mypop, by = "stratum") |>
    # adding also a pi based on the total area for each strata
  mutate(
    pi_area = n * area / A_h,
      # add pi with n = 40
    pi_h2 = rep(40, 4) / N_h,
      # add unique ID's to each population unit
    ID = row_number(),
      # add mean values for response variables
    res1_mean = case_when(
      stratum == "A" ~ 5,
      stratum == "B" ~ 10,
      stratum == "C" ~ 15,
      stratum == "D" ~ 20
    ),
    res2_mean = case_when(
      stratum == "A" ~ 50,
      stratum == "B" ~ 10,
      stratum == "C" ~ 30,
      stratum == "D" ~ 14,
    )
  ) |>
    # create dummy response variables
  group_by(stratum) |>
  mutate(
    res1 = rnorm(n = n(), mean = res1_mean, sd = 1),
    res2 = rnorm(n = n(), mean = res2_mean, sd = 3),
    # add dummy auxiliary variables
    aux1 = seq(0.5, 5.5, length.out = n())^2,
    aux2 = rgamma(n = n(), shape = .4)
  ) |>
  ungroup() |>
    # round values (but don't round the probabilities)
  mutate(across(starts_with(c("aux", "res")), \(x) round(x, 3))) |>
  dplyr::select(-res1_mean, -res2_mean)

usethis::use_data(mypop, overwrite = TRUE)
