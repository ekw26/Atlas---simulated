library(numbers)
library(zoo)
library(ggplot2)
library(MASS)

# n patients
n_pats <- 1000
ids <- seq(1, n_pats)

# patient-level clustering
# this is random noise included to make the default assumptions 'wrong'
p_cluster <- rnorm(n_pats, 0, 0.02)


#want each patient to have a random start month and year
#more likely to be later years
start_year <- 2019 - rgeom(n_pats, 0.25)
start_month <- sample(1:12, n_pats, replace = T)


# calculate each patients max follow up
# assume dec 2019 is final possible month
max_follow_up <- (2020 - start_year)*12 - start_month + 1

# randomly assign each person follow up length between 1 and their maximum
rfu <- function(m) sample(1:m, 1)
rand_follow_up <- sapply(max_follow_up, rfu, simplify = T)

data <- data.frame(id = rep(ids, rand_follow_up), p_cluster = rep(p_cluster, rand_follow_up), months_pre_diag = sequence(rand_follow_up), diag_year = rep(start_year, rand_follow_up), diag_month = rep(start_month, rand_follow_up))

#work out current year and current month for each row
# add into simulation formula

get_current_date <- function(diag_month, diag_year, n_months) {
  current_month <- diag_month
  current_year <- diag_year
  for (i in 1:n_months) {
    current_month <- current_month - 1
    if (current_month == 0) {
      current_month <- 12
      current_year <- current_year - 1
    }
  }
  
  return(list("year" = current_year, "month" = current_month))
}

current_dates <- as.data.frame(t(mapply(get_current_date, data$diag_month, data$diag_year, data$months_pre_diag)))
data$current_month <- as.numeric(current_dates$month)
data$current_year <- as.numeric(current_dates$year)
data$current_date <- as.yearmon(paste(data$current_year, data$current_month), "%Y %m")

# generate a new variable which tells us when the month of inflection point is
inf_point <- 4
data$alt_month <- 0
data[data$months_pre_diag <= inf_point, 'alt_month'] <- inf_point - data[data$months_pre_diag <= inf_point, 'months_pre_diag'] + 1

# generate a person-month variable for consistency
data$py <- 1

# simulate some data based on round-ish numbers
# follow a poisson distribution based around the patient-level mean defined by exp of variables
data$n_consultations <- rpois(nrow(data), exp(
  -1
  + 1e-3*0.05*data$current_year
  - 0.02*data$months_pre_diag
  + 0.2*data$alt_month
  + rnorm(nrow(data), 0, 0.05)
  + data$p_cluster
  + log(data$py)
))
# add age 

#plot it
g <- ggplot(data, aes(months_pre_diag, n_consultations))
g + geom_jitter(col = "darkgreen", alpha = 0.3) + 
  geom_smooth(method = "glm", method.args = list(family = "poisson")) + 
  geom_smooth(method = "glm.nb", col = "purple") +
  labs(title = "Number of consultations per month", x = "Months before diagnosis", y = "Number of consultations") + 
  scale_x_reverse()

ggplot(data, aes(current_date, n_consultations)) + geom_jitter(col = "darkred", alpha = 0.3)

pois <- glm(n_consultations ~ current_year + months_pre_diag, data = data, family = "poisson")
neg_bin <- glm.nb(n_consultations ~ current_year + months_pre_diag, data = data)
pframe <- with(data, expand.grid(months_pre_diag = seq(1, 300), current_year = seq(1968, 2019)))
pframe$n_consultations_pois <- predict(pois, newdata = pframe, type = "response")
pframe$n_consultations_neg_bin <- predict(neg_bin, newdata = pframe, type = "response")

ggplot(data, aes(months_pre_diag, n_consultations)) + geom_jitter(col = "lightgreen", alpha = 0.3) + geom_line(aes(months_pre_diag, n_consultations_pois, col = current_year), data = pframe)
