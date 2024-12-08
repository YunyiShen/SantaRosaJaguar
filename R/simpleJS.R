library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()/2)
raw_data = read.csv("./data/DATOS_MANOLO.csv")

sex_detected = (raw_data$sex == "Female") * 1
raw_y = raw_data[,-(1:2)] |> as.matrix()

targeted_total_inid = 50 # targeted number of individuals for each sex
n_female_aug = targeted_total_inid - sum(sex_detected)
n_male_aug = targeted_total_inid - sum(1 - sex_detected)


n_aug = n_female_aug + n_male_aug
aug_y = matrix(0, n_aug, ncol(raw_y))
sex_aug = rep(0, n_aug)
sex_aug[1:n_female_aug] = 1

y = rbind(raw_y, aug_y)
sex = c(sex_detected, sex_aug)
stan_data = list(T = ncol(y), 
                 N = nrow(y), 
                 sex = sex, 
                 y = y)



m_init <- stan_model("./stan/simple_js.stan")
set.seed(42)
m_fit <- sampling(m_init,  data = stan_data,chains = 4, iter = 2000, 
                  init = function() list(gamma = c(0.2, 0.2), psi = c(0.15, 0.15), 
                                         phi = c(0.85, 0.85), p=c(0.2, 0.2) ), 
                  verbose = TRUE)

save(m_fit, stan_data,
     file="./res/simple_js_stan_fit.rda")


