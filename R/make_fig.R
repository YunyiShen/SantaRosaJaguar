library(ggplot2)
library(rstan)
library(ggpattern)

########## load data and results ###################
raw_data = read.csv("./data/DATOS_MANOLO.csv")
years = 2009:2023
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

load("./res/simple_js_stan_fit.rda")

############## make plot ################
library(reshape2)
z = extract(m_fit, "z")$z
females_z = z[,sex == 1,]
males_z = z[,sex == 0,]

total_N = apply(z, c(1,3), function(w){sum(w==2)}) |> as.data.frame()
colnames(total_N) = years
total_N= melt(total_N)

females_N = apply(females_z, c(1,3), function(w){sum(w==2)}) |> as.data.frame() # total number each year
colnames(females_N) = years
females_N = melt(females_N)
females_N$sex = "Female"

males_N = apply(males_z, c(1,3), function(w){sum(w==2)}) |> as.data.frame() # # total number each year
colnames(males_N) = years
males_N = melt(males_N)
males_N$sex = "Male"

NN = rbind(females_N, males_N)
NN_mean <- aggregate(value~variable + sex, data = NN, FUN = median)

total_N_mean <- aggregate(value~variable, data = total_N, FUN = median)

jpeg("./res/Figs/size.jpg", width = 6, height = 3, units = "in",res = 500)
ggplot(NN, aes(x = variable, y = value, pattern = sex)) + 
  geom_boxplot_pattern(pattern_density = 0.7, # Adjust density of the pattern
                       pattern_angle = 45, 
                       pattern_size = .05,
                       pattern_scale = 0.1,
                       outlier.size = .5) + 
  scale_pattern_manual(values = c("stripe", "circle")) + #
  #geom_boxplot() +
  theme_classic() + 
  xlab("") + 
  ylab("# individuals") + 
  geom_line(aes(group = sex),data = NN_mean) + 
  theme(legend.position = "top")
dev.off()


jpeg("./res/Figs/total_size.jpg", width = 6, height = 3, units = "in",res = 500)
ggplot(total_N, aes(x = variable, y = value)) + 
  geom_boxplot(outlier.size = .5) + 
  #geom_boxplot() +
  theme_classic() + 
  xlab("") + 
  ylab("# individuals") 
dev.off()

### vital rates ###
params = extract(m_fit, c("psi", "phi", "gamma", "p"))

jpeg("./res/Figs/vital_rates.jpg", width = 6, height = 3, units = "in",res = 500)

par(mfcol = c(2,4),mar = c(2.5,2.5,1.7,.5), mgp = c(1.5, 0.5, 0))
plot(density(params$psi[,1]*50), main = "Initial recruitment", xlab = "", ylab = "Female")
polygon(density(params$psi[,1]*50), col = "#9b9b9b")
plot(density(params$psi[,2]*50), main = "", xlab = "", ylab = "Male")
polygon(density(params$psi[,2]*50), col = "#9b9b9b")

plot(density(params$gamma[,1]*50), main = "Recruitment", xlab = "", ylab = "")
polygon(density(params$gamma[,1]*50), col = "#9b9b9b")
plot(density(params$gamma[,2]*50), main = "", xlab = "", ylab = "")
polygon(density(params$gamma[,2]*50), col = "#9b9b9b")

plot(density(params$phi[,1]), main = "Survival", xlab = "", ylab = "")
polygon(density(params$phi[,1]), col = "#9b9b9b")
plot(density(params$phi[,2]), main = "", xlab = "", ylab = "")
polygon(density(params$phi[,2]), col = "#9b9b9b")

plot(density(params$p[,1]), main = "Detection", xlab = "", ylab = "")
polygon(density(params$p[,1]), col = "#9b9b9b")
plot(density(params$p[,2]), main = "", xlab = "", ylab = "")
polygon(density(params$p[,2]), col = "#9b9b9b")

dev.off()





