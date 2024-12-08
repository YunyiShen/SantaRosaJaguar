functions {
  vector log_sum_exp_helper(matrix x, int N){
    vector[N] res;
    for (i in 1:N){
      res[i] = log_sum_exp(x[i,]);
    }
    return res;
  }
}


data {
  int<lower=1> T; // number of time steps
  int<lower=1> N; // number of individuals
  matrix<lower=0, upper=1>[N, T] y; // observations
  vector<lower=0, upper=1>[N] sex; // sex of the individuals
}

parameters {
    real<lower = 0, upper = 1> psi[2]; // initial recruitment for both sexes
    real<lower = 0, upper = 1> phi[2]; // survival probability
    real<lower = 0, upper = 1> gamma[2]; // probability being recruited
    real<lower = 0, upper = 1> p[2]; //detection probability
}

transformed parameters { // calculate likelihood, with forward 
    matrix[N,T] log_p_not_recruited; // likelihood for individual being not recruited
    matrix[N,T] log_p_alive; // likelihood for individual being alive
    matrix[N,T] log_p_dead;  // likelihood for individual being dead
    matrix[N,3] logLikl; // at the end likelihood for individual being in each state
    matrix[N,2] acc2;
    matrix[N,2] acc3;
    
    real negINF = -1e10;
    vector[N] logdetprob = log(p[1] * sex + p[2] * (1-sex));
    vector[N] log1mdetprob = log1m(p[1] * sex + p[2] * (1-sex));

    vector[N] logrecruitprob = log(gamma[1] * sex + gamma[2] * (1-sex));
    vector[N] log1mrecruitprob = log1m(gamma[1] * sex + gamma[2] * (1-sex));

    vector[N] logsurviveprob = log(phi[1] * sex + phi[2] * (1-sex));
    vector[N] log1msurviveprob = log1m(phi[1] * sex + phi[2] * (1-sex));

    // handel first time step
    for (i in 1:N){
    log_p_not_recruited[, 1] = log1m(psi[1]) * sex + log1m(psi[2]) * (1 - sex);
    log_p_not_recruited[, 1] = log_p_not_recruited[, 1] + negINF * (y[, 1]);
    log_p_alive[, 1] = log(psi[1]) * sex + log(psi[2]) * (1-sex);
    log_p_alive[, 1] = log_p_alive[, 1] + logdetprob .* y[, 1] + log1mdetprob .* (1 - y[, 1]);
    log_p_dead[,1] = negINF * (1+sex);

    for (tt in 2:T) {
        // land at not recruited
        log_p_not_recruited[, tt] = log1mrecruitprob + log_p_not_recruited[, tt-1] + negINF * (y[, tt]);
        // land at alive
            // being recruited
        acc2[,1] = logrecruitprob + log_p_not_recruited[, tt-1] + logdetprob .* y[, tt] + log1mdetprob .* (1 - y[, tt]);
            // surviving
        acc2[,2] = logsurviveprob + log_p_alive[, tt-1] + logdetprob .* y[, tt] + log1mdetprob .* (1 - y[, tt]);
        log_p_alive[, tt] = log_sum_exp_helper(acc2, N);

        // land at dead
            // dead
        acc3[,1] = log1msurviveprob + log_p_alive[, tt-1] + negINF * y[, tt];
            // dead already last time
        acc3[,2] = log_p_dead[, tt-1] + negINF * y[, tt];
        log_p_dead[, tt] = log_sum_exp_helper(acc3,N);
    }

      logLikl[i, 1] = log_p_not_recruited[i, T];
      logLikl[i, 2] = log_p_alive[i, T];
      logLikl[i, 3] = log_p_dead[i, T];
    }


}

model {
    // priors
    psi ~ beta(1, 1);
    phi ~ beta(1, 1);
    gamma ~ beta(1, 1);
    p ~ beta(1, 1);
    // likelihood
    for (n in 1:N) {
        target += log_sum_exp(logLikl[n]);
    }
}

generated quantities {
    matrix[N, T] z;// hidden states
    // forward filtering backward sampling
    {
        int t;
        vector[3] bkw;
        for (i in 1:N){ // each individual
            z[i,T] = categorical_rng(softmax(logLikl[i, ]')); // last time step
            for (tt in 2:T){ // backward sampling
                t = T-tt + 1;
                if (z[i,t+1] == 1){
                    bkw[1] = log1mrecruitprob[i];
                    bkw[2] = negINF;
                    bkw[3] = negINF;
                } else if (z[i,t+1] == 2){
                    bkw[1] = logrecruitprob[i];
                    bkw[2] = logsurviveprob[i];
                    bkw[3] = negINF;
                } else {
                    bkw[1] = negINF;
                    bkw[2] = log1msurviveprob[i];
                    bkw[3] = 0;
                }
                z[i,t] = categorical_rng(softmax(bkw));
            
            }

        }


    }


}
