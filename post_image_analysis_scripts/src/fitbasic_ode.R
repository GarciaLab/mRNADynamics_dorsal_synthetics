functions {
  
  vector fraction_fun(real[] x, real c, real kd, real tcycle) {
    //y6State = get_y6(c, kd, tcycle)
    //real timeStepsPerMin = 10;
    //endOfCycle = tcycle*timeStepsPerMin;
    real n = 5;
    int nx = num_elements(x);
    real alpha = c*kd*tcycle;
    vector[nx] out;
    for (i in 1:nx) {
      real arg = alpha ./ (x[i] + kd);
      out[i] = 1 - gamma_q(n, arg);
    }
    return out;
  }
  
  vector onset_fun(real[] x, real c, real kd, real tcycle) {
    real n = 5;
    int nx = num_elements(x);
    vector[nx] out;
    real alpha = c*kd*tcycle;
    for (i in 1:nx) {
      real arg = alpha ./(x[i]+kd);
      out[i] = ( (x[i]+kd) ./ (c .*x[i]) ) .* (tgamma(n+1) .*(1-gamma_q(n+1, arg))) ./ (tgamma(n) .*(1-gamma_q(n, arg)));
    }
    return out;
  }
}


data {
  int<lower=0> N; // observation counter
  int<lower=0> K; //param counter
  real x[N];
  real fraction[N];
  real onset[N];
  //params: c, kd, tcycle
  real p0[K];
  real lb[K];
  real ub[K];
}



parameters {
  real<lower=lb[1], upper=ub[1]>c;
  real<lower=lb[2], upper=ub[2]>kd;
  real<lower=lb[3], upper=ub[3]>tcycle;
  real<lower=0> sigma_fraction;
  real<lower=0> sigma_onset;
}



model{
  //priors
  c ~ normal(p0[1], ub[1]);
  kd ~ normal(p0[2], ub[2]);
  tcycle ~ normal(p0[3], ub[3]);
  
  // sigma_fraction ~ cauchy(0, 1);
  // sigma_onset ~ cauchy(0, 1);
  
  //likelihoods
  fraction ~ normal(fraction_fun(x, c, kd, tcycle), sigma_fraction);
  onset ~ normal(onset_fun(x, c, kd, tcycle), sigma_onset);
}


generated quantities{ 
  
  vector[N] fraction_mean; 
  vector[N] fraction_pred; 
  vector[N] onset_mean; 
  vector[N] onset_pred; 
  
  fraction_mean = fraction_fun(x, c, kd, tcycle);
  onset_mean = onset_fun(x, c, kd, tcycle);
  
  for(i in 1:N){ 
    fraction_pred[i] = normal_rng(fraction_mean[i], sigma_fraction); 
    onset_pred[i] = normal_rng(onset_mean[i], sigma_onset); 
  }
}
