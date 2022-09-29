data {
int<lower=0> N; // observation counter
int<lower=0> K; //param counter
real x[N];
real Y[N];
real p0[K];
real lb[K];
real ub[K];
}
parameters {
real<lower=0> w;
real<lower=0>KD;
real<lower=0> tau;
}
transformed parameters {
real sigma;
real m[N];
for (i in 1:N) 
m[i] =(((x[i]/KD)*w) /(1+ (x[i] /KD) + ((x[i] /KD) *w)));
sigma = 1 / sqrt(tau); 
}
model {
// priors 
w ~ normal(0.0, 10); 
KD ~ normal(500, 10000); 
// likelihood
Y ~ normal(m, sigma); 
}
generated quantities{ 
real Y_mean[N]; 
real Y_pred[N]; 
for(i in 1:N){ 
// Posterior parameter distribution of the mean 
Y_mean[i] =(((x[i]/KD)*w) /(1+ (x[i] /KD) + ((x[i] /KD) *w)));
// Posterior predictive distribution 
Y_pred[i] = normal_rng(Y_mean[i], sigma);   
}
}

