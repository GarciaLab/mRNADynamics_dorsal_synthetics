functions {
    real simpleweak(real x, real KD, real w, real R) {
        return R*(((x ./ KD) .* w) ./ (1+ (x  ./ KD) + ((x ./ KD) .* w)));
    }
    vector simpleweakvec(real[] x, int[] dsid, vector KDs, real w, real R) {
        int nx = num_elements(x);
    //    int n = 1;
        vector[nx] out;
    //    int dsid[nx];
    //    dsid[1] = 1;
    //    for (i in 2:nx){
    //      if (x[i] < x[i-1]){
    //        n = n + 1;
    //      }
    //      dsid[i] = n;
    //    }
    for (i in 1:nx) {
      int dsidi = dsid[i];
      out[i] = R *(((x[i] ./ KDs[dsidi]) * w) ./ (1+ (x[i]  ./ KDs[dsidi]) + ((x[i] ./ KDs[dsidi]) * w)));
    }  
    return out;
    }
}
data {
int<lower=0> N; // observation counter
int<lower=0> K; //param counter
real x[N];
real Y[N];
int dsid[N];
//params: w, kd1, kd2, kd3, kd4, kd5, kd6, kd7, R
real p0[K];
real lb[K];
real ub[K];
}
parameters {
real<lower=lb[1], upper=ub[1]>w;
real<lower=lb[2], upper=ub[2]>KD1;
real<lower=lb[3], upper=ub[3]>KD2;
real<lower=lb[4], upper=ub[4]>KD3;
real<lower=lb[5], upper=ub[5]>KD4;
real<lower=lb[6], upper=ub[6]>KD5;
real<lower=lb[7], upper=ub[7]>KD6;
real<lower=lb[8], upper=ub[8]>KD7;
real<lower=lb[K], upper=ub[K]>R;
real<lower=0> sigma;
}
model {
vector[7] KDs;
KDs[1] = KD1;
KDs[2]=KD2;
KDs[3]=KD3;
KDs[4]=KD4;
KDs[5]=KD5;
KDs[6]=KD6;
KDs[7]=KD7;
// priors 
w ~ normal(p0[1], p0[1]); 
KD1 ~ normal(p0[2], p0[2]);
KD2 ~ normal(p0[3], p0[3]); 
KD3 ~ normal(p0[4], p0[4]); 
KD4 ~ normal(p0[5], p0[5]); 
KD5 ~ normal(p0[6], p0[6]); 
KD6 ~ normal(p0[7], p0[7]); 
KD7 ~ normal(p0[8], p0[8]); 
R ~ normal(p0[K], p0[K]); 
// likelihood
Y ~ normal(simpleweakvec(x, dsid, KDs, w, R), sigma); 
}

generated quantities{ 
vector[N] Y_mean; 
vector[N] Y_pred; 
vector[7] KDs;
KDs[1] = KD1;
KDs[2]=KD2;
KDs[3]=KD3;
KDs[4]=KD4;
KDs[5]=KD5;
KDs[6]=KD6;
KDs[7]=KD7;
Y_mean = simpleweakvec(x,dsid, KDs, w, R);

for(i in 1:N){ 
// Posterior parameter distribution of the mean 
//Y_mean[i] = R*(((x[i]/KD)*w) /(1+ (x[i] /KD) + ((x[i] /KD) *w)));
// Posterior predictive distribution 
Y_pred[i] = normal_rng(Y_mean[i], sigma);   
}
}