functions {
  
  vector diff(vector v){
    
    int n = num_elements(v);
    vector[n-1] out;
    
    for (i in 1:n-1){
      out[i] = v[i+1] - v[i];
    }
    
    return out;
  }
  
  
  matrix getFraction_Onset(vector dorsalVals, real c, real kd,  real tcycle) {
    
    int numCells = 100;   // the total number of nuclei in the simulation
    real transcriptionStart = 0.01; //minutes, delayed start of the transcriptional window
    real dt = tcycle/80; //this 80 seems sufficient for all purposes. sorry for hardcoding.
    int NOffStates = 4;   //number of off states
    int NInactiveStates = 0; //number of inactive states
    int NLinStates = NOffStates+NInactiveStates;
    int nt = 79;
    int time_vec[nt]; // just a vector to loop over
    vector[nt+1] time_vec_2; // just a vector to loop over
    int transcriptionStartRow = 1;
    int nM2 = NLinStates+1;
    matrix[nt+1, nM2] M;
    int nd = num_elements(dorsalVals);
    vector[nd] kdt_off;
    vector[nt+1] yEnd;
    vector[nd] onset;
    vector[nd] fraction;
    matrix[nd, 2] onset_fraction;
    
    for (i in 1:nt){
      time_vec[i] = i+1;
    }
    
    for (i in 0:nt){
      time_vec_2[i+1] = dt*i;
    }
    
    for (i in 1:nt+1){
      for (j in 1:nM2){
        M[i,j] = 0; //initialize it to zero everywhere
      }
    }
    //Initial conditions: everyone is at the first state at the first dt.
    //note that the system enters state one with a delay specified by "transcriptionStart"
    M[transcriptionStartRow,1]=numCells; //
    
    
    kdt_off = (c*(dorsalVals ./kd) ./ (1 + dorsalVals ./kd))*dt;
    
    for (d in 1:nd){
      for (t in time_vec[transcriptionStartRow:nt]){ // loop over time steps        
      //dls = dorsalTraceFluo[idx[t-1]]; // each dorsal bin has a corresponding concentration time trace
      //dls = dorsalVals(d) + diff(dorsalVals(1:2)); // this is in case we want constant Dorsal       
      //  kdt_off = (c*(dls./kd) ./ (1 + dls./kd))*dt; // transition rate between off states
      
      // [~,nearestBin] = min(abs(modelOpts.middleBinValues - dorsalVals(d)));
      //    dorsalTraceFluo = modelOpts.TimeVariantDorsalValues(:,nearestBin);
      //Calculate the first state
      
      M[t,1] = (1-kdt_off[d])*M[t-1,1]; // it transitions with a rate of kdt_off
      
      
      // loop over the rest of the off states
      for (s in 2:NLinStates){
        M[t,s] = (1-kdt_off[d])*M[t-1,s] + kdt_off[d]*M[t-1,s-1];
      }
      
      //Calculate the last state
      M[t,nM2] = M[t-1,nM2] + kdt_off[d]*M[t-1,nM2-1];
      
      }
      
      
      yEnd = M[:,nM2]; //number of nuclei in the last state as a function of time
      
      onset_fraction[d, 1] = sum(diff(yEnd) .* time_vec_2[1:nt]) ./ sum(diff(yEnd)); //expected value
      
      onset_fraction[d, 2] = M[nM2,nM2]/numCells;
      
    }
    
    return onset_fraction;
  }
  
  
}


data {
  int<lower=0> N; // observation counter
  int<lower=0> K; //param counter
 // real x[N];
  vector[N] x;
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
  
 fraction ~ normal(getFraction_Onset(x, c, kd, tcycle)[:, 1], sigma_fraction);
 onset ~ normal(getFraction_Onset(x, c, kd, tcycle)[:, 1], sigma_fraction);
}


generated quantities{ 
  
  vector[N] fraction_mean; 
  vector[N] fraction_pred; 
  vector[N] onset_mean; 
  vector[N] onset_pred; 
  
  fraction_mean = getFraction_Onset(x, c, kd, tcycle)[:, 1];
  onset_mean = getFraction_Onset(x, c, kd, tcycle)[:, 2];
  
  for(i in 1:N){ 
    fraction_pred[i] = normal_rng(fraction_mean[i], sigma_fraction); 
    onset_pred[i] = normal_rng(onset_mean[i], sigma_onset); 
  }
}
