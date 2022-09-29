simChain <- function(nSims, c, kd, nentries, moffs, pientry, piexit) {
  
set.seed(1);
# nSims = 1E3;
t_cycle = 8; #min
exitOnlyDuringOffStates = TRUE;

dls = seq(10, 4000, by=250)

# nentries = 5;
# moffs = 5;
nSilentStates = 1;
nOffEntryStates = moffs + nentries;
nStates = nentries + moffs+ 1 + nSilentStates;

# pientry = 2;
# piexit = 2;
# c = 10;
# kd = 1000;

tau_entry = rexp(nSims*nentries,pientry)
dim(tau_entry) <- c(nentries, nSims)

tau_exit = rexp(nSims*(nStates-1),piexit)
dim(tau_exit) <- c(nStates-1, nSims)

factive = rep(NA,length(dls) )
onset = rep(NA,length(dls))

for (i in 1:length(dls)){
  
  pioff = c*(dls[i]/kd) / (1 + (dls[i]/kd))
  
  tau_off = rexp(nSims*(moffs+1),pioff)
  dim(tau_off) <- c((moffs+1), nSims)
  
  tau_entry_off = rbind(tau_entry, tau_off)
  
  if(nSilentStates == 1){
    whichTransition = tau_entry_off < tau_exit 
    if(exitOnlyDuringOffStates){
      whichTransition[1:nentries,] = TRUE;
    }
  } else {
    whichTransition = rep(TRUE, length(tau_entry_off));
  }
  
  reachedOn = colSums(whichTransition[1:nOffEntryStates,]) == nOffEntryStates;
 
  onsets_sim = colSums(tau_entry_off[1:nOffEntryStates,])
  
  trunc = onsets_sim < t_cycle;
  
  factive[i] = sum(reachedOn[trunc])/nSims;
  onset[i] =  mean(onsets_sim[trunc]);
  
}

return(cbind(factive, onset));

}