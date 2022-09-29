l = 0.5; 
nSims = 10^9;
gen = {'twister', 'simdTwister', 'combRecursive', 'multFibonacci', 'philox', 'threefry'};
dt = [];

a = [];
for k = 1:length(gen)
    rng(1, gen{k})
    tict
    exprnd(0.5, [1 nSims]); 
    dt(k) = toc
end
