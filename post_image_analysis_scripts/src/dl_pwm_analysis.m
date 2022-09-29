function dl_pwm_analysis()

dbs623 = [4 4 4 1 1 1 1 2 2 2];
deltaG_dbs623 = getEnergy(dbs623);
kd623 = exp(-deltaG_dbs623)

dbs581 = [4 4 4 1 1 1 3 2 2 2];
deltaG_581 = getEnergy(dbs581);
kd581= exp(-deltaG_581)

dbs429 = [1 4 4 1 1 1 3 2 2 1];
deltaG_429 = getEnergy(dbs429);
kd429 = exp(-deltaG_429)

end

function deltaG = getEnergy(seq)

d = [0.071429     0.21429    0.095238     0.80952     0.83333     0.88095     0.59524     0.11905     0.21429     0.21429
0.21429     0.11905     0.02381    0.047619    0.071429     0.02381    0.047619     0.30952     0.57143     0.59524
0.16667     0.02381    0.047619     0.02381    0.071429    0.047619     0.28571     0.47619     0.14286    0.071429
0.54762     0.64286     0.83333     0.11905     0.02381    0.047619    0.071429    0.095238    0.071429     0.11905];

pb = [.285 .215 .285 .215]; %a c t g droso genome probs

for j = 1:size(d, 1)
    for k = 1:size(d, 2)
        H(j, k) = d(j, k) ./ pb(j);
    end
end


f = [];
for k = 1:length(seq)
    f (k)=  H (seq(k), k );
end

deltaG = sum(f);

end