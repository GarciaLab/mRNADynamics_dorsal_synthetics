function fixSchnitzInconsistencies(prefix)

%     prefix = '2020-06-27-1DgS2_2xDl_2'; %test prefix
    TrackNuclei(prefix, 'retrack')
    TrackmRNADynamics(prefix)
    CompileParticles(prefix,  'minBinSize', 0, 'MinParticles', 0,...
            'yToManualAlignmentPrompt');
    alignCompiledParticlesByAnaphase(prefix);
    
end
