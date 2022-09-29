function coeffText = getDistributionText(mdl)

%extract fitting outputs to make the title
    names = string(mdl.ParameterNames);
    vals = mdl.ParameterValues';
    

%     coeffText = string(coeffnames(mdl)) + " = " +...
%         round(coeffvalues(mdl)',2,'significant') +...
%         " (" + cis_lower + " " + cis_upper + ")";
    
    coeffText = names + " = " +...
        round(vals',2,'significant');