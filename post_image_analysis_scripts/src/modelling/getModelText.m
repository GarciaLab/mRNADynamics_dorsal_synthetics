function [coeffText, fits] = getModelText(mdl, fits, k)

    %extract fitting outputs to make the title
    names = string(coeffnames(mdl));
    vals = coeffvalues(mdl)';
    cis = confint(mdl);
    fits(k).names = names;
    fits(k).vals = vals;
    fits(k).cis = cis;
    cis_upper = string(round(cis(2,:)', 2, 'significant'));
    cis_lower = string(round(cis(1,:)', 2, 'significant'));
    cis_upper(ismissing(cis_upper)) = "NaN";
    cis_lower(ismissing(cis_lower)) = "NaN";

%     coeffText = string(coeffnames(mdl)) + " = " +...
%         round(coeffvalues(mdl)',2,'significant') +...
%         " (" + cis_lower + " " + cis_upper + ")";
    
    coeffText = string(coeffnames(mdl)) + " = " +...
        round(coeffvalues(mdl)',2,'significant');
    
end