function y = entryAnalytical(dl, theta, t_cycle)

c = theta(1);
kd = theta(2);
nentries = theta(3);
moffs = theta(4);
pientry = theta(5);

occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );

pioff = c*occupancy(dl, kd); 

onset_pdf = @(t, p, rho, n, m)  exp(-p.*t) .* (p.^m) .* t.^(m+n-1) .* rho.^n .*...
    (gamma(m+n)).^(-1) .*...
    hypergeom( n, m+n, t.*(p-rho) );

factive = integral(@(t) onset_pdf(t, pientry, pioff, nentries, moffs), 0, t_cycle, 'ArrayValued', true);

onset_pdf_trunc = @(t, p, rho, n, m, T) onset_pdf(t, p, rho, n, m) ./ factive;

MFPT_trunc = integral( @(t) t.*onset_pdf_trunc(t, pientry, pioff, nentries, moffs, t_cycle) , 0, t_cycle, 'ArrayValued', true);


y = [factive, MFPT_trunc];

end