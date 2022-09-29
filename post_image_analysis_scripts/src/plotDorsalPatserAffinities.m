names = {'Dg', 'DgS', 'DgAW', 'DgW', 'DgSVW', 'DgVW', 'DgVVW'}';
flyregScores = [6.23, 5.81, 5.13, 5.39, 4.80, 4.29, 4.73]';
t = sortrows(table(names, flyregScores), 'flyregScores')
barh(t.flyregScores);
xlim([4, 7])
yticklabels(t.names);
xlabel('Patser score')
title('Affinity of dorsal sites in synthetics, two PWM sources')