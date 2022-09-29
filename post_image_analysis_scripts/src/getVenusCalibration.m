function [ap, venus] = getVenusCalibration(schnitzcells)

ap = [];
venus = []; 
GREGOR_TIME_POINT = 15; %min

% load absolute Bcd concentration data from Gregor 2007
bcd_abs_path = 'S:\Armando\Dropbox\DorsalSyntheticsDropbox\venus_calibration\';
bkg_data = readtable([bcd_abs_path 'Gregor2007Black.csv']);
bkg_data_clean.AP = bkg_data.Var1;
bkg_data_clean.nM = bkg_data.Var2;
bcd_data01 = readtable([bcd_abs_path 'Gregor2007BcdRed.csv']);
bcd_data02 = readtable([bcd_abs_path 'Gregor2007BcdBlue.csv']);
bcd_data_clean.AP = [bcd_data01.Var1 ; bcd_data02.Var1];
bcd_data_clean.nM = [bcd_data01.Var2 ; bcd_data02.Var2];
bcd_data_clean.setID = [repelem(1,numel(bcd_data01.Var1)), repelem(2,numel(bcd_data02.Var1))];

for s = 1:length(schnitzcells)
    if schnitzcells(s).cycle == 14
        t = schnitzcells(s).timeSinceAnaphase;
        for f = 1:length(t)
        
            if GREGOR_TIME_POINT - t(f) <= .2 &&...
                    length(schnitzcells(s).FluoTimeTrace) ==...
                    length(schnitzcells(s).timeSinceAnaphase)
                
                ap = [ap, schnitzcells(s).APPos(f)];
                venus = [venus, schnitzcells(s).FluoTimeTrace(f)];
                
            end
        end
            
        
    end
end

figure;
yyaxis left;
scatter(ap, venus, 'g')
ylabel('Bcd-Venus (au)')
hold on
yyaxis right;
scatter(bcd_data_clean.AP, bcd_data_clean.nM, 'k')
ylabel('Bcd concentration (nM)')
xlim([.1, .25])
xlabel('Fraction AP');
ap = double(ap);
venus = double(venus);

dat = table(ap, venus);
dat = sortrows(dat);
p = polyfit(dat.ap, dat.venus, 1)

end