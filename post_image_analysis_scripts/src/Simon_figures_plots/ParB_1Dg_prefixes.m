Ch02Prefixes = {
 '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_4ch2filt','2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_5ch2filt',...
  '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_6ch2filt','2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_8ch2filt',...
  '2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_9ch2filt', '2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_10ch2filt',...
   '2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_11ch2filt','2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_18ch2filt',...
'2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_19ch2filt','2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_20ch2filt',...
 '2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_21ch2filt','2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_22ch2filt',...
'2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_26ch2filt','2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_27ch2filt',...
  '2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_30ch2filt'};

% these don't have MS2 spots:
%'2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_28ch2filt','2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_29ch2filt'
% '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_1ch2filt',...
%  '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_2ch2filt', '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_3ch2filt',...

% for p = 1:length(Prefixes)
%     Prefix = Ch02Prefixes{p};
%     
%     %CheckNucleiSegmentation(Prefix,'preLoadMovie')
%     %chooseAnaphaseFrames(Prefix)
%     %segmentSpots(Prefix,5700, 'nWorkers', 8, 'Shadows', 0)
% end

%% now Ch01 
%ResultsFolder = 'S:\Simon\Dropbox\DorsalSyntheticsDropbox\';

for p = 4:length(Ch02Prefixes)
    Prefix = Ch02Prefixes{p};
    Ch01Prefix = Prefix(1:end-7);
    CheckParticleTracking(Ch01Prefix,'multiview', 'preLoadMovie')
%    segmentSpots(Ch01Prefix,[5700, 5700],'nWorkers', 8,'segmentChannel',1, 'Shadows',0)
%     PrefixResultsFolder = [ResultsFolder '\' Prefix];
%     Ch01PrefixResultsFolder = [ResultsFolder '\' Ch01Prefix];
%      load([Ch01PrefixResultsFolder '/Ellipses.mat'])
%      save([PrefixResultsFolder '/Ellipses.mat'],'Ellipses')

    %Prefix = Prefix(1:end-7);
    %CheckParticleTracking(Prefix,'multiview')
    %CheckNucleiSegmentation(Prefix,'preLoadMovie')
    %chooseAnaphaseFrames(Prefix)
end



