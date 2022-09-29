function xrange = getXRange(enhancers, expmnt)
    
    xrange = nan(length(enhancers), 2);
%     xrange(1,:)  = [1000, 2250];  %1Dg
   xrange(1,:)  = [500, 2750];  %1Dg
    if strcmpi(expmnt, 'affinities')
        xrange(2, :) = [500, 2250]; %upper limit for %1DgS
        xrange(3, 1) = 500; %lower limit for 1DgW
%         xrange(4, 2) = 1500; %upper limit for %1DgAW
    xrange(4, 2) = 2250; %upper limit for %1DgAW
    end
    
end

