function [p0, lb, ub, name] = getInits(expmnt, md, metric, x_max, y_max, nSets, varargin)

minKD = 200;
maxKD = 1E4;
minw = 1E-2; %1E-2
maxw = 1E1; %1E2
w0 = 1;
minR = 10;
maxR = 1E3;
subset = false;

nph = 3;
naff = 7;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


if expmnt == "affinities" && contains(md, "simpleweak") && metric=="fraction" 
    %simplebindingweak_fraction- wDP kd1..kdn
    p0 = [w0; [x_max/2;x_max;maxKD.*ones(1, nSets-2)']];
    lb = [minw; minKD.*ones(1, nSets)'];
    ub = [maxw; maxKD*ones(1, nSets)'];
    name = ["w"; "KD" + num2str((1:nSets)')];
elseif expmnt == "affinities" && contains(md, "simpleweak") && metric=="fluo"
    %simplebindingweak_fluo- wDP kd1..kdn amp off
    if ~subset
        p0 = [w0; [x_max/2;x_max;maxKD.*ones(1, nSets-2)']; 500; 10];
    else
        p0 = [w0; [maxKD.*ones(1, nSets)']; 500; 10];
    end
    lb = [minw; minKD.*ones(1, nSets)'; minR; 10];
    ub = [maxw; maxKD*ones(1, nSets)'; maxR; 1E3];
    name = ["w"; "KD" + num2str((1:nSets)'); "R"; "offset"];
elseif expmnt == "phases" && contains(md, "simpleweak") && metric=="fraction"
     %simplebindingweak_fraction- wDP1...wDPn kd
    p0 = [x_max/2; w0*ones(1, nSets)'];
    lb = [minKD; minw*ones(1, nSets)'];
    ub = [maxKD; maxw*ones(1, nSets)'];
    name = ["KD"; "w" + num2str((1:nSets)')];
elseif expmnt == "phases" && contains(md, "simpleweak") && metric=="fluo"
     %simplebindingweak_fraction- wDP1...wDPn kd amp off
    p0 = [x_max/2; w0*ones(1, nSets)'; 500; 10];
    lb = [minKD; minw*ones(1, nSets)'; minR; 10];
    ub = [maxKD; maxw*ones(1, nSets)'; maxR; 1E3];
    name = ["KD"; "w" + num2str((1:nSets)'); "R"; "offset"];
elseif expmnt == "phaff" && contains(md, "simpleweak") && metric=="fraction"
    %enhancers- 1dg...1dgvw 1dg-5, 1dg-8
    % w1, w2, w3, kd1, kd2...kd7    
    p0 = [w0*ones(1, nph)'; [x_max/2;x_max;maxKD.*ones(1, naff-2)']];
    lb =  [minw*ones(1, nph)'; minKD.*ones(1, naff)'];
    ub = [maxw*ones(1, nph)'; maxKD*ones(1, naff)'];
    name = ["w" + num2str((1:nph)'); "KD" + num2str((1:naff)')];
elseif expmnt == "phaff" && contains(md, "simpleweak") && metric=="fluo"
    %enhancers- 1dg...1dgvw 1dg-5, 1dg-8
    % w1, w2, w3, kd1, kd2...kd7, amp, off
    p0 = [w0*ones(1, nph)'; [x_max/2;x_max;maxKD.*ones(1, naff-2)']; 500; 10];
    lb = [minw*ones(1, nph)'; minKD.*ones(1, naff)'; minR; 10];
    ub = [maxw*ones(1, nph)'; maxKD*ones(1, naff)'; maxR; 1E3];
    name = ["w" + num2str((1:nph)'); "KD" + num2str((1:naff)');  "R"; "offset"];

end

if contains(md, "dimer")
    k = [1000, 1E-3, 1E5];
    p0 = [p0; k(1)];
    lb = [lb; k(2)];
    ub = [ub; k(3)];
    name = [name; "k"];
end

%simple binding- amp kd1 kd2 kd3 p/kdp offset wdp 
% p0 = [y_max; [x_max/2;x_max;1E5*ones(1, nSets-2)']; 1 ;0; 1];
% lb = [0; 200*ones(1, nSets)'; 0; 0; 1];
% ub = [min(1, y_max*2); Inf*ones(1, nSets)'; Inf; y_max*10; Inf];

%hill- amp kd1..kdn n offset
% p0 = [y_max;  [x_max/2;x_max;1E4*ones(1, nSets-2)'] ; 1 ; 1];
% lb = [0; 200*ones(1, nSets)';1; -y_max];
% ub = [y_max*2; 1E4*ones(1, nSets)'; 8; y_max*10];


%simple binding weak- amp kd1 kd2 kd3 offset wdp
% p0 = [y_max; [x_max/2;x_max;1E4*ones(1, nSets-2)']; 0; .5];
% lb = [0; 200*ones(1, nSets)'; 0; 0];
% ub = [min(1, y_max*2); Inf*ones(1, nSets)'; y_max*10; Inf];

end