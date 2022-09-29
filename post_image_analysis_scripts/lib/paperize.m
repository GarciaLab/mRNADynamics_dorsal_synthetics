function [ ] = paperize(figHandle,width,height,varargin)
% This function modifies a figure to make it print-sized for saving into a
% manuscript file (such as EPS). It also resizes the figure on the screen
% so the user can see how it will look when printed, changes all fonts
% to TNR, and scales fonts/lines/markers proportionally. It is easy to add
% user-defined parameters by copying the style of the set(findall(...))
% commands.
%
% REQUIRED ARGUMENTS:
% figHandle - the handle of the figure to be modified.
% width - 'a4Half' and 'a4Full' are for single-column and cross-column
%           (8.2/17.2 cm) figures on a portrait a4 sheet. 'letterHalf' and
%           'letterFull' are same, for letter paper (3 9/16" and 7 1/2").
%           User can also specify width numerically, in cm.
% height - height of the figure in centimeters.
%
% OPTIONAL ARGUMENTS:
% percent scaling for plot or subplots (e.g., 110 -> 10% growth)
% this may be useful for shrinking margins when desired.
%
% (c) 2017 - Daniel D. Borup
%

if(nargin > 3)
    margScale = varargin{1};
end
figHandle.Units = 'centimeters';
figHandle.PaperUnits = 'centimeters';

if(strcmp(width,'a4Half'))
    width = 8.2;
    fontSize = 10;
    lineThk = 1.5;
    markerSize = 4;
    
elseif(strcmp(width,'a4Full'))
    width = 17.2;
    fontSize = 11;
    lineThk = 2;
    markerSize = 6;
    
elseif(strcmp(width,'letterFull'))
    width = 2.54*(7.5);
    fontSize = 11;
    lineThk = 2;
    markerSize = 6;
    
elseif(strcmp(width,'letterHalf'))
    width = 2.54*(3.5625);
    fontSize = 10;
    lineThk = 1.5;
    markerSize = 4;
    
elseif(isnumeric(width) && width > 3)
    slope = (width - 8.4)/9;
    fontSize = 10 + slope;
    lineThk = 1.5 + slope/2;
    markerSize = 4 + 2*slope;

else
    disp(['Error - size must be a valid paper specifier',...
        'or a width in centimeters (>3).']);
    return
end

figHandle.Position = [10 10 width height];
figHandle.PaperPosition = [10 10 width height];
figCenter = [width, height]./2;
set(findall(figHandle,'-property','Interpreter'),'Interpreter','LaTeX');
set(findall(figHandle,'-property','FontSize'),'FontSize',fontSize);
set(findall(figHandle,'-property','FontName'),'FontName','Arial');
set(findall(figHandle,'-property','LineWidth'),'LineWidth',lineThk);
set(findall(figHandle,'-property','MarkerSize'),'MarkerSize',markerSize);
if(nargin > 3)
    subHandles = get(figHandle,'children');
    for i=1:length(subHandles)
        sh = subHandles(i);
        if(strcmp(sh.Tag,''))
            sh.Units = 'centimeters';
            currPos = sh.Position;
            currPos(1:2) = currPos(1:2) - figCenter;
            scaleMat = currPos * (margScale/100 - 1);
            sh.Position = sh.Position + scaleMat.*[.5 .5 1 1];
        end
    end
end

end