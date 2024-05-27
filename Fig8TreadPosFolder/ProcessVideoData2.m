function [t, x, y, angle] = ProcessVideoData2(file,diode,p)
%diode: 0 - luminance, 1 - doublediode (processed by a separate function),
%2 - red, 3 - green
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Extracted X
fieldSelection(3) = 1; % Extracted Y
fieldSelection(4) = 1; % Extracted Angel
fieldSelection(5) = 1; % Targets
fieldSelection(6) = 1; % Points
extractHeader = 0; % Do we return header 1 = Yes, 0 = No.
extractMode = 1; % Extract all data, 5 different extraction modes, see help file for Nlx2MatVt
        
[t,x,y,angles,targets,points] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
% this code isheritaged from room 5, where x is -x, and y is -y
% here. The code recognize left, right & reward etc set in rrom
% 5, to simply use it, I set x = -x, y = -y here.
% Though I can change code to recognize L,R, Base, reward etc
% by this room4 direction, no meaning to do that since there is
% easy fix
% Li Yuan, Aug-06-2020, USCD
x=-x;
y=-y;
t = t*10^-6;

% Convert data
handles.targets = uint32(targets);

[dTargets,trackingColour] = decodeTargets(handles.targets);

if ~exist('diode','var')||isempty(diode)
    diode = 0;
end
switch diode
    case 1 % double diode, process with a separate function
        [rx,ry,~,~] = extractPosition2(dTargets,2);
        [gx,gy,~,~] = extractPosition2(dTargets,3);
        % Merge colors to single track
        x = (rx+gx)/2; x = x';
        y = (ry+gy)/2; y = y';
    case 2 % Use red signal
        [x,y,~,~] = extractPosition2(dTargets,2);
        x = x'; y = y'; color = ('red');
    case 3 % Use Green Signal
        [x,y,~,~] = extractPosition2(dTargets,3);
        x = x'; y = y'; color = ('green');
    case 0 % Use luminance
        [x,y,~,~] = extractPosition2(dTargets,1);
        x = x'; y = y'; color = ('lum');  
   % case 4 is directly using x y 
   
    case 5 %use points
        % Convert data
        points = uint32(points);      
        [x, y] = Points2XY(points);
end

% If videoFix screwed up timestamps, (makes 10^10 instead 10^4)
if t(1) < 0
errordlg(sprintf('Negative timestamps in: %s', file));
end
if t(1) > 10^14
   errordlg('Fixing Video timestamps')
   t = t * 10^-6;
end   

% Suppress postions at [0,0] because this is caused by bad tracking
%[x,y,t] = suppressZeros(x,y,t);
    
ind = x==0;
%t(ind) = NaN; %[];
x(ind) = NaN; %[];
y(ind) = NaN; %[];

% plot original track
plot(x,y)
title('None smoothed path')

% smooth position
[x,y] = meanpath(x,y,p);
          
% Add angle
angle(1:15) = NaN;
for i = 16:(length(x)-15)
    %disp(i);
    if (x(i+15) == x(i-15)) && (y(i+15) == y(i-15))
        angle = NaN;
    else
        [a, d] = cart2pol(x(i+15)-x(i-15),y(i+15)-y(i-15));
        angle(i) = a;
    end
end 
angle((length(x)-14):length(x)) = NaN;


function [x,y] = meanpath(x,y,p)

temp_x = x;
temp_y = y;
for cc=p.smoothSize+1:length(x)-p.smoothSize
  x_window = x(cc-p.smoothSize:cc+p.smoothSize); 
  y_window = y(cc-p.smoothSize:cc+p.smoothSize);
  temp_x(cc) = nanmean(x_window); 
  temp_y(cc) = nanmean(y_window);
end

x = temp_x;
y = temp_y;
