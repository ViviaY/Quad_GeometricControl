classdef ColorsData < handle
    properties
        colors = [0 0.4470 0.7410;...       % Blue (Default MATLAB Blue)
              0.8500 0.3250 0.0980;...  % Orange (MATLAB Default)
              0.9290 0.6940 0.1250;...  % Yellow
              0.4940 0.1840 0.5560;...  % Purple
              0.4660 0.6740 0.1880;...  % Green
              0.3010 0.7450 0.9330;...  % Light Blue
              0 0 0;...                 % Black
              0.75 0.75 0.75;...        % Gray
              1 0 1;...                 % Magenta
              0 1 1;...                 % Cyan
              0 0.5 0;...               % Dark Green (Olive)
              0 0.4470 0.7410;...       % Deep Blue (Repeated to ensure visibility)
              0.5 0 0.5;...             % Dark Purple
              0.3 0.3 0.3;...           % Dark Gray
              0.6350 0.0780 0.1840;...  % Dark Brown
              0.2 0.5 0.3;...           % Forest Green
              0.2 0.2 0.6;...           % Dark Navy Blue
              0.3 0.6 0.8;...           % Sky Blue
              0.8 0.4 0.0;...           % Burnt Orange
              0.29,0.0,0.51;...         % Indigo
              0.8 0.6 0.7];             % Soft Pink
        red =[1,0,0];                   % used for bounds                   
    end
end