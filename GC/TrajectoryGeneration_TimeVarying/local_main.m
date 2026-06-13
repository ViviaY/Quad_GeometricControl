%% clear all caches
clc;           % Clear the command window
clear;         % Clear all variables
clear all;     % Clear all variables, functions, and class definitions (more thorough)
close all;     % Close all open figure windows
fclose('all'); % Close all open files
rehash;        % Refresh MATLAB's search path and function cache


%% common paramters
e3=[0;0;1];

% inertia matrix
%J=diag([1.43e-5,1.43e-5,2.89e-5]); % crazyfly
param.J = diag([0.0820, 0.0845, 0.1377]);% CDC
%J=diag([0.43,0.43, 1.02])*0.01;
%J=diag([0.007, 0.007, 0.012]);
%J=diag([0.03, 0.03, 0.03]);
%J=diag([0.0196,0.0196,0.0264]);
param.J_min = min(diag(param.J));
param.J_max = max(diag(param.J)); 


param.g = 9.81;
%mass
%m=0.033; % crazyfly
param.m = 4.34; % CDC
%m=0.755;
%m=1.4;
%m=0.8;

% acceleration bound

% CrazyFly gain parameters and parameters for the theoretical bound
%  vm=[5;5;5];
%  am=[0.1;0.1;10];
%  kp=2;
%  kv=0.1;
%  kR=2;
%  kw=0.1;
%  V1_0=0.05; % V1(0)
%  psi_bar=0.02;
%  alpha_psi=0.9;
% c1=0.5*min([sqrt(kp*m), 4*m*kp*kv/(kv^2+4*m*kp)]);
% c2=0.3*min([sqrt(kR*lambda_m_J), 4*(lambda_m_J)*kR*kw/((kw^2)+4*(lambda_m_J)*kR)]);


% % CDCPaper parameters
anneal_options.vm = [2;2;2];
anneal_options.am = [1;1;10];
anneal_options.psi_bar = 0.04;
anneal_options.alpha_psi = 0.4;
anneal_options.V1_0 = 1; 

annealing_output.opt_k = [16.5065, 5.3233, 24.5341, 1.4898, 0.5494, 0.6409];
% if you want to skip ObtainingBoundsTimeVarying and already have
% TimeVaryingBoundsData.mat, you can comment ObtainingBoundsTimeVarying and
% uncomment the following lines
% Bounds = load("TimeVaryingBoundsData.mat")
plotting = false;
Bounds = ObtainingBoundsTimeVarying(annealing_output, anneal_options, param, plotting); 

generatedTraj = Trajectory_Synthesis_TimeVarying(param, anneal_options, Bounds);
% Same as before, if you want to skip Trajectory_Synthesis_TimeVarying and already have
% TimeVaryingBoundsData.mat, you can uncomment the following line
% generatedTraj = load("GeneratedTrajectoryData.mat");

GeneratingPlotsAndOutPuts_TimeVarying(param, anneal_options, generatedTraj, Bounds);



