function mu = nominal_lqr_controller(z, x_r, K)
% Generate nominal control input using LQR controller
% Inputs:
%   z - current state [6x1] (position, velocity)
%   x_r - reference state [9x1] (position, Euler angles, velocity)
%   K - LQR gain matrix
% Output:
%   mu - nominal control input [3x1] (virtual acceleration)

% Extract position and velocity errors
pos_err = z(1:3) - x_r(1:3);
vel_err = z(4:6) - x_r(7:9);
state_err = [pos_err; vel_err];

% Apply LQR control law
mu = -K * state_err;

end