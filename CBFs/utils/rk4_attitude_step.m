function [R_next, omega_next] = rk4_attitude_step(R, omega, tau, J, dt)
% RK4 integration for rotational dynamics
% Inputs:
%   R - current rotation matrix [3x3]
%   omega - current angular velocity [3x1]
%   tau - external torque [3x1]
%   J - inertia matrix [3x3]
%   dt - time step [scalar]
% Outputs:
%   R_next - next rotation matrix [3x3]
%   omega_next - next angular velocity [3x1]

% Define attitude derivative function
function dR = dR_dt(R_curr, omega_curr)
    dR = R_curr * hat_map(omega_curr);
end

% Define angular velocity derivative function
function domega = domega_dt(omega_curr, tau_curr)
    domega = J \ (tau_curr - cross(omega_curr, J * omega_curr));
end

% RK4 for rotation matrix R
k1_R = dR_dt(R, omega);
k2_R = dR_dt(R + 0.5 * dt * k1_R, omega);
k3_R = dR_dt(R + 0.5 * dt * k2_R, omega);
k4_R = dR_dt(R + dt * k3_R, omega);
R_next = R + (dt / 6.0) * (k1_R + 2 * k2_R + 2 * k3_R + k4_R);

% Project back to SO(3) via polar decomposition
[U, ~, V] = svd(R_next);
R_next = U * V';

% RK4 for angular velocity omega
k1_w = domega_dt(omega, tau);
k2_w = domega_dt(omega + 0.5 * dt * k1_w, tau);
k3_w = domega_dt(omega + 0.5 * dt * k2_w, tau);
k4_w = domega_dt(omega + dt * k3_w, tau);
omega_next = omega + (dt / 6.0) * (k1_w + 2 * k2_w + 2 * k3_w + k4_w);

end

function S = hat_map(v)
% Skew-symmetric matrix from vector (so(3))
S = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end