function z_next = rk4_step(z, acc, dt)
% RK4 integration for translational dynamics
% Inputs:
%   z - current state [6x1] (position and velocity)
%   acc - acceleration [3x1]
%   dt - time step [scalar]
% Output:
%   z_next - next state [6x1]

% Define state derivative function
function dz = dz_dt(z_curr, acc_curr)
    dz = zeros(size(z_curr));
    dz(1:3) = z_curr(4:6);  % position derivative = velocity
    dz(4:6) = acc_curr;     % velocity derivative = acceleration
end

% RK4 steps
k1 = dz_dt(z, acc);
k2 = dz_dt(z + 0.5*dt*k1, acc);
k3 = dz_dt(z + 0.5*dt*k2, acc);
k4 = dz_dt(z + dt*k3, acc);

% Update state
z_next = z + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4);

end