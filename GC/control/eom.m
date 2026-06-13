function Xdot = eom(t, X, annealing_output, param, Points_Array, tau)

e3 = [0, 0, 1]';
m = param.m;
J = param.J;


[~, v, W, R] = split_to_states(X);

desired = DesiredTrajectory(t, Points_Array, tau, param);
[f, M, ~, ~] = position_control(X, desired, annealing_output, param);

xdot = v;
vdot = - param.g * e3 + f / m * R * e3;
Wdot = J \ (-hat(W) * J * W + M);
Rdot = R * hat(W);

Xdot=[xdot; vdot; Wdot; reshape(Rdot,9,1)];

end