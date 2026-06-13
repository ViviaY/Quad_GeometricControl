function [phi, theta, psi] = euler_from_R_zyx(R)
% Extract ZYX Euler angles (phi, theta, psi) from rotation matrix R
% R = Rz(psi) @ Ry(theta) @ Rx(phi)
% Returns (phi, theta, psi) as floats in radians

% Ensure R is a matrix
R = double(R);

% Protect against numerical issues
r20 = R(3,1);
theta = atan2(-r20, sqrt(R(3,2)^2 + R(3,3)^2));
phi = atan2(R(3,2), R(3,3));
psi = atan2(R(2,1), R(1,1));

end