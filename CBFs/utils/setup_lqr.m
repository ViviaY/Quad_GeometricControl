function K = setup_lqr()

% LQR state and control weights
Q = diag([10, 10, 10, 1, 1, 1]);  % position and velocity weights
R = eye(3);                       % control input weights

% Continuous-time system matrices
A = [zeros(3,3), eye(3);
     zeros(3,6)];
B = [zeros(3,3);
     eye(3)];
     

[K, ~, ~] = lqr(A, B, Q, R);

end