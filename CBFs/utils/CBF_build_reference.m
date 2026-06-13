function [p_d, v_d, a_d, R_d, f_d, w_d, tau_d, yaw_d] = CBF_build_reference(bezier_samples, quadparam, track_tspan)
% Construct reference trajectory from Bezier sampled points
%
% Inputs:
%   bezier_samples: struct with fields
%       - position: [N x 3]
%       - velocity: [N x 3]
%       - acceleration: [N x 3]
%   quadparam: struct with fields
%       - m: mass
%       - g: gravity constant
%       - J: inertia matrix [3x3]
%   track_tspan: [N x 1] time vector
%
% Outputs:
%   p_d: [N x 3] desired positions
%   v_d: [N x 3] desired velocities
%   a_d: [N x 3] desired accelerations
%   R_d: {N x 1} desired rotation matrices (3x3 each)
%   f_d: [N x 1] desired thrusts
%   w_d: [N x 3] desired angular velocities
%   tau_d: [N x 3] desired torques
%   yaw_d: [N x 1] desired yaw angles

    % Extract Bezier sampled data
    p_d = bezier_samples.position;
    v_d = bezier_samples.velocity;
    a_d = bezier_samples.acceleration;

    num_samples = size(p_d, 1);

    % Initialize results
    R_d = cell(num_samples,1);
    f_d = zeros(num_samples,1);
    w_d = zeros(num_samples,3);
    tau_d = zeros(num_samples,3);
    yaw_d = zeros(num_samples,1);

    % Loop through each time step
    for i = 1:num_samples
        acc = a_d(i,:);
        gravity = [0, 0, quadparam.g];
        
        % Desired z-axis of body frame
        acc_with_g = acc + gravity;
        norm_acc = norm(acc_with_g);
        if norm_acc < 1e-6
            z_B = [0, 0, 1];
        else
            z_B = acc_with_g / norm_acc;
        end
        
        % Construct desired rotation matrix
        y_C = [-sin(yaw_d(i)), cos(yaw_d(i)), 0];
        x_B = cross(y_C, z_B);
        x_B_norm = norm(x_B);
        if x_B_norm < 1e-6
            y_C = [0, 1, 0];
            x_B = cross(y_C, z_B);
            x_B_norm = norm(x_B);
        end
        x_B = x_B / x_B_norm;
        y_B = cross(z_B, x_B);
        R_di = [x_B(:), y_B(:), z_B(:)];
        R_d{i} = R_di;

        % Desired thrust
        f_d(i) = quadparam.m * dot(acc_with_g, z_B);

        % Angular velocity and torque
        if i == 1
            w_d(i,:) = [0, 0, 0];
            tau_d(i,:) = [0, 0, 0];
        else
            dt = track_tspan(i) - track_tspan(i-1);
            R_prev = R_d{i-1};
            R_curr = R_di;

            % Approximate R_dot
            R_dot = (R_curr - R_prev) / dt;
            omega_hat = R_prev' * R_dot;

            % Extract angular velocity
            omega = [omega_hat(3,2), omega_hat(1,3), omega_hat(2,1)];
            w_d(i,:) = omega;

            % Angular acceleration
            if i == 2
                omega_dot = omega / dt;
            else
                omega_dot = (omega - w_d(i-1,:)) / dt;
            end

            % Torque: J*omega_dot + omega x (J*omega)
            J_omega = (quadparam.J * omega.').';
            tau_di = (quadparam.J * omega_dot.')' + cross(omega, J_omega);
            tau_d(i,:) = tau_di;
        end
    end

    % Ensure first element is consistent
    if num_samples > 1
        w_d(1,:) = w_d(2,:);
        tau_d(1,:) = tau_d(2,:);
    end
end
