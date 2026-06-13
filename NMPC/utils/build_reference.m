function x_ref = build_reference(bezier_samples, track_tspan, quadparam)
    % Construct reference states for a quadrotor from Bezier samples.
    % Extract samples
    p_d = bezier_samples.position;
    v_d = bezier_samples.velocity;
    a_d = bezier_samples.acceleration;
    v_thresh = 1e-3;

    N = size(p_d, 1);
    if size(v_d,1) ~= N || size(a_d,1) ~= N
        error('position/velocity/acceleration must have the same number of rows.');
    end

    if numel(track_tspan) ~= N
        error('track_tspan length must match number of samples N.');
    end

    % Allocate outputs
    R_d   = zeros(3,3,N);
    xi_d  = zeros(N,3);
    f_d   = zeros(N,1);
    w_d   = zeros(N,3);
    tau_d = zeros(N,3);
    yaw0 = 0;  
    yaw_d = compute_yaw_from_velocity(v_d, a_d, yaw0, v_thresh);

    m = quadparam.m;
    g = quadparam.g;
    J = quadparam.J;

    % Main loop
    for i = 1:N
        % Desired body z-axis from acceleration + gravity
        gravity    = [0; 0; g];
        acc = a_d(i, :).' + gravity;
        norm_acc = norm(acc);
        if norm_acc < 1e-9
            z_B = [0;0;1];
        else
            z_B = acc / norm_acc;
        end

        % Construct x_B using yaw=0 frame's y_C
        % y_C for yaw = 0 is [0; 1; 0]
        y_C = [-sin(yaw_d(i)); cos(yaw_d(i)); 0]; 
        x_B = cross(y_C, z_B);
        if norm(x_B) < 1e-9
            y_C = [0; 1; 0];
            x_B = cross(y_C, z_B);
        end
        x_B = x_B / max(norm(x_B), 1e-12);
        y_B = cross(z_B, x_B);

        % Assemble rotation matrix and (optionally) re-orthonormalize
        R = [x_B, y_B, z_B];
        R = orthonormalize(R);
        R_d(:,:,i) = R;
        xi_d(i, :) =  so3_log(R);

        % Desired thrust along z_B direction
        f_d(i) = m * dot(acc, z_B);

        % Angular velocity and torque by finite differences (skip i=1)
        if i == 1
            w_d(i,:)   = [0,0,0];
            tau_d(i,:) = [0,0,0];
        else
            dt = max(track_tspan(i) - track_tspan(i-1), 1e-12);

            % Time derivative of rotation (simple finite difference)
            R_prev = R_d(:,:,i-1);
            R_curr = R_d(:,:,i);
            R_dot  = (R_curr - R_prev) / dt;

            % omega_hat in body frame: R_prev' * R_dot
            omega_hat = R_prev.' * R_dot;

            % Use skew-symmetric part to be robust to numerical noise
            omega_hat = 0.5 * (omega_hat - omega_hat.');

            % Extract angular velocity vector (vee operator)
            omega = vee3(omega_hat).'; % 1x3
            w_d(i,:) = omega;

            % Angular acceleration by finite difference
            if i == 2
                omega_dot = (w_d(i,:) - [0,0,0]) / dt;
            else
                omega_dot = (w_d(i,:) - w_d(i-1,:)) / dt;
            end

            % Torque: J*omega_dot + omega x (J*omega)
            Jomega   = (J * omega.').';
            tau_d(i,:) = (J * omega_dot.').'+ cross(omega, Jomega);
        end
    end

    % Backfill the first sample for w_d and tau_d if N>1
    if N > 1
        w_d(1,:)   = w_d(2,:);
        tau_d(1,:) = tau_d(2,:);
    end

    x_ref = [p_d, v_d, xi_d, w_d];
end


function yaw_d = compute_yaw_from_velocity(v_d, a_d, yaw0, v_thresh)
    % Compute yaw from horizontal velocity with robust fallbacks and unwrap.
    if nargin < 3 || isempty(yaw0),    yaw0 = 0;      end
    if nargin < 4 || isempty(v_thresh), v_thresh = 1e-3; end

    % Ensure sizes
    N = size(v_d,1);
    v_d = double(v_d);  % in case inputs are casadi DM
    if nargin < 2 || isempty(a_d)
        a_d = zeros(N,3);
    else
        a_d = double(a_d);
        if size(a_d,1) ~= N
            error('a_d must have the same number of rows as v_d.');
        end
    end

    % Normalize yaw0 to [-pi, pi] without Mapping Toolbox
    yaw0 = atan2(sin(yaw0), cos(yaw0));

    yaw_d      = zeros(N,1);
    yaw_d(1)   = yaw0;

    for i = 1:N
        vx = v_d(i,1); vy = v_d(i,2);
        if ~isfinite(vx) || ~isfinite(vy), vx = 0; vy = 0; end

        sp = hypot(vx, vy);
        if sp >= v_thresh
            psi = atan2(vy, vx);
        else
            ax = a_d(i,1); ay = a_d(i,2);
            if ~isfinite(ax) || ~isfinite(ay), ax = 0; ay = 0; end
            asp = hypot(ax, ay);
            if asp >= v_thresh
                psi = atan2(ay, ax);
            else
                % Hold previous yaw
                if i == 1
                    psi = yaw0;
                else
                    psi = yaw_d(i-1);
                end
            end
        end

        % Write and unwrap vs previous sample
        yaw_d(i) = psi;
        if i > 1
            d = yaw_d(i) - yaw_d(i-1);
            if d >  pi, yaw_d(i) = yaw_d(i) - 2*pi; end
            if d < -pi, yaw_d(i) = yaw_d(i) + 2*pi; end
        end
    end
end



