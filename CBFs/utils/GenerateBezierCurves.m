function sampled_trajectory = GenerateBezierCurves(bezier_data, track_tspan)
    % Sample uniformly spaced trajectory points from MATLAB Bezier data
    %
    % Inputs:
    %   bezier_data: struct with fields:
    %       - Points_Array: [3 x N_p x N_seg]
    %       - t_Array: [1 x (N_seg+1)] segment boundary times
    %   track_tspan: vector of sample times
    %
    % Outputs:
    %   sampled_trajectory: struct with position, velocity, acceleration

    % Original layout: [3 x N_p x N_seg]
    P = bezier_data.Points_Array;
    t_Array = bezier_data.t_Array;

    % Reorder to [N_seg x N_p x 3]
    control_points = permute(P, [3,2,1]);  

    fprintf('Control points size (after permute): %s\n', mat2str(size(control_points)));
    fprintf('t_Array size: %s\n', mat2str(size(t_Array)));

    num_curves = size(control_points, 1);
    if length(t_Array) ~= num_curves+1
        error('Number of time stamps (%d) does not match number of curve segments (%d)', ...
              length(t_Array), num_curves);
    end
    fprintf('Number of Bezier segments: %d\n', num_curves);
    fprintf('Total trajectory time: %.3f\n', t_Array(end));

    num_samples = length(track_tspan);

    % Preallocate arrays
    position_samples = zeros(num_samples, 3);
    velocity_samples = zeros(num_samples, 3);
    acceleration_samples = zeros(num_samples, 3);

    % Loop through each sample time
    for i = 1:num_samples
        t = track_tspan(i);

        % Segment index
        segment_idx = find(t_Array <= t, 1, 'last');
        if isempty(segment_idx)
            segment_idx = 1;
        elseif segment_idx >= num_curves+1
            segment_idx = num_curves;
        end

        % Local normalized time
        segment_start_time = t_Array(segment_idx);
        segment_duration = t_Array(segment_idx+1) - segment_start_time;
        tau = (t - segment_start_time) / segment_duration;
        tau = max(0, min(tau, 1.0));

        % Extract control points of current segment: [N_p x 3]
        segment_control_points = squeeze(control_points(segment_idx,:,:));

        % Evaluate Bezier
        [pos, vel, acc] = evaluate_bezier(segment_control_points, tau, segment_duration);
        position_samples(i,:) = pos;
        velocity_samples(i,:) = vel;
        acceleration_samples(i,:) = acc;
    end

    % Pack results
    sampled_trajectory.position = position_samples;
    sampled_trajectory.velocity = velocity_samples;
    sampled_trajectory.acceleration = acceleration_samples;
    sampled_trajectory.t_Array = t_Array;
    sampled_trajectory.original_control_points = control_points;
end


%% ===== Subfunction: Evaluate Bezier curve =====
function [position, velocity, acceleration] = evaluate_bezier(control_points, tau, duration)
    % Evaluate Bezier curve position, velocity, and acceleration
    % control_points: [N_p x 3] (row = control point, col = x,y,z)
    % tau: normalized time in [0,1]
    % duration: segment duration

    n = size(control_points,1) - 1; % degree of the curve

    % Position
    position = zeros(1,3);
    for j = 0:n
        coeff = nchoosek(n,j) * (1-tau)^(n-j) * tau^j;
        position = position + coeff * control_points(j+1,:);
    end

    % Velocity
    velocity = zeros(1,3);
    for j = 0:n-1
        coeff = n * (nchoosek(n-1,j) * (1-tau)^(n-1-j) * tau^j);
        diff_pt = control_points(j+2,:) - control_points(j+1,:);
        velocity = velocity + coeff * diff_pt;
    end
    velocity = velocity / duration;

    % Acceleration
    acceleration = zeros(1,3);
    for j = 0:n-2
        coeff = n*(n-1) * (nchoosek(n-2,j) * (1-tau)^(n-2-j) * tau^j);
        diff2_pt = control_points(j+3,:) - 2*control_points(j+2,:) + control_points(j+1,:);
        acceleration = acceleration + coeff * diff2_pt;
    end
    acceleration = acceleration / (duration^2);
end
