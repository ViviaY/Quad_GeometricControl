function [p0_multi, v0_multi, R0_multi, w0_multi] = CBF_InitialPoints(initial_points_file)
    % Load initial points data or generate random initial points
    try
        init_data = load(initial_points_file);
        fprintf('loading the initial points: %s\n', initial_points_file);
    catch ME
        fprintf('failed to load initial points: %s\n', ME.message);
        p0_multi = []; v0_multi = []; R0_multi = []; w0_multi = [];
        return;
    end

    % Access the struct fields
    initial = init_data.initial;
    initial_points_data = initial.initial_points;

    fprintf('initial points data size: %s\n', mat2str(size(initial_points_data)));

    num_trajectories = size(initial_points_data, 2);

    % Extract positions, velocities, and angular velocities
    p0_multi = initial_points_data(1:3, :).';   % transpose to match Python .T
    v0_multi = initial_points_data(4:6, :).';
    w0_multi = initial_points_data(7:9, :).';

    % Rotation matrices
    R0_data = initial_points_data(10:18, :);
    R0_multi = zeros(num_trajectories, 3, 3);

    for i = 1:num_trajectories
        R_flat = R0_data(:, i);
        R0_multi(i, :, :) = reshape(R_flat, [3, 3]); 
    end
end
