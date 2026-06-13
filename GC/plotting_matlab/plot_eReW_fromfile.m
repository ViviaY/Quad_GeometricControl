function plot_eReW_fromfile(matfile)

    % === 1. Load data ===
    data = load(matfile, "t_span", "eR_list", "eW_list");
    t = data.t_span;
    eR_list = data.eR_list;
    eW_list = data.eW_list;

    % === 2. Create figure with 2 subplots ===
    figure;

    % ----- Subplot 1: ||e_R(t)|| -----
    subplot(1,2,1);
    plot(t, eR_list', 'LineWidth', 1.0);
    xlabel("Time (s)", 'Interpreter','latex');
    ylabel("$\|e_R(t)\|$", 'Interpreter','latex');
    title("Rotation Error Norm", 'Interpreter','latex');
    grid on;

    % ----- Subplot 2: ||e_\omega(t)|| -----
    subplot(1,2,2);
    plot(t, eW_list', 'LineWidth', 1.0);
    xlabel("Time (s)", 'Interpreter','latex');
    ylabel("$\|e_\omega(t)\|$", 'Interpreter','latex');
    title("Angular Velocity Error Norm", 'Interpreter','latex');
    grid on;

    sgtitle('Error Norm Trajectories (100 runs)', 'Interpreter','latex');

    % === 3. Compute initial errors ===
    eR_init = eR_list(:, 1);
    eW_init = eW_list(:, 1);

    % === 4. Print statistics ===
    fprintf("Initial ||e_R(0)||: min = %.4f, mean = %.4f, max = %.4f\n", ...
        min(eR_init), mean(eR_init), max(eR_init));

    fprintf("Initial ||e_\\omega(0)||: min = %.4f, mean = %.4f, max = %.4f\n", ...
        min(eW_init), mean(eW_init), max(eW_init));

end
