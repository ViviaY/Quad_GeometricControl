function plot_error_norms(t, error_data, Lpt, Lvt, init_n, crop)
% Function to plot the norm of position, velocity, and attitude errors
%
% Inputs:
%   error_data - Cell array containing error structures for each trajectory
%   init_n - Number of trajectories to plot
%
% Example usage:
%   plot_error_norms(error_data, 20)
       

    colors = [0 0.4470 0.7410;...       % Blue (Default MATLAB Blue)
              0.8500 0.3250 0.0980;...  % Orange (MATLAB Default)
              0.9290 0.6940 0.1250;...  % Yellow
              0.4940 0.1840 0.5560;...  % Purple
              0.4660 0.6740 0.1880;...  % Green
              0.3010 0.7450 0.9330;...  % Light Blue
              0 0 0;...                 % Black
              0.75 0.75 0.75;...        % Gray
              1 0 1;...                 % Magenta
              0 1 1;...                 % Cyan
              0 0.5 0;...               % Dark Green (Olive)
              0 0.4470 0.7410;...       % Deep Blue (Repeated to ensure visibility)
              0.5 0 0.5;...             % Dark Purple
              0.3 0.3 0.3;...           % Dark Gray
              0.6350 0.0780 0.1840;...  % Dark Brown
              0.2 0.5 0.3;...           % Forest Green
              0.2 0.2 0.6;...           % Dark Navy Blue
              0.3 0.6 0.8;...           % Sky Blue
              0.8 0.4 0.0;...           % Burnt Orange
              0.8 0.6 0.7];             % Soft Pink
    
   red = [1,0,0]; % only used for bounds

   tiledlayout(3,1,'TileSpacing','tight','Padding','tight')
   LineWidth = 2;
   title_fontsize = 12;
   lable_fontsize = 14;
   

    
    for j = 1:init_n
        t_j = error_data{j}.t;   % Time vector
        p_j = error_data{j}.p;   % Position error (N×3)
        v_j = error_data{j}.v;   % Velocity error (N×3)

        % Compute the norm (magnitude) of each error
        norm_p = vecnorm(p_j, 2, 2);  % Compute the 2-norm for each row of p_j
        norm_v = vecnorm(v_j, 2, 2);  % Compute the 2-norm for each row of v_j

        subplot(3,1,1);
        hold on;
        yline(bounds.Lu, 'Color', [0.8, 0.2, 0.4], 'LineWidth', 1.5);
        
        for i = 1:size(lyap_V, 1)
            plot(t(1:crop_index), lyap_V(i,1:crop_index), 'Color', colors(i, :), 'LineWidth', 1.5);
        end
        % title("Lyapunov");
        % plot(V1);
        % plot(V2);
        % plot(V_bound);
        % yline(uniform_V_bound);
        legend("$\mathcal{L}_{u}^2(\overline{\mathcal{V}}_{1},\overline{\mathcal{V}}_{2})$",'interpreter','latex');

        % Plot the norm of position error
        subplot(3,1,2);
        hold on;
        plot(t_j, norm_p, 'Color', colors(j, :), 'LineWidth', LineWidth);
        title('Position Errors', 'FontSize',title_fontsize);
        % xlabel('$t$ [s]','interpreter','latex', 'FontSize',lable_fontsize);
        ylabel('$\|e_p\|$ [m]','interpreter','latex', 'FontSize',lable_fontsize);

        % Plot the norm of velocity error
        subplot(3,1,3);
        hold on;
        plot(t_j, norm_v, 'Color', colors(j, :), 'LineWidth', LineWidth);
        
        title('Velocity Errors', 'FontSize',title_fontsize);
        % xlabel('$t$ [s]','interpreter','latex', 'FontSize',lable_fontsize);
        ylabel('$\|e_v\|$ [m/s]','interpreter','latex', 'FontSize',lable_fontsize);

        % Plot the norm of attitude error

        plot(t_j, norm_R, 'Color', colors(j, :), 'LineWidth', LineWidth);
        title('Attitude Errors', 'FontSize', title_fontsize);
        xlabel('$t$ [s]','interpreter','latex', 'FontSize',lable_fontsize);
        ylabel('$\|e_R\|$','interpreter','latex', 'FontSize',lable_fontsize);
    end


    % **Add Bounds curve after the loop**
    subplot(3,1,2);
    hold on;
    plot_Lpt = plot(t, Lpt, 'Color', [1, 0, 0], 'LineWidth', LineWidth, ...
                  'DisplayName', '$\mathcal{L}_{p}(\overline{\mathcal{V}}_{1},\overline{\mathcal{V}}_{2})$');  % Red color

    legend(plot_Lpt, 'Interpreter', 'latex', 'Location', 'best'); 
    
    subplot(3,1,3);
    plot_Lvt = plot(t, Lvt, 'Color', [1, 0, 0], 'LineWidth', LineWidth, ...
                  'DisplayName', '$\mathcal{L}_{v}(\overline{\mathcal{V}}_{1},\overline{\mathcal{V}}_{2})$');  % Red color
    legend(plot_Lvt, 'Interpreter', 'latex', 'Location', 'best'); 





    % Add legends to distinguish different trajectories
    subplot(3,1,1); grid on; %legend('show');
    subplot(3,1,2); grid on; %legend('show');
    subplot(3,1,3); grid on; %legend('show');

end
