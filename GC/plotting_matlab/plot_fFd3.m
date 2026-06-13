function plot_fFd3(t, f_list, Fd3_list, initial, bounds, crop, fig_dir)

fig = figure;
colors_setting = ColorsData();
colors = colors_setting.colors;
red = colors_setting.red; % for bounds

if initial.init_n > 20
   colors =  jet(initial.init_n+1);
end

crop_index = floor(size(f_list, 2) / crop);
tiledlayout(2,1,'TileSpacing','tight','Padding','tight')

%% plot1
% nexttile([1 2])
nexttile

hold on;
yline(bounds.F_bound, 'Color', red, 'LineWidth', 1.5);
for i = 1:initial.init_n
    fi = abs(squeeze(f_list(i,:)));
    plot(t, fi, 'Color', colors(i, :), 'LineWidth', 1.5);
end

xlim([0 t(length(t))]);
ylim([0.95*min(f_list,[],"all") 1.01*bounds.F_bound]);
% Add labels and title
xlabel('$t$ [s]','interpreter','latex');
ylabel('$|f(t)|$ [N]','interpreter','latex');
% title('Norm of position error vs time','interpreter','latex');

legend('$\overline{\mathcal{F}}$','interpreter','latex');

% Adjust plot appearance
box on;
set(gca,'ticklabelinterpreter','latex');
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.2);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02, 0.02]);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'ZMinorTick', 'on');
set(gca,'xtick',[]);
hold off;


%% plot2
% nexttile([1 2])
nexttile

hold on;
for i = 1:initial.init_n
    fi = squeeze(Fd3_list(i,:));
    % plot(t(1:200), error_p(1:200), 'Color', colors(i, :), 'LineWidth', 1.5);
    plot(t, fi, 'Color', colors(i, :), 'LineWidth', 1.5);
end

% Add labels and title
xlabel('$t$ [s]','interpreter','latex');
ylabel('$F_{d,3}(t)$ [N]','interpreter','latex');
% title('Norm of position error vs time','interpreter','latex');

xlim([0 t(length(t))]);

% Adjust plot appearance
box on;
% set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.2);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02, 0.02]);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca,'TickLabelInterpreter','latex');
% set(gca, 'YScale', 'log');
hold off;

set(fig,'renderer','Painters')
exportgraphics(fig, fullfile(fig_dir, 'Fig8_F.png'), 'Resolution', 600);
print(fig, fullfile(fig_dir, 'Fig8_F.eps'), '-depsc');

end