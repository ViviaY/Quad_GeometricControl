function plot_eReW(t, eR_list, ew_list, initial, crop, fig_dir)

colors_setting = ColorsData();
colors = colors_setting.colors;
red = colors_setting.red; % only used for bounds

crop_index = floor(size(eR_list, 2) / crop);

% tiledlayout(3,1,'TileSpacing','tight','Padding','tight')
tiledlayout(2,1,'TileSpacing','tight')

%% plot1
nexttile
hold on;
if initial.init_n > 20
   colors =  jet(initial.init_n+1);
end
% plotting errors e_R(t)
for i = 1:initial.init_n
    error_R = squeeze(eR_list(i,1:crop_index));
    plot(t(1:crop_index), error_R(1:crop_index), 'Color', colors(i, :), 'LineWidth', 1.5);
end

xlim([0 t(crop_index)]+0.1);

% Add labels and title
% xlabel('$t$ [s]','interpreter','latex');
ylabel('$\|e_R(t)\|$ [rad]','interpreter','latex');
% title('Norm of position error','interpreter','latex');
% legend('$\mathcal{L}_{p}(\overline{\mathcal{V}}_{1},\overline{\mathcal{V}}_{2}, t)$','interpreter','latex');
% Adjust plot appearance
box on;
% set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.2);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02, 0.02]);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca,'xtick',[]);
set(gca,'TickLabelInterpreter','latex');
% set(gca, 'YScale', 'log');
hold off;

%% plot2
nexttile
hold on;

% plotting errors e_v(t)
for i = 1:initial.init_n
    eW = squeeze(ew_list(i,:));
    plot(t(1:crop_index), eW(1:crop_index), 'Color', colors(i, :), 'LineWidth', 1.5);
end

xlim([0 t(crop_index)]+0.1);

% Add labels and title
% xlabel('$t$ [s]','interpreter','latex');
ylabel('$\|e_w(t)\|$ [rad/s]','interpreter','latex');
% title('Norm of position error','interpreter','latex');
% legend('$\mathcal{L}_{v}(\overline{\mathcal{V}}_{1},\overline{\mathcal{V}}_{2}, t)$','interpreter','latex');
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

% % Save PNG
exportgraphics(gcf, fullfile(fig_dir, 'Fig_eR_eW.png'), 'Resolution', 600);
% Save EPS
print(gcf, fullfile(fig_dir, 'Fig_eR_eW.eps'), '-depsc');

end