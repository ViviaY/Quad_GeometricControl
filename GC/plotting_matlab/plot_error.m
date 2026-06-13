function plot_error(t, ep_list, ev_list, lyap_V, Lpt, Lvt, L, ...
                    initial, crop, fig_dir)


base_blue = [0.5, 0.7, 1.0]; 
noise_level = 0.5;          

colors = base_blue + noise_level * randn(100,3);
colors = max(min(colors,1),0); % clamp to [0,1]
red = [1, 0, 0];
crop_index = floor(size(ep_list, 2) / crop);
tiledlayout(3,1,'TileSpacing','tight')

%% plot1
nexttile
hold on;

% plotting the bound of Lp(t)
plot(t(1:crop_index), Lpt(1:crop_index), 'Color', red, 'LineWidth', 3.5);

% plotting errors e_p(t)
for i = 1:initial.init_n
    error_p = squeeze(ep_list(i,1:crop_index));
    plot(t(1:crop_index), error_p(1:crop_index), 'Color', colors(i, :), 'LineWidth', 1.5);
end

xlim([0 t(crop_index)]+0.1);
ylim([0 max(error_p(1:crop_index))+0.35]);

% Add labels and title
% xlabel('$t$ [s]','interpreter','latex');
ylabel('$\|e_p(t)\|$ [m]','interpreter','latex');
% title('Norm of position error','interpreter','latex');
l1 =legend('$\mathcal{L}_{p}(\overline{\mathcal{V}}_{1},\overline{\mathcal{V}}_{2}, t)$','interpreter','latex');
% l1 = legend('show', 'Location', 'northeast');
% Get the current position of the legend
pos1 = l1.Position;
% Adjust the position, for example, moving the legend downward
l1.Position = [pos1(1) - 0.02, pos1(2) - 0.01, pos1(3), pos1(4)];

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
set(gca, 'YAxisLocation', 'left','XTickLabel',[], 'box', 'off');
set(gca,'TickLabelInterpreter','latex');
% set(gca, 'YScale', 'log');

hold off;

%% plot2
nexttile
hold on;

% plotting the bound of Lv(t)
plot(t(1:crop_index), Lvt(1:crop_index), 'Color', red, 'LineWidth', 3.5);

% plotting errors e_v(t)
for i = 1:initial.init_n
    vi = squeeze(ev_list(i,:));
    plot(t(1:crop_index), vi(1:crop_index), 'Color', colors(i, :), 'LineWidth', 1.5);
end

xlim([0 t(crop_index)]+0.1);
ylim([0 max(Lvt(1:crop_index))+0.1]);

% Add labels and title
% xlabel('$t$ [s]','interpreter','latex');
ylabel('$\|e_v(t)\|$ [m/s]','interpreter','latex');
% title('Norm of position error','interpreter','latex');
l2 = legend('$\mathcal{L}_{v}(\overline{\mathcal{V}}_{1},\overline{\mathcal{V}}_{2}, t)$','interpreter','latex');
% Adjust plot appearance
% l2 = legend('show', 'Location', 'northeast');
% Get the current position of the legend
pos2 = l2.Position;
% Adjust the position, for example, moving the legend downward
l2.Position = [pos2(1) - 0.02, pos2(2) - 0.01, pos2(3), pos2(4)];

box on;
set(gca,'ticklabelinterpreter','latex');
set(gca, 'YAxisLocation', 'left','XTickLabel',[], 'box', 'off');
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.2);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02, 0.02]);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca,'xtick',[]);
hold off;



%% plot3
% nexttile([1 2])
nexttile
hold on;

% plotting sqrt_V_Uniform_Bound^2
% yline(sqrt_V_Uniform_Bound^2, 'Color', red, 'LineWidth', 1.5);
plot(t(1:crop_index), L(1:crop_index), 'Color', red, 'LineWidth', 3.5);

% plotting lyap_V
for i = 1:size(lyap_V, 1)
    plot(t(1:crop_index), lyap_V(i,1:crop_index), 'Color', colors(i, :), 'LineWidth', 1.5);
end

% title("Lyapunov");
l3 = legend("$\mathcal{L}_{u}^2(\overline{\mathcal{V}}_{1},\overline{\mathcal{V}}_{2})$",'interpreter','latex');
% l3 = legend('show', 'Location', 'northeast');
% Get the current position of the legend
pos3 = l3.Position;
% Adjust the position, for example, moving the legend downward
l3.Position = [pos3(1) - 0.02, pos3(2) - 0.01, pos3(3), pos3(4)];

xlabel('$t$ [s]','interpreter','latex');
ylabel('$V(t)$','interpreter','latex');
xlim([0 t(crop_index)+0.1]);
ylim([0 max(L(1:crop_index))+0.1]);

box on;
% set(gca, 'FontName', 'Arial');
% set(gca, 'FontSize', 12);
% set(gca, 'LineWidth', 1.2);
% set(gca, 'TickDir', 'out');
% set(gca, 'TickLength', [0.02, 0.02]);
% set(gca, 'XMinorTick', 'on');
% set(gca, 'YMinorTick', 'on');
% set(gca,'TickLabelInterpreter','latex');
% set(gca, 'YScale', 'log');

set(gca,'ticklabelinterpreter','latex');
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.2);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02, 0.02]);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'YAxisLocation', 'left', 'box', 'off');
hold off;

% Save PNG
exportgraphics(gcf, fullfile(fig_dir, 'Fig6_ep_ev_LBound.png'), 'Resolution', 600);
% Save EPS
print(gcf, fullfile(fig_dir, 'Fig6_ep_ev_LBound.eps'), '-depsc');

end


