function annealing_output = annealing(anneal_options, param)
rng("default");

initial_solution = anneal_options.initial_solution; % Initial solution guess
lower_bound = anneal_options.lower_bound; % Lower bounds for variables
upper_bound = anneal_options.upper_bound; % Upper bounds for variables

total_runtime = 0;
runtime_n = anneal_options.runtime_n;

optimal.x = zeros(runtime_n, size(initial_solution,2));
optimal.fval = zeros(runtime_n, 1);

options = optimoptions(@simulannealbnd, 'MaxFunctionEvaluations', 50000);

f = @(x)position_bound(x, anneal_options, param);

% Call the Simulated Annealing function
for i = 1:runtime_n
    tic;
    [x, fval, ~, ~] = simulannealbnd(f, initial_solution, lower_bound, upper_bound, options);
    runtime = toc;
    total_runtime = total_runtime + runtime;

    optimal.x(i,:) = x;
    optimal.fval(i) = fval;
end

ave_runtime = total_runtime / runtime_n;
opt_fval = min(optimal.fval);
opt_k = optimal.x(find(optimal.fval == min(optimal.fval)),:);
opt_bounds = evaluate(opt_k, anneal_options, param);

others.annealing_ave_runtime = ave_runtime;
others.annealing_min_fval = opt_fval;

annealing_output.opt_k = opt_k;
annealing_output.opt_bounds = opt_bounds;
annealing_output.others = others;
% sorted_opt_fval = sort(optimal.fval);
% for i=1:length(sorted_opt_fval)
%     opt_k_sort = optimal.x(find(optimal.fval == sorted_opt_fval(i)),:);
% 
%     disp(['sorted_opt_feval = ', num2str(sorted_opt_fval(i), '%.6f'), ', opt_k_sort = ', mat2str(opt_k_sort)]);
% 
% end

% Display results
fprintf('Total time: %.3f, Average runtime of optimization:%.3f \n', total_runtime, ave_runtime);
fprintf('Minimum value among those n trials: %.3f \n', opt_fval);
disp('Optimal gains:');
disp(opt_k);
disp('Optimal bounds and parameters: Lp, Lv, Lf, F_bound, c1, c2');
disp([opt_bounds.Lp, opt_bounds.Lv, opt_bounds.Lf, opt_bounds.F_bound, opt_bounds.c1, opt_bounds.c2]);

end
