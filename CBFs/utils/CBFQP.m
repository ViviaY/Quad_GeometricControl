classdef CBFQP
    properties
        k       % CBF safety parameters (k0, k1)
        delta   % Allowable deviation when tracking [m]
        solver  % CASADI QP solver
        
        % State and input variables
        z   % Current reduced state vector
        v   % Nominal reduced input vector
        x_r % Reference state vector
        u_r % Reference input vector
        
        % CASADI optimization variables
        casadi_mu  % Decision variable
        casadi_qp  % QP solver object
    end
    
    methods
        function obj = CBFQP()
            % Constructor - Initialize the CBF-QP controller with CASADI
            import casadi.*
            
            % CBF-QP Parameters            
            omega_n = 6;      % Natural frequency (2-8), larger -> faster convergence
            zeta    = 1.5;    % Damping ratio (0.7-1.2)
            a1 = 2*zeta*omega_n;
            a2 = omega_n^2;

            obj.k = [a1, a2];
       
            % Initialize signals
            obj.z = zeros(6, 1);   % Current reduced state vector
            obj.v = zeros(4, 1);   % Nominal reduced input vector
            obj.x_r = zeros(9, 1); % Reference state vector
            obj.u_r = zeros(4, 1); % Reference input vector
            
            % Initialize CASADI solver for QP
            obj.casadi_mu = MX.sym('mu', 3, 1);  % Decision variable (3D acceleration)
            
            % Dummy parameters (will be updated at runtime)
            mu_nom = MX.sym('mu_nom', 3, 1);   % Nominal virtual input
            b_cbf  = MX.sym('b_cbf', 6, 1);    % CBF constraint RHS
            
            % Cost: || mu - mu_nom ||^2
            obj_func = sum((obj.casadi_mu - mu_nom).^2);
            
            % Constraints: A_cbf * mu <= b_cbf
            % In this formulation, g = A_cbf*mu + b_cbf, with ubg = 0.
            A_cbf = [ eye(3);  -eye(3) ];
            cbf_constraints = A_cbf * obj.casadi_mu + b_cbf;
            
            % Create QP solver object
            qp_prob = struct('x', obj.casadi_mu, ...
                             'f', obj_func, ...
                             'g', cbf_constraints, ...
                             'p', [mu_nom; b_cbf]);
            qp_opts = struct('print_time', false, 'verbose', false);
            
            % Use OSQP solver through CASADI
            try
                obj.casadi_qp = qpsol('qp_solver', 'osqp', qp_prob, qp_opts);
                fprintf('[CBFQP] Using QP solver: OSQP\n');
            catch ME1
                fprintf('[CBFQP] OSQP unavailable: %s\n', ME1.message);
                try
                    obj.casadi_qp = qpsol('qp_solver', 'qpoases', qp_prob, qp_opts);
                    fprintf('[CBFQP] Using QP solver: qpOASES\n');
                catch ME2
                    fprintf('[CBFQP] qpOASES unavailable: %s\n', ME2.message);
                    obj.casadi_qp = qpsol('qp_solver', 'ipopt', qp_prob, qp_opts);
                    fprintf('[CBFQP] Using QP solver: IPOPT via qpsol\n');
                    warning('Using IPOPT for QP solving. This might be slower than specialized QP solvers.');
                end
            end
        end
        
        function r_ddot = quad_acc(obj, T, xi)
            % Compute and return the acceleration of the quadcopter
            % Inputs:
            %   T  - normalized thrust [scalar]
            %   xi - vector of euler angles [3x1] (phi, theta, psi)
            % Output:
            %   r_ddot - acceleration vector [3x1]
            
            GRAV = 9.81;  % Gravity constant
            
            phi   = xi(1);
            theta = xi(2);
            psi   = xi(3);
            
            r_ddot = zeros(3, 1);
            r_ddot(1) = T * (sin(phi) * sin(psi) + cos(phi) * sin(theta) * cos(psi));
            r_ddot(2) = T * (cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi));
            r_ddot(3) = T * cos(phi) * cos(theta) - GRAV;
        end
        
        function mu = virtualize(obj, v)
            % Virtualize reduced input mu = Psi(v)
            % Input:
            %   v - reduced input vector [4x1] (T, phi, theta, psi)
            % Output:
            %   mu - virtual input vector [3x1] (acceleration)
            
            T  = v(1);
            xi = v(2:4);
            mu = obj.quad_acc(T, xi);
        end
        
        function v = recover(obj, mu, psi)
            % Recover reduced input v = Psi^-1(mu, psi)
            % Inputs:
            %   mu  - virtual input [3x1] (acceleration)
            %   psi - associated yaw angle [scalar]
            % Output:
            %   v - reduced input [4x1] (T, phi, theta, psi)
            
            GRAV = 9.81;  % Gravity constant
            
            v = zeros(4, 1);
            v(4) = psi;
            v(3) = atan2(cos(psi) * mu(1) + sin(psi) * mu(2), mu(3) + GRAV);  % theta
            v(2) = atan2(sin(psi) * mu(1) - cos(psi) * mu(2), ...
                   (mu(3) + GRAV) * cos(v(3)));  % phi
            v(1) = (mu(3) + GRAV) / (cos(v(2)) * cos(v(3)));  % T
        end
        
        function v_s = filter(obj, z, v, x_r, u_r, delta_i)
            % Apply safety filter to obtain safe reduced input using CASADI
            % Inputs:
            %   z       - Reduced state vector [6x1] (position and velocity)
            %   v       - Reduced input vector [4x1] (T, phi, theta, psi)
            %   x_r     - Reference state vector [9x1]
            %   u_r     - Reference input vector [4x1]
            %   delta_i - Allowable tracking deviation [scalar]
            % Output:
            %   v_s - Safe reduced input [4x1]
            
            % Profiling counters. Persistent variables are used because CBFQP is
            % a MATLAB value class unless declared as < handle>.
            persistent n_call n_fail
            persistent t_virtualize t_cbf t_pack t_solve t_extract t_recover
            persistent qp_time_total qp_call_count qp_min_time qp_max_time qp_time_sq_total
            persistent printed_stats_fields qp_iter_total qp_iter_count

            if isempty(n_call)
                n_call = 0;
                n_fail = 0;
                t_virtualize = 0;
                t_cbf = 0;
                t_pack = 0;
                t_solve = 0;
                t_extract = 0;
                t_recover = 0;
                qp_time_total = 0;
                qp_call_count = 0;
                qp_min_time = inf;
                qp_max_time = 0;
                qp_time_sq_total = 0;
                printed_stats_fields = false;
                qp_iter_total = 0;
                qp_iter_count = 0;
            end

            n_call = n_call + 1;
            
            %% 1. Compute virtual nominal input
            tt = tic;
            mu_nom = obj.virtualize(v);
            t_virtualize = t_virtualize + toc(tt);
            
            %% 2. Build CBF constraints
            tt = tic;

            % Safety constraints: K[h; h_dot] + ddh >= 0
            % Build K such that K @ [h; h_dot] = k0 * h + k1 * h_dot
            K = [obj.k(1) * eye(6), obj.k(2) * eye(6)];
            
            % Keep original assignment for compatibility with existing code.
            obj.delta = delta_i;

            h = [z(1:3)-x_r(1:3); x_r(1:3)-z(1:3)] - delta_i;
            h_dot = [z(4:6) - x_r(7:9); x_r(7:9) - z(4:6)];
            eta = [h; h_dot];
            
            % Compute reference acceleration term
            acc_ref = obj.quad_acc(u_r(1), x_r(4:6));
            d2hdt2 = [-acc_ref; acc_ref];
            
            % CBF constraint RHS
            b_cbf = K * eta + d2hdt2;

            t_cbf = t_cbf + toc(tt);
            
            %% 3. Pack parameters and bounds
            tt = tic;
            params = [mu_nom; b_cbf];
            lbg = -inf(6,1);
            ubg = zeros(6,1);
            t_pack = t_pack + toc(tt);
            
            %% 4. Solve QP
            try
                tt = tic;
                sol = obj.casadi_qp('x0', mu_nom, ...
                                    'p', params, ...
                                    'lbg', lbg, ...
                                    'ubg', ubg);
                qp_elapsed = toc(tt);

                % Correctly accumulate QP solve time for filter profile.
                t_solve = t_solve + qp_elapsed;

                % Separate QP-only statistics.
                qp_time_total = qp_time_total + qp_elapsed;
                qp_call_count = qp_call_count + 1;
                qp_min_time = min(qp_min_time, qp_elapsed);
                qp_max_time = max(qp_max_time, qp_elapsed);
                qp_time_sq_total = qp_time_sq_total + qp_elapsed^2;

                % Print stats fields once, and try to accumulate iteration count
                % if the CasADi plugin exposes one.
                st = obj.casadi_qp.stats();
                if ~printed_stats_fields
                    fprintf('\n[CBFQP] Solver stats fields after first solve:\n');
                    disp(fieldnames(st));
                    printed_stats_fields = true;
                end

                iter = -1;
                if isfield(st, 'iter_count')
                    iter = st.iter_count;
                elseif isfield(st, 'iter')
                    iter = st.iter;
                elseif isfield(st, 'iterations')
                    iter = st.iterations;
                end
                if isnumeric(iter) && isscalar(iter) && iter >= 0
                    qp_iter_total = qp_iter_total + iter;
                    qp_iter_count = qp_iter_count + 1;
                end

                %% 5. Extract solution
                tt = tic;
                mu_opt = full(sol.x);
                t_extract = t_extract + toc(tt);

            catch ME
                n_fail = n_fail + 1;
                % Keep original fallback behavior: use nominal virtual input.
                % Avoid warning spam during profiling because warnings are slow.
                % warning('CBFQP:QPFailed', 'CBFQP QP failed: %s', ME.message);
                mu_opt = mu_nom;
            end
            
            %% 6. Recover safe input
            tt = tic;
            v_s = obj.recover(mu_opt, v(4));
            t_recover = t_recover + toc(tt);

            %% 7. Print profiling periodically
            if mod(n_call, 1000) == 0
                profile_total = t_virtualize + t_cbf + t_pack + t_solve + t_extract + t_recover;

                if qp_call_count > 1
                    qp_mean = qp_time_total / qp_call_count;
                    qp_var = max(qp_time_sq_total / qp_call_count - qp_mean^2, 0);
                    qp_std = sqrt(qp_var);
                elseif qp_call_count == 1
                    qp_mean = qp_time_total;
                    qp_std = 0;
                else
                    qp_mean = NaN;
                    qp_std = NaN;
                end

                if qp_iter_count > 0
                    avg_iter = qp_iter_total / qp_iter_count;
                else
                    avg_iter = -1;
                end

                fprintf('\n[CBFQP filter profile] calls = %d, fails = %d\n', n_call, n_fail);
                fprintf('  virtualize: %.6f ms/call\n', 1000 * t_virtualize / n_call);
                fprintf('  build CBF:  %.6f ms/call\n', 1000 * t_cbf / n_call);
                fprintf('  pack args:  %.6f ms/call\n', 1000 * t_pack / n_call);
                fprintf('  QP solve:   %.6f ms/call\n', 1000 * t_solve / max(qp_call_count, 1));
                fprintf('  extract:    %.6f ms/call\n', 1000 * t_extract / max(qp_call_count, 1));
                fprintf('  recover:    %.6f ms/call\n', 1000 * t_recover / n_call);
                fprintf('  total:      %.6f ms/call\n', 1000 * profile_total / n_call);
                fprintf('  avg iter:   %.3f iterations/call\n', avg_iter);

                fprintf('\n[QP Solver Statistics (OSQP via CasADi)]\n');
                fprintf('  Total solves:        %d\n', qp_call_count);
                fprintf('  Total QP time:       %.4fs\n', qp_time_total);
                fprintf('  Avg per solve:       %.3f ms\n', 1000 * qp_mean);
                fprintf('  Min/Max per solve:   %.3f / %.3f ms\n', 1000 * qp_min_time, 1000 * qp_max_time);
                fprintf('  Std dev:             %.3f ms\n', 1000 * qp_std);
                fprintf('  QP fail count:       %d\n', n_fail);
                fprintf('\n');
            end
        end
    end
end
