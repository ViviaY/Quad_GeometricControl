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
            omega_n = 6;      % Natural frequency (2–8), larger -> faster convergence
            zeta    = 1.5;    % Damping ratio (0.7–1.2)
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
            b_cbf = MX.sym('b_cbf', 6, 1);     % CBF constraint RHS
            
            % Cost: || mu - mu_nom ||^2
            obj_func = sum((obj.casadi_mu - mu_nom).^2);
            
            % Constraints: A_cbf * mu <= b_cbf
            A_cbf = [ eye(3);  -eye(3) ];
            cbf_constraints = A_cbf * obj.casadi_mu + b_cbf;
            
            % Create QP solver object
            qp_prob = struct('x', obj.casadi_mu, 'f', obj_func, 'g', cbf_constraints, 'p', [mu_nom; b_cbf]);
            qp_opts = struct('print_time', false, 'verbose', false);
            
            % Use OSQP solver through CASADI
            try
                % First try OSQP if available
                obj.casadi_qp = qpsol('qp_solver', 'osqp', qp_prob, qp_opts);
                 fprintf('[CBFQP] Using QP solver: OSQP\n');
            catch
                % Fall back to qpoases
                try
                    obj.casadi_qp = qpsol('qp_solver', 'qpoases', qp_prob, qp_opts);
                    fprintf('[CBFQP] Using QP solver: qpOASES\n');
                catch
                    % Final fallback to built-in IPOPT
                    obj.casadi_qp = qpsol('qp_solver', 'ipopt', qp_prob, qp_opts);
                    warning('Using IPOPT for QP solving. This might be slower than specialized QP solvers.');
                end
            end
        end
        
        function r_ddot = quad_acc(obj, T, xi)
            % Compute and return the acceleration of the quadcopter
            % Inputs:
            %   T - normalized thrust [scalar]
            %   xi - vector of euler angles [3x1] (phi, theta, psi)
            % Output:
            %   r_ddot - acceleration vector [3x1]
            
            GRAV = 9.81;  % Gravity constant
            
            phi = xi(1);
            theta = xi(2);
            psi = xi(3);
            
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
            
            T = v(1);
            xi = v(2:4);
            mu = obj.quad_acc(T, xi);
        end
        
        function v = recover(obj, mu, psi)
            % Recover reduced input v = Psi^-1(mu, psi)
            % Inputs:
            %   mu - virtual input [3x1] (acceleration)
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
            %   z - Reduced state vector [6x1] (position and velocity)
            %   v - Reduced input vector [4x1] (T, phi, theta, psi)
            %   x_r - Reference state vector [9x1]
            %   u_r - Reference input vector [4x1]
            % Output:
            %   v_s - Safe reduced input [4x1]
            
            import casadi.*
            
            % Compute virtualh nominal input
            mu_nom = obj.virtualize(v);
            
            % Safety constraints: K[h; h_dot] + ddh ≥ 0
            % Build K such that K @ [h; h_dot] = k0 * h + k1 * h_dot
            K = [obj.k(1) * eye(6), obj.k(2) * eye(6)];
            
            % Compute CBF value and derivative
            obj.delta = delta_i;
            h =  [z(1:3)-x_r(1:3);x_r(1:3)-z(1:3)] - delta_i;
            h_dot = [z(4:6) - x_r(7:9); x_r(7:9) - z(4:6)];
            eta = [h; h_dot];
            
            % Compute reference acceleration term
            acc_ref = obj.quad_acc(u_r(1), x_r(4:6));
            d2hdt2 = [-acc_ref; acc_ref];
            
            % CBF constraint RHS: K * eta + d2hdt2
            b_cbf = (K * eta + d2hdt2);
            
            % Set up and solve the QP using CASADI
            try
                % Set parameters for the solver
                params = [mu_nom; b_cbf];
                
                % Solve the QP
                % sol = obj.casadi_qp('x0', mu_nom, 'p', params);
                % bounds: g = A*mu - b <= 0  ==> ubg = 0;  lbg = -inf
                lbg = -inf(6,1);
                ubg = zeros(6,1);
            
                sol = obj.casadi_qp('x0', mu_nom, 'p', params, 'lbg', lbg, 'ubg', ubg);
                
                % Extract solution
                mu_opt = full(sol.x);
                catch ME
                    % warning('CBFQP QP failed: %s (ID: %s)', ME.message, ME.identifier);
                    warning(Me.message, "id: %s", ME.identifier);
                    mu_opt = mu_nom;
            end
            
            % Recover safe input
            v_s = obj.recover(mu_opt, v(4));
        end
    end
end