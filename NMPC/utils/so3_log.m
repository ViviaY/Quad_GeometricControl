function xi = so3_log(R)
% SO(3) -> so(3) logarithm map
% Works for double and casadi.SX/MX (AD-friendly)
    is_sym = isa(R,'casadi.SX') || isa(R,'casadi.MX');
    epsd   = 1e-12;

    if is_sym
        import casadi.*
        % clamp trace → c ∈ [-1,1]
        tr = R(1,1)+R(2,2)+R(3,3);
        c  = (tr - 1)/2;
        c  = min(1, max(-1, c));
        theta = acos(c);

        S   = R - R.';
        vee = [S(3,2); S(1,3); S(2,1)];

        % small-angle branch
        xi_small = 0.5 * vee;

        % generic branch
        xi_gen = (theta ./ (2*(sin(theta)+epsd))) .* vee;

        % near-pi branch: extract axis robustly
        near_pi = (1 + c) < 1e-6;
        vx = sqrt(max(0, (R(1,1)+1)/2));
        vy = sqrt(max(0, (R(2,2)+1)/2));
        vz = sqrt(max(0, (R(3,3)+1)/2));

        choose_x = vx >= vy;
        ax_x = [ vx;
                 (R(1,2)+R(2,1))/(4*vx + epsd);
                 (R(1,3)+R(3,1))/(4*vx + epsd) ];
        ax_y = [ (R(1,2)+R(2,1))/(4*vy + epsd);
                  vy;
                 (R(2,3)+R(3,2))/(4*vy + epsd) ];
        ax_xy = if_else(choose_x, ax_x, ax_y);
        choose_z = (if_else(choose_x, vx, vy)) >= vz;
        ax = if_else(choose_z, ax_xy, ...
             [ (R(1,3)+R(3,1))/(4*vz + epsd);
               (R(2,3)+R(3,2))/(4*vz + epsd);
                vz ]);

        xi_pi = pi * ax;

        small = theta < 1e-6;
        xi = if_else(small, xi_small, if_else(near_pi, xi_pi, xi_gen));

    else
        % ---- numeric branch ----
        c = (trace(R)-1)/2; c = max(-1,min(1,c));
        theta = acos(c);
        S = R - R.'; vee = [S(3,2); S(1,3); S(2,1)];

        if theta < 1e-8
            xi = 0.5 * vee;
        elseif (1 + c) < 1e-6  % near pi
            vx = sqrt(max(0,(R(1,1)+1)/2));
            vy = sqrt(max(0,(R(2,2)+1)/2));
            vz = sqrt(max(0,(R(3,3)+1)/2));
            if vx >= vy && vx >= vz
                ax = [ vx;
                      (R(1,2)+R(2,1))/(4*vx + epsd);
                      (R(1,3)+R(3,1))/(4*vx + epsd) ];
            elseif vy >= vz
                ax = [ (R(1,2)+R(2,1))/(4*vy + epsd);
                       vy;
                      (R(2,3)+R(3,2))/(4*vy + epsd) ];
            else
                ax = [ (R(1,3)+R(3,1))/(4*vz + epsd);
                      (R(2,3)+R(3,2))/(4*vz + epsd);
                       vz ];
            end
            xi = pi * ax;
        else
            xi = (theta / (2*(sin(theta)+epsd))) * vee;
        end
    end
end
