function [x, u, f] = eom_so3(QuadConfig)
    import casadi.*

    x = MX.sym('x', 12);  % [p; v; xi; w]
    u = MX.sym('u', 4);   % [f; tau]

    p = x(1:3);
    v = x(4:6);
    xi = x(7:9);          % Lie algebra rotation vector
    w = x(10:12);

    % Exponential map to SO(3)
    R = so3_exp(xi);

    a = (u(1)/QuadConfig.m) * (R * QuadConfig.e3(:)) - QuadConfig.g * QuadConfig.e3(:);
    w_dot = QuadConfig.J_inv * (u(2:4) - cross(w, QuadConfig.J * w));

    dx = [v;
          a;
          w;
          w_dot];

    f = dx;
end
