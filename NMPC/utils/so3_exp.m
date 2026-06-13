function R = so3_exp(xi)
% Exponential map so(3) -> SO(3)
% Works for double and CasADi SX/MX
    import casadi.*  % 仅在使用 CasADi 时生效，数值分支不依赖它

    xi = xi(:);
    I  = eye(3);
    wx = [   0    -xi(3)  xi(2);
           xi(3)    0    -xi(1);
          -xi(2)  xi(1)    0  ];

    theta2 = xi.'*xi;  % 对 double、SX、MX 都安全
    is_sym = isa(theta2,'casadi.SX') || isa(theta2,'casadi.MX');

    if is_sym
        % —— CasADi 分支：用 if_else，两个分支都必须数值安全 ——
        theta = sqrt(theta2);
        epsd  = 1e-12;    % 防除零
        t0    = 1e-4;     % 小角级数切换阈值

        A_small = 1 - theta.^2/6 + theta.^4/120;
        B_small = 0.5 - theta.^2/24 + theta.^4/720;

        A_full  = sin(theta)./(theta + epsd);
        B_full  = (1 - cos(theta))./(theta.^2 + epsd);

        use_small = theta < t0;    % 这是 CasADi 的逻辑表达式
        A = if_else(use_small, A_small, A_full);
        B = if_else(use_small, B_small, B_full);
    else
        % —— 数值分支：普通 if/else，无需 if_else ——
        th = sqrt(double(theta2));
        if th < 1e-8
            A = 1 - th^2/6;
            B = 0.5 - th^2/24;
        else
            A = sin(th)/th;
            B = (1 - cos(th))/(th^2);
        end
    end

    R = I + A*wx + B*(wx*wx);
end
