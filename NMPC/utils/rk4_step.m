function x_next = rk4_step(f_dyn, x, u, dt)
% One RK4 step consistent with model f_dyn
    k1 = full(f_dyn(x,                 u));
    k2 = full(f_dyn(x + 0.5*dt*k1,     u));
    k3 = full(f_dyn(x + 0.5*dt*k2,     u));
    k4 = full(f_dyn(x +     dt*k3,     u));
    x_next = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end