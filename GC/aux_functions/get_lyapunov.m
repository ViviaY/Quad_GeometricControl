function lyap = get_lyapunov(annealing_output, N, param, e, R, d)
% This script is used to solve the real lyapunov values based on the
% errors. Not the bounds of L, since the bounds L just depends on bar{V}_1 and bar{V}_2. 
% k: the control gain;
% N: the number of the states; which is the len(t) where t is the sampling time;
% e: tracking errors;
% R: tracking rotation matrix
% d: desired trajectory 

% Unpack control gains
k.x = annealing_output.opt_k(1);  
k.v = annealing_output.opt_k(2);  
k.R = annealing_output.opt_k(3);  
k.W = annealing_output.opt_k(4); 

psi_list = zeros(1,N);
for i = 1:N
    psi_list(i) = get_psi(R(:,:,i), d.R(:,:,i));
end

V1 = zeros(1,N);
V2 = zeros(1,N);
for i = 1:N
    psi_i = get_psi(R(:,:,i), d.R(:,:,i));
    [V1(i), V2(i)] = lyapunov(param, k, e.x(:,i), e.v(:,i), e.R(:,i), e.W(:,i), psi_i);
end

V = V1 + V2;
% uniform_V_bound = (Lu(V1(1), V2(1), param))^2;

lyap.V1 = V1;
lyap.V2 = V2;
lyap.V = V;
% lyap.uniform_V_bound = uniform_V_bound;

end