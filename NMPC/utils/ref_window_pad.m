function Xref_win = ref_window_pad(ref, t_idx, H, dt)
    % ref: 12×T
    NX = size(ref,1); T = size(ref,2);
    if t_idx+H <= T
        Xref_win = ref(:, t_idx:t_idx+H);
        return;
    end
    
    need = t_idx + H - T; % Number of steps to be padded
    Xref_win = [ref(:, t_idx:T), zeros(NX, need)];

    % Constant-velocity / zero-acceleration extrapolation
    p  = ref(1:3, T);   v  = ref(4:6, T);
    xi = ref(7:9, T);   w  = ref(10:12, T);
    for i = 1:need
        p  = p + v*dt;   
        % Velocity is kept constant
        Xref_win(:, T - t_idx + 1 + i) = [p; v; xi; w];
    end
end