function X0_list = InitialPoints(X0)
% Convert 18xM states [p; v; w; vec(R)] into 12xM states [p; v; xi; w],
% where xi = log(R) ∈ so(3).

    % Ensure size is 18xM (transpose if provided as Mx18)
    if size(X0,1) ~= 18 && size(X0,2) == 18
        X0 = X0.';  % -> 18xM
    end
    assert(size(X0,1) == 18, 'X0 must be 18xM.');

    X0 = double(X0);              % work with numeric doubles
    M  = size(X0,2);

    % Slice blocks
    p0 = X0(1:3,   :);           
    v0 = X0(4:6,   :);            
    w0 = X0(7:9,   :);            
    Rv = X0(10:18, :);            

    % Convert R -> xi = log(R)
    xi0 = zeros(3, M);
    for j = 1:M
        % Rebuild R from its column-stacked vector
        R = reshape(Rv(:,j), 3, 3);   % column-major -> 3x3

        % Project numerically back onto SO(3) for robustness
        [U,~,V] = svd(R);
        R = U*V.';
        if det(R) < 0
            U(:,3) = -U(:,3);
            R = U*V.';
        end
        
        xi0(:,j) = so3_log(R);
    end


    X0_list = [p0; v0; xi0; w0];
end
