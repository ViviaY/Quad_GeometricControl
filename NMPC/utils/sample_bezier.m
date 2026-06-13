function Bezier = sample_bezier(Bezier_data, track_tspan, ifPlot)
% Load Bezier data from a MATLAB .mat file and sample position, velocity,
% and acceleration along the trajectory.

    Bezier = struct();
    Points_Array = Bezier_data.Points_Array;
    t_Array = Bezier_data.t_Array(:).';   % 1 x (K+1)

    % Reorder control points to (segments, n+1, 3)
    sz = size(Points_Array);
    if numel(sz) ~= 3
        error('Points_Array must be a 3D array, got size %s', mat2str(sz));
    end

    if sz(1) == 3
        % Case: (3, n+1, K) -> (K, n+1, 3)
        Points_Array = permute(Points_Array, [3, 2, 1]);
    elseif sz(3) == 3
        % Case: already (K, n+1, 3)
    else
        error('Unexpected Points_Array format: %s', mat2str(size(Points_Array)));
    end

    [K, n1, d] = size(Points_Array);
    if d ~= 3
        error('Last dimension of Points_Array must be 3, got %d', d);
    end
    if numel(t_Array) ~= K+1
        error('t_Array length must be segments+1, got %d (segments=%d)', numel(t_Array), K);
    end

    % Prepare sampling times 
    tt = track_tspan(:);  % N x 1
    N = numel(tt);

    % Map global time to segment index and local time
    [segIdx, tLocal] = map_time_to_segment(tt, t_Array);
    dur = t_Array(segIdx+1).' - t_Array(segIdx).';
    dur = max(dur, 1e-12);   % avoid zero division

    % Evaluate Bezier segments 
    position     = zeros(N, 3);
    velocity     = zeros(N, 3);
    acceleration = zeros(N, 3);

    for i = 1:N
        Pk = squeeze(Points_Array(segIdx(i), :, :));  % (n+1) x 3
        u  = tLocal(i);
        [p, v, a] = bezier_with_derivatives(Pk, u);

        % Scale derivatives by duration
        position(i, :)     = p;
        velocity(i, :)     = v ./ dur(i);
        acceleration(i, :) = a ./ (dur(i)^2);
    end

    % Optional plot
    if ifPlot
        figure; hold on; grid on; axis equal;
        title('Bezier trajectory sampling');
        xlabel('x'); ylabel('y'); zlabel('z');

        % Plot control polygons
        for k = 1:K
            Pk = squeeze(Points_Array(k, :, :));
            plot3(Pk(:,1), Pk(:,2), Pk(:,3), 'o-');
        end
        % Plot sampled trajectory
        plot3(position(:,1), position(:,2), position(:,3), 'LineWidth', 2);
        legend('Control polygons','Sampled trajectory');
        view(3);
    end
    
    % Output
    Bezier.position     = position;
    Bezier.velocity     = velocity;
    Bezier.acceleration = acceleration;
end

%% Map global time to segment index and normalized local time
function [segIdx, tLocal] = map_time_to_segment(t, tArray)
% Inputs:
%   tt     : N×1 sample times
%   tArray : 1×(K+1) segment boundaries
% Outputs:
%   segIdx : N×1, segment index in [1..K]
%   tLocal : N×1, normalized local time in [0,1]

    tArray = tArray(:).';   % 1×(K+1)
    K = numel(tArray) - 1;

    segIdx = discretize(t, tArray, 'IncludedEdge','right');  % 1..K or NaN

    % Handle out-of-bound times
    segIdx(isnan(segIdx) & t <  tArray(1))  = 1;
    segIdx(isnan(segIdx) & t >= tArray(end)) = K;
    segIdx = max(1, min(segIdx, K));

    tStart = tArray(segIdx).';
    tEnd   = tArray(segIdx+1).';
    dur    = max(tEnd - tStart, 1e-12);

    tLocal = (t - tStart) ./ dur;
    tLocal = min(max(tLocal, 0), 1);
end

%% Evaluate Bezier curve and its derivatives at parameter u
function [p, v, a] = bezier_with_derivatives(P, t)
% Inputs:
%   P : (n+1)×3 control points
%   u : scalar in [0,1]
% Outputs:
%   p : 1×3 position
%   v : 1×3 derivative w.r.t u
%   a : 1×3 second derivative w.r.t u

    n = size(P,1) - 1;

    p = bernstein(P, t);

    if n >= 1
        P1 = n * (P(2:end,:) - P(1:end-1,:));
        v  = bernstein(P1, t);
    else
        v = zeros(1,3);
    end

    if n >= 2
        P2 = n*(n-1) * (P(3:end,:) - 2*P(2:end-1,:) + P(1:end-2,:));
        a  = bernstein(P2, t);
    else
        a = zeros(1,3);
    end
end

%% Evaluate Bezier curve (Bernstein form) at parameter u
function val = bernstein(P, t)
% Inputs:
%   P : (n+1)×d control points
%   u : scalar in [0,1]
% Outputs:
%   val : 1×d evaluated point
    m = size(P,1);
    n = m - 1;
    t  = min(max(t,0),1);
    um = 1 - t;

    B = zeros(m,1);
    for k = 0:n
        B(k+1) = nchoosek(n,k) * (um)^(n-k) * (t)^k;
    end
    val = (B.' * P);
end
