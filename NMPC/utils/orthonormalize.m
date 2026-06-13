function R = orthonormalize(R)
% Orthonormalize a 3x3 matrix via polar decomposition (one-step).
% Ensures R is a proper rotation with det ~ +1.
    [U,~,V] = svd(R);
    R = U*V.';
    if det(R) < 0
        % Fix improper rotation
        U(:,3) = -U(:,3);
        R = U*V.';
    end
end