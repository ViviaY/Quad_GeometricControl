function v = vee3(S)
% Map a 3x3 skew-symmetric matrix to R^3.
% S = [  0  -v3  v2
%       v3    0 -v1
%      -v2   v1   0]
    v = [S(3,2); S(1,3); S(2,1)];
end