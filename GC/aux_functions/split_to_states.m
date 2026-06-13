function [x, v, W, R] = split_to_states(X)

x = X(1:3);
v = X(4:6);
W = X(7:9);
R = reshape(X(10:18), 3, 3);
end