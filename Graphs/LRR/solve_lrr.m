function [Z,E] = solve_lrr(X,lambda,display)
Q = orth(X');
A = X*Q;
[Z,E] = lrra(X,A,lambda,display);
Z = Q*Z;