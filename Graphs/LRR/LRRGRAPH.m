function W = LRRGRAPH(feature)
    lambda = 0.18;
    Z = solve_lrr(feature,lambda,0);
%     Z = RKLRR(feature,0.18);
%     Z = IBDLR(feature*feature',3,3,5,0);
    [U,S,V] = svd(Z,'econ');
    S = diag(S);
    r = sum(S>1e-4*S(1));
    U = U(:,1:r);S = S(1:r);
    U = U*diag(sqrt(S));
    W = (U*U').^4;

end