function W = LRR_ORGENGRAPH(feature,rho)

    lam = 0.35;
    Z = solve_lrr(feature,lam,0);
    
    nu0 = 10;  lambda = 0.95;  Nsample = 20;
    if Nsample >= size(Z,1), Nsample = size(Z,1)-1; end
    X = cnormalize_inplace(Z');
    EN_solver =  @(X, y, lambda, nu) rfss( X, y, lambda / nu, (1-lambda) / nu );
    R = ORGEN_mat_func(X, EN_solver, 'nu0', nu0, 'nu_method', 'nonzero',...
        'lambda', lambda, 'Nsample', Nsample, 'maxiter', 10, 'outflag', false);
    W = BuildAdjacency(thrC(R,rho));