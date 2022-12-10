function W = ORGENGRAPH_PAMA(feature,nu0,lambda,Nsample,rho)
    feature = feature';
    if Nsample >= size(feature,1), Nsample = size(feature,1)-1; end
    X = cnormalize_inplace(feature');
    EN_solver =  @(X, y, lambda, nu) rfss( X, y, lambda / nu, (1-lambda) / nu );
    R = ORGEN_mat_func(X, EN_solver, 'nu0', nu0, 'nu_method', 'fixed',...
        'lambda', lambda, 'Nsample', Nsample, 'maxiter', 10, 'outflag', false);
    W = BuildAdjacency(thrC(R,rho));
    
end