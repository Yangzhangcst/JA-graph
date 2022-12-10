function W = ORPCA_ORGENGRAPH(feature,rho)

    epochs = 10; d = 50;
    Z = feature';
    [pp, nn] = size(Z);
    lambda1 = 1/sqrt(pp);
    lambda2 = 1/sqrt(pp);
    
    LL = randn(pp, d);
    RR = zeros(nn, d);
    AA = zeros(d, d);
    BB = zeros(pp, d);
    for ep=1:epochs
        for t=1:nn
            z = Z(:, t);
            [r, e] = solve_proj2(z, LL, lambda1, lambda2);
            AA = AA + r * r';
            BB = BB + (z-e) * r';
            LL = update_col_orpca(LL, AA, BB, lambda1);
            RR(t, :) = r';
        end
    end
    XX = LL * RR';

    nu0 = 10;  lambda = 0.95;  Nsample = 20;
    if Nsample >= size(XX,1), Nsample = size(XX,1)-1; end
    X = cnormalize_inplace(XX');
    EN_solver =  @(X, y, lambda, nu) rfss( X, y, lambda / nu, (1-lambda) / nu );
    R = ORGEN_mat_func(X, EN_solver, 'nu0', nu0, 'nu_method', 'nonzero',...
        'lambda', lambda, 'Nsample', Nsample, 'maxiter', 10, 'outflag', false);
    W = BuildAdjacency(thrC(R,rho));
end