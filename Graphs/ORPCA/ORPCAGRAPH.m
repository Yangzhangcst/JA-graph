function W = ORPCAGRAPH(feature)
    
    epochs = 10; d = 50; %default 60
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
    [U,S,V] = svds(XX,d);
    S = diag(S);
    rr = sum(S>1e-4*S(1));
    U = U(:,1:rr);S = S(1:rr);
    U = U*diag(sqrt(S));
    W = (U*U').^4;

end