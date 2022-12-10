function W = OLRSCGRAPH(feature,rho)

    epochs = 10;
    d = 40;
    Z = feature;
    [pp, nn] = size(Z);
    lambda1 = 1;
    lambda2 = 1/sqrt(pp);
    lambda3_base = 1/sqrt(pp);
    MM = zeros(pp, d);
    AA = zeros(d, d);
    BB = zeros(pp, d);
    UU = zeros(nn, d);
    VV = zeros(nn, d);
    DD = randn(pp, d);
    for ep=1:epochs
        for t=1:nn
            z = Z(:, t);
            lambda3 = sqrt(t) * lambda3_base;
            [v, e] = OLRR_solve_ve(z, DD, lambda1, lambda2);
            normz = norm(z);
            u = (DD - MM)' * z / (normz * normz + 1/lambda3);
            MM = MM + z * u';
            AA = AA + v * v';
            BB = BB + (z-e) * v';
            DD = OLRR_solve_D(DD, MM, AA, BB, lambda1, lambda3);
            UU(t, :) = u';
            VV(t, :) = v';
        end
        MM = zeros(pp, d);
    end
    X = UU*VV';
    W= BuildAdjacency(thrC(X,rho));

end