function W = L1GRAPH(feature,rho)

    alpha = 20;
    affine = false;
    CMat = admmOutlier_mat_func(feature,affine,alpha);
%     CMat = admmLasso_mat_func(feature',affine,alpha);
    NN = size(feature,2);  
    C = CMat(1:NN,:);
    W = BuildAdjacency(thrC(C,rho));

end