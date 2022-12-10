function W = L2GRAPH(feature,rho)
    lambda =0.7; K=7;
    coef = []; 
    dat = feature;
    tmp = dat'*dat; 
    if lambda == 0
        Proj_M = pinv(tmp);
    else
        Proj_M = tmp + (lambda) *eye(size(tmp));
        Proj_M =  inv(Proj_M);    
    end
    Q = Proj_M*dat'; clear tmp;
    % --- get coeffient matrixs
    for ii  = 1:size(dat,2)
        stdOrthbasis = zeros(size(dat,2),1);
        stdOrthbasis(ii) = 1;
        tmp1 = stdOrthbasis'* Q *dat(:,ii);
        tmp2 = pinv(stdOrthbasis'* Proj_M * stdOrthbasis);
        tmp = Proj_M * (dat'*dat(:,ii) - (tmp1*tmp2)*stdOrthbasis) ;
        coef = [coef tmp];
    end
    clear ii tmp1 tmp2 tmp;   
    clear Proj_M;
    coef = coef - eye(size(coef)).*coef;
    coef = coef./( repmat(sqrt(sum(coef.*coef)), [size(coef, 1),1]) );
    W = BuildAdjacency(thrC(coef,rho),K);
    
end