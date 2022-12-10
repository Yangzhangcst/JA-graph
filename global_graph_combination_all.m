function W = global_graph_combination_all(feature,labels_img,para)

% for each over-segmentation
status = regionprops(labels_img,'basic');
Area = cat(1,status.Area);
centroid = cat(1,status.Centroid);


% please choose different kinds of global graph combination
switch para.gmode
    case 'knn'
        W = KNNGRAPH(feature,para);
    case 'knn+orgen'
        W1 = KNNGRAPH(feature,para);
        W2 = ORGENGRAPH(feature,para.rho);
%         W = W1+W2; 
        W = sqrt(W1.*W1+W2.*W2);
    case 'L0'
        W = L0GRAPH(feature,para.L,centroid,Area);
%         W = fun_l0(feature,para.L);
    case 'L1'
        W = L1GRAPH(feature,para.rho);
    case 'L1+orgen'
        W1 = L1GRAPH(feature,para.rho);
        W2 = ORGENGRAPH(feature,para.rho);
%         W = W1+W2; 
        W = sqrt(W1.*W1+W2.*W2);
    case 'L2'
        W = L2GRAPH(feature,para.rho);
    case 'L2+orgen'
        W1 = L2GRAPH(feature,para.rho);
        W2 = ORGENGRAPH(feature,para.rho);
%         W = W1+W2; 
        W = sqrt(W1.*W1+W2.*W2);
    case 'LRR'
        W = LRRGRAPH(feature);     
    case 'LRR+ORGEN'
        W1 = LRRGRAPH(feature);
        W2 = ORGENGRAPH(feature,para.rho);
        %         W = W1+W2; 
        W = sqrt(W1.*W1+W2.*W2);
    case 'LRR-ORGEN'
        W = LRR_ORGENGRAPH(feature,para.rho);
    case 'ORPCA'
        W = ORPCAGRAPH(feature);
    case 'ORPCA+ORGEN'
        W1 = ORPCAGRAPH(feature);
        W2 = ORGENGRAPH(feature,para.rho); 
        W = sqrt(W1.*W1+W2.*W2);
    case 'OLRSC'
        W = OLRSCGRAPH(feature,para.rho);    
    case 'OLRSC+ORGEN'
        W1 = OLRSCGRAPH(feature,para.rho); 
        W2 = ORGENGRAPH(feature,para.rho); 
        W = sqrt(W1.*W1+W2.*W2);
    case 'ORGEN'
        W = ORGENGRAPH(feature,para.rho);
    case 'ORGEN_nu10_lab0.95_n20_fixed'
        W = ORGENGRAPH_PAMA(feature,10,0.95,20,para.rho);
    case 'ORGEN_nu2_lab0.95_n20_fixed'
        W = ORGENGRAPH_PAMA(feature,2,0.95,20,para.rho);
    case 'ORGEN_nu5_lab0.95_n20_fixed'
        W = ORGENGRAPH_PAMA(feature,5,0.95,20,para.rho);
    case 'ORGEN_nu15_lab0.95_n20_fixed'
        W = ORGENGRAPH_PAMA(feature,15,0.95,20,para.rho);
    case 'ORGEN_nu20_lab0.95_n20_fixed'
        W = ORGENGRAPH_PAMA(feature,20,0.95,20,para.rho);
    case 'ORGEN_nu25_lab0.95_n20_fixed'
        W = ORGENGRAPH_PAMA(feature,25,0.95,20,para.rho);
    case 'ORGEN_nu30_lab0.95_n20_fixed'
        W = ORGENGRAPH_PAMA(feature,30,0.95,20,para.rho);
    case 'ORGEN_nu35_lab0.95_n20_fixed'
        W = ORGENGRAPH_PAMA(feature,35,0.95,20,para.rho);
    case 'ORPCA-ORGEN'
        W = ORPCA_ORGENGRAPH(feature,para.rho);         
end
