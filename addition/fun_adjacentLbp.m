function [ updata_matrix_www] = fun_adjacentLbp(lbpFeatures,labels_img,seg_edges,seg_lab_vals,para)
P_layNum = size(lbpFeatures,1);
for j = 1: P_layNum
    lbpFeatures_lay = lbpFeatures(j,:);
    lbpFeatures_lay(cellfun(@isempty,lbpFeatures_lay))=[];
    S_layNum = size(lbpFeatures_lay,2);
    TrainData_orig = cell2mat(lbpFeatures_lay);
    TrainData_orig = reshape(TrainData_orig,para.LBP,S_layNum);
%     [coef,score,latent,t2] = pca(TrainData_orig');%pca
%     Contribution_rate  = cumsum(latent)./sum(latent);
%     [ind,~]=find(Contribution_rate<=0.91);
%     pca_TrainData = TrainData_orig'*coef(:,1:ind(end));
    TrainData =TrainData_orig;
    
    for i = 1: S_layNum 
        TrainData(:,i) = TrainData(:,i)/(sqrt(sum(TrainData(:,i).^2)));
        %         TrainData(:,i) = TrainData(:,i)/(max(TrainData(:,i)));
    end
    
    adjacent_matrix = zeros(S_layNum,S_layNum);
    for i = 1:S_layNum
        seg_edges_test = seg_edges(1,j);
        seg_edges_matrix = cell2mat(seg_edges_test);
        %        se_edges_test = reshape(se_edges_test,length(se_edges_test)/2,2);
        Index = find(seg_edges_matrix(:,2)==i);
        adjacent_matrix(i,seg_edges_matrix(Index,1)) = 1;
    end
    adjacent_matrix_cell{j} = adjacent_matrix'+adjacent_matrix + eye(size(adjacent_matrix));
   
    adjacent_matrix_value = ones(S_layNum,S_layNum);
    for i = 1:S_layNum
        ad_test = adjacent_matrix_cell{j};
        ad_test = ad_test - eye (size(ad_test,2),size(ad_test,2));
        index_ad = find(ad_test(i,:) == 1);
        TrainData_test = TrainData(:,index_ad);
        TrainData_sample = TrainData(:,i);
        
        coefficient = lsqnonneg(TrainData_test,TrainData_sample);

       
        for k = 1:length(index_ad)
           
            adjacent_matrix_value(i,index_ad(k)) = sqrt(sum((TrainData(:,i)-coefficient(k).*TrainData_test(:,k)).^2));
        end
        
        
    end
    adjacent_matrix_value = adjacent_matrix_value;%Rij
    WWW_LBP_ad = 1-(adjacent_matrix_value'+adjacent_matrix_value)./2;%weight
    
    if strcmp(para.graphs,'EnRgraph+Dgraph')
         WW = global_graph_combination_all(TrainData,labels_img{j},para);
         WWW_LBP_update = WWW_LBP_ad + WW;
    elseif strcmp(para.graphs,'EnRgraph')
         WW = global_graph_combination_all(TrainData,labels_img{j},para);
         WWW_LBP_update = WW;
    elseif  strcmp(para.graphs,'Dgraph')
         WWW_LBP_update = WWW_LBP_ad;
    elseif  strcmp(para.graphs,'UDgraph')
         w = makeweights(seg_edges{j},seg_lab_vals{j},para.beta);
         WWW_LBP_update = adjacency(seg_edges{j},w);
    elseif strcmp(para.graphs,'EnRgraph+UDgraph')
         w = makeweights(seg_edges{j},seg_lab_vals{j},para.beta);
         W = adjacency(seg_edges{j},w);
         WW = global_graph_combination_all(TrainData,labels_img{j},para);
         WWW_LBP_update = W + WW;
    elseif strcmp(para.graphs,'EnRgraph+Dgraph+UDgraph')
         w = makeweights(seg_edges{j},seg_lab_vals{j},para.beta);
         W = adjacency(seg_edges{j},w);
         WW = global_graph_combination_all(TrainData,labels_img{j},para);
         WWW_LBP_update = W + WW + WWW_LBP_ad;
    end

%     WW = global_graph_combination_all(TrainData,labels_img{j},para);
%     WWW_LBP_update = WWW_LBP_ad + WW;
    updata_matrix_www{j} = WWW_LBP_update;%weight
   

end
