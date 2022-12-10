function [W_update] = adjacent(img_loc,seg_lab_vals,labels_img,seg_edges,para)
%%LAB
%
LAB_features = extractLAB(seg_lab_vals,labels_img);
[update_LAB_www] = fun_adjacentLab(LAB_features,labels_img,seg_edges,seg_lab_vals,para);%weight

% for i = 1:length(adjacent_LAB)
%     WWW_LAB{i} = 1-(adjacent_LAB{1,i}'+adjacent_LAB{1,i})./2;%weight
% end



%% LBP
LBP_features = extractLbp(img_loc,labels_img,para.LBP);
[ update_LBP_www] = fun_adjacentLbp(LBP_features,labels_img,seg_edges,seg_lab_vals,para);
% for i = 1:length(adjacent_LBP)
%     WWW_LBP{i} = 1-(adjacent_LBP{1,i}'+adjacent_LBP{1,i})./2;
% end

if strcmp(para.fea,'Lab+LBP')
    W_update =  cal_feature_weight(update_LAB_www,update_LBP_www,labels_img);
elseif strcmp(para.fea,'Lab')
     W_update = update_LAB_www;
elseif  strcmp(para.fea,'LBP')
     W_update = update_LBP_www;
end
