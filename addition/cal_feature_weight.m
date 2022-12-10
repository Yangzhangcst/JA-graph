function  www_Features= cal_feature_weight(WWW_LAB,WWW_LBP,labels_img)
for i = 1:length(labels_img)

 www_Features{i}  = sqrt(WWW_LAB{1,i}.*WWW_LAB{1,i}+WWW_LBP{1,i}.*WWW_LBP{1,i});
end
