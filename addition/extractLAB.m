function     LAB_features = extractLAB(seg_lab_vals,labels_img)
for m =1:length(labels_img)
    [labels_num(m) ,~ ] = max(max(labels_img{m})); 
     img_seg_lab= cell2mat(seg_lab_vals(m));
     A_label = isnan(img_seg_lab);
     img_seg_lab(A_label==1) =0;
    for i = 1: labels_num(m)  
        
        LAB_features(m,i)= {img_seg_lab(i,:)};
%         LAB_features = LAB_features();
    end
end

end