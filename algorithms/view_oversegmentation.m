function view_oversegmentation(label_img,seg_img,out_path,only_name,iswrite)

if ~exist(fullfile(out_path,'segs'),'dir'), mkdir(fullfile(out_path,'segs')); end

%% make the resulted image with red boundaries
for i=1:size(label_img,2)
    [imgMasks,segOutline,imgMarkup]=segoutput(seg_img{i},double(label_img{i}));
    if iswrite, imwrite(imgMarkup,fullfile(out_path,'segs', [only_name, '_', int2str(i) '.bmp'])); end
    clear imgMasks segOutline imgMarkup;
end