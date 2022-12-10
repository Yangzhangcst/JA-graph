function LbpFeatures = extractLbp(img_loc,labels_img,nbins ) 

origImage = imread(img_loc);
grayImage = rgb2gray(origImage);
[row , col] = size(grayImage);
grayImage0 = zeros(row+2, col+2);
grayImage0(2:row+1,2:col+1) = grayImage;
SP=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
lbpImage = LBP(grayImage0,SP,0,'i');
% lbpImage =extractLBPFeatures (grayImage0,SP,0,'i');
% imshow(I12); 
A = [];
for m =1:length(labels_img)
    [labels_num(m) ,~ ] = max(max(labels_img{m}));
    for i = 1: labels_num(m)
        [x , y ] = find(labels_img{m}(:,:) == i);
        %  [ lbp_region ]= lbpImage(x(:),y(:));
        for j = 1: size(x)
          B = lbpImage(x(j),y(j));
          A =[A; B];
        end
         [LBPvalue_tongji,~]= hist(double(A),nbins);
         LbpFeatures(m,i) ={LBPvalue_tongji./sum(LBPvalue_tongji)};
        x = 0;
        y = 0;
        A = [ ];
    end
end