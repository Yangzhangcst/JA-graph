function B = build_bipartite_graph(img_loc,para,seg,WW,seg_lab_vals,seg_edges)

for i= 1:length(WW)
   [row_WW,~] = size(WW{1,i});
   WWW = WW{1,i}-eye(row_WW,row_WW);
   WWW(WWW < 0.0001) = 0;
   WW_sparse{i} = sparse(WWW);
end
img = imread(img_loc); 
[X,Y,~] = size(img); 
Np = X*Y; 

% get the overall number of superpixels
Nsp = 0;
for k = 1:length(seg)
    Nsp = Nsp + size(seg{k},2); 
end

W_Y = sparse(Nsp,Nsp); 
edgesXY = []; 
j = 1;
W = WW_sparse;
for k = 1:length(seg) % for each over-segmentation
    
    % affinities between superpixels
    Nk = size(seg{k},2); % number of superpixels in over-segmentation k
    W_Y(j:j+Nk-1,j:j+Nk-1) = prune_knn(W{k},para.nb);
    
    % affinities between pixels and superpixels
    for i = 1:Nk
        idxp = seg{k}{i}; % pixel indices in superpixel i
        Nki = length(idxp); 
        idxsp = j + zeros(Nki,1);
        edgesXY = [edgesXY; [idxp, idxsp]];
        j = j + 1;
    end
end

W_XY = sparse(edgesXY(:,1),edgesXY(:,2),para.alpha,Np,Nsp);

% affinity between a superpixel and itself is set to be the maximum 1.
W_Y(1:Nsp+1:end) = 1;

B = [W_XY;W_Y];

