function label_img = Tcut(B,Nseg,img_size,para)

% B - |X|-by-|Y|, cross-affinity-matrix
% note that |X| = |Y| + |I|

[Nx,Ny] = size(B);
if Ny < Nseg
    error('Need more superpixels!'); 
end

%%% build the superpixel graph
dx = sum(B,2);
Dx = sparse(1:Nx,1:Nx,1./dx);
Wy = B'*Dx*B;

%%% compute Ncut eigenvectors
% normalized affinity matrix
d = sum(Wy,2);
D = sparse(1:Ny,1:Ny,1./sqrt(d));
nWy = D*Wy*D;
nWy = (nWy+nWy')/2;

nWy(isnan(nWy))=0;

% computer eigenvectors
[evec,eval] = eig(full(nWy)); % use eigs for large superpixel graphs  
[~,idx] = sort(diag(eval),'descend');
Ncut_evec = D*evec(:,idx(1:Nseg));

%%% compute the Ncut eigenvectors on the entire bipartite graph (transfer!)
evec = Dx * B * Ncut_evec;

%%% k-means clustering
% extract spectral representations for pixels
evec = evec(1:prod(img_size),:);

% normalize each row to unit norm
evec = bsxfun( @rdivide, evec, sqrt(sum(evec.*evec,2)) + 1e-10 );

% k-means
if strcmp(para.cluster,'kmeans')
    labels = k_means(evec',Nseg);
elseif strcmp(para.cluster,'kmeanspp')
    labels = k_meanspp(evec', Nseg);
elseif  strcmp(para.cluster,'litekmeans')
    Center = clusteringcenter(evec',Nseg);
    labels = litekmeans(evec,Nseg,'Start',Center);
end

label_img = reshape(labels,img_size);
