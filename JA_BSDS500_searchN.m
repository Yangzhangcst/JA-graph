clc;clear ; close all;
%% set path
addpath 'addition'
addpath 'algorithms'
addpath 'evals'
addpath (genpath('Graphs/'))
addpath 'Graph_based_segment'
addpath 'msseg'
addpath 'others'

%% set parameters for bipartite graph
para.alpha = 0.01; % affinity between pixels and superpixels
para.beta  =  20;   % scale factor in superpixel affinity
para.nb = 1; % number of neighbors for superpixels
para.rho = 1;
para.K = 4; % number of neighbors for KNN
para.L = 3; % number of neighbors for L0-graph
para.Nimgs = 500; % number of images in BSDS300/500
para.LAB = 3;
para.LBPs = 128;
para.graphs = {'EnRgraph+Dgraph'}; % four sets 'EnRgraph+Dgraph+UDgraph', 'EnRgraph+Dgraph', 'EnRgraph', 'Dgraph'
para.feas = {'Lab+LBP'}; % three sets 'Lab', 'LBP', 'Lab+LBP'
para.gmodes={'ORGEN'};%'knn','L0','L1','L2','LRR','ORPCA','OLRSC'
para.savefig = false;
para.cluster = {'kmeans'};

%% read numbers of segments used in the paper
bsdsRoot ='./database/BSDS500';
saveRoot = 'results_Lab_LBP';
fid = fopen(fullfile('Nsegs_500.txt'),'r');
[BSDS_INFO] = fscanf(fid,'%d %d \n',[2,para.Nimgs]);
fclose(fid);

%% PRI,VoI,GCE,BDE.
PRI_all = zeros(para.Nimgs,1);
VoI_all = zeros(para.Nimgs,1);
GCE_all = zeros(para.Nimgs,1);
BDE_all = zeros(para.Nimgs,1);

%%
Nseg_save = [];
for lab = 1:numel(para.LBPs)
    para.LBP = para.LBPs(lab);
    for n =1:length(para.feas)
        para.fea = para.feas{n};
        for m =1:length(para.gmodes)
            para.gmode = para.gmodes{m};
            for idxI =1:para.Nimgs 
                % read number of segments
                tic;Nseg = BSDS_INFO(2,idxI);
                % locate image
                img_name = int2str(BSDS_INFO(1,idxI));
                img_loc = fullfile(bsdsRoot,'images','test',[img_name,'.jpg']);
                present = 'test';
                if ~exist(img_loc,'file')
                    img_loc = fullfile(bsdsRoot,'images','train',[img_name,'.jpg']);
                    present = 'train';
                    if ~exist(img_loc,'file')
                        img_loc = fullfile(bsdsRoot,'images','val',[img_name,'.jpg']);
                        present = 'val';
                    end
                end
                img = im2double(imread(img_loc));[X,Y,~] = size(img);
                out_path = fullfile(saveRoot,'BSDS500',para.gmode,para.fea,num2str(para.LBP),img_name);
                if ~exist(out_path,'dir'), mkdir(out_path); end
                % generate superpixels
                [para_MS, para_FH] = set_parameters_oversegmentation(img_loc);
                [seg,labels_img,seg_vals,seg_lab_vals,seg_edges,seg_img] = make_superpixels(img_loc,para_MS,para_FH);

                % generate affinity matrix
                [WWW] = adjacent(img_loc,seg_lab_vals,labels_img,seg_edges,para);

                % save over-segmentations
                view_oversegmentation(labels_img,seg_img,out_path,img_name,para.savefig);
                clear labels_img seg_img;

                % build bipartite graph
                B = build_bipartite_graph(img_loc,para,seg,WWW,seg_lab_vals,seg_edges);
                clear seg seg_lab_vals seg_edges;
                img_size = [X,Y];

                out_path_gt= fullfile(saveRoot,'BSDS500',para.gmode,para.fea,num2str(para.LBP),img_name);
                if ~exist(out_path_gt,'dir'), mkdir(out_path_gt); end
                [gt_imgs, gt_cnt] = view_gt_segmentation(bsdsRoot,img,present,out_path_gt,img_name,para);
                nclusters=1:40;  E=zeros(length(nclusters),5);  segs=cell(1,length(nclusters));
                for ncluster=1:length(nclusters)
                    label_img = Tcut(B,nclusters(ncluster),[X,Y],para);
                    % display the result
                    view_segmentation(img,label_img(:),out_path,img_name,para.savefig);
                    % Evaluation and save result
                    out_vals = eval_segmentation(label_img,gt_imgs);
                    E(ncluster,:)=[nclusters(ncluster),out_vals.PRI,out_vals.VoI,out_vals.GCE,out_vals.BDE];
                    segs{ncluster}=uint16(label_img);
                end
                out_seg_path = fullfile(saveRoot,'BSDS500',para.gmode,para.fea,num2str(para.LBP),'segmentation');
                if ~exist(out_seg_path,'dir'), mkdir(out_seg_path); end
                out_seg = fullfile(out_seg_path,[img_name, '.mat']);
                save('-v7',out_seg, 'segs');
                outname = fullfile(out_path,[img_name, '.mat']);
                save('-v7',outname, 'E');
        %         load(outname);
                % Evaluation and save result
                [maxE,idx] = max(E(:,2));ti = toc;
                fprintf('%6d %6s: %2d %9.6f, %9.6f, %9.6f, %9.6f %.4fs\n', idxI,...
                    img_name,E(idx,1), E(idx,2), E(idx,3), E(idx,4), E(idx,5),ti);
                Neg_all(idxI) = E(idx,1);
                PRI_all(idxI) = E(idx,2);  VoI_all(idxI) = E(idx,3);
                GCE_all(idxI) = E(idx,4);  BDE_all(idxI) = E(idx,5);
            end

            fprintf('Mean: %14.6f, %9.6f, %9.6f, %9.6f \n', mean(PRI_all), mean(VoI_all), mean(GCE_all), mean(BDE_all));
            fid_out = fopen(fullfile(saveRoot,'BSDS500',para.gmode,para.fea,num2str(para.LBP),'evaluation.txt'),'w');
            for idxI=1:para.Nimgs
                fprintf(fid_out,'%6d %9.6f, %9.6f, %9.6f, %9.6f \n', BSDS_INFO(1,idxI), PRI_all(idxI), VoI_all(idxI), GCE_all(idxI), BDE_all(idxI));
            end
            fprintf(fid_out,'Mean: %10.6f, %9.6f, %9.6f, %9.6f \n', mean(PRI_all), mean(VoI_all), mean(GCE_all), mean(BDE_all));
            fclose(fid_out);

            fid_out2 = fopen(fullfile(saveRoot,'BSDS500',para.gmode,para.fea,num2str(para.LBP),'Nsegs.txt'),'w');
            for idxI=1:para.Nimgs
                fprintf(fid_out2,'%6d %d \n', BSDS_INFO(1,idxI),Neg_all(idxI));
            end
            fclose(fid_out2);
        end
    end
end