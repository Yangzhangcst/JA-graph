function W = KNNGRAPH(feature,para)
    feature = feature';
    [~, edges]=knn(feature,para.K);
    w = makeweights(edges,feature,para.beta);
    W = adjacency(edges,w);
    
end