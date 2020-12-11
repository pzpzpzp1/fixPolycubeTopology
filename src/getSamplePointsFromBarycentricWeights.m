function samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T)
    n = size(ws,1);
    
    % 3 x T x 3
    V = zeros(3,size(T,1)*3);
    V(1,:) = reshape(X(T(:,1),:),1,[]);
    V(2,:) = reshape(X(T(:,2),:),1,[]);
    V(3,:) = reshape(X(T(:,3),:),1,[]);
    samplePoints = reshape(ws*V,n,size(T,1),3);
    
end