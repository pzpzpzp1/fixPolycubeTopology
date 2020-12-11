function labels = getLabelsFromXT(X,T)
    v1 = X(T(:,1),:); v2 = X(T(:,2),:); v3 = X(T(:,3),:); 
    BTnormals = cross(v1-v2,v2-v3); BTnormals = BTnormals./vecnorm(BTnormals,2,2);
    load labelDirs.mat;
    [~,labels] = max((BTnormals*labelDirs)');
end