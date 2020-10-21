files = dir('*.vtk');
figure;
for i=1:numel(files)
    fname = files(i).name;
    [~,name,ext] = fileparts(fname);
    fid = fopen(fname);
    [X,T] = loadVTKTET(fid);
    
    data = getTetDataRT(T,X,1,0,0);
    BTs = data.triangles(data.boundaryTriangles,:);
    [BX,BT] = minimizeMesh(X,BTs);
%     writeOBJ([name '.obj'],BX,BT);
    
    v1 = BX(BT(:,1),:); v2 = BX(BT(:,2),:); v3 = BX(BT(:,3),:); 
    BTnormals = cross(v1-v2,v2-v3); BTnormals = BTnormals./vecnorm(BTnormals,2,2);
    [~,labels] = max((BTnormals*[eye(3) -eye(3)])');
    
    clf; hold all; axis equal; rotate3d on; fa = 1; title(num2str(i));
    patch('faces',BT(labels==1,:),'vertices',BX,'facecolor','red','facealpha',fa,'edgecolor','none');
    patch('faces',BT(labels==2,:),'vertices',BX,'facecolor','green','facealpha',fa,'edgecolor','none');
    patch('faces',BT(labels==3,:),'vertices',BX,'facecolor','blue','facealpha',fa,'edgecolor','none');
    patch('faces',BT(labels==4,:),'vertices',BX,'facecolor','red','facealpha',fa,'edgecolor','none');
    patch('faces',BT(labels==5,:),'vertices',BX,'facecolor','green','facealpha',fa,'edgecolor','none');
    patch('faces',BT(labels==6,:),'vertices',BX,'facecolor','blue','facealpha',fa,'edgecolor','none');
    
    pause;
end