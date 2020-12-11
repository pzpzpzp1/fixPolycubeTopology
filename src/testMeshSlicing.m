clear all; close all;

% fname = 'testMeshes/bunny.obj'; 
fname = 'testMeshes/cubunny.obj'; 
% fname = 'testMeshes/bunny_coarse.obj'; 
% fname = 'testMeshes/bumpy_torus.obj'; 
% fname = 'testMeshes/coarse_sphere.obj'; 
[X, T] = readOBJ(fname); 
X = X*axang2rotm([randn(1,3) 1e-5]); % make generic. makes me nervous.
mesh = processXTtri(X,T);

minZ = min(X(:,3));
maxZ = max(X(:,3));

figure; hold all; rotate3d on; axis equal;
patch('faces',T,'vertices',X,'facealpha',.1,'facecolor',[.9 1.9 .9]/2,'edgecolor','none');
nslices = 20;
zs = linspace(minZ, maxZ, nslices);
% i = round(nslices/2);
% i = nslices-1;
for i = 2:numel(zs)-1
    z = zs(i);
    [U,E,J] = slice_triangles(X,T,[0 0 1 -z]);
    
%     patch('faces',E(:,[1 2 1]),'vertices',U,'facealpha',1,'facecolor','none','edgecolor','cyan','linewidth',2);
%     patch('faces',T(J,:),'vertices',X,'facealpha',1,'facecolor','red','edgecolor','black');
%     figure; hold all; rotate3d on; axis equal;
%     patch('faces',E(:,[1 2 1]),'vertices',U,'facealpha',1,'facecolor','none','edgecolor','black');
%     scatter3(U(:,1),U(:,2),U(:,3),'k.')
    
    [cycles] = getCycles(E);
    [res_digraph, cycles] = orientCycles(cycles, E, U, mesh, J);
    % p = plot(res_digraph,'XData',U(:,1),'YData',U(:,2),'ZData',U(:,3),'ArrowSize',7);
    
    VisualizePC = 0;
    for j=1:numel(cycles)
        cycle = cycles{j};
        cycleData = markCycles(cycle, U, VisualizePC, 0);
        plot3(cycleData(:,1), cycleData(:,2), cycleData(:,3), 'r-');
    end
    
end

