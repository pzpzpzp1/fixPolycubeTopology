% clear all; close all;

% parameters
n = 10; % barycentric points
% face corner edge
% scl = [.8,   1,   1]; % isolated faces. edge and corners
% scl = [.9,   1.1,   1]; % isolated edges
scl = [.9,   1.2,   1.05]; % corners centric
% scl = [0,   1,   0]; % corners only

%% load mesh
fname = 'testMeshes/coarse_sphere.obj'; % sphere for gauss map vis
% fname = 'testMeshes/ex1.obj'; % simple ramp on block
% fname = 'testMeshes/ex17.obj'; % actual polycube
% fname = 'testMeshes/ex18.obj'; % wierd vertex
% fname = 'testMeshes/ex6.obj'; % p
% fname = 'testMeshes/ex7.obj'; % p
[X, T] = readOBJ(fname);

%% process mesh
mesh = processXTtri(X,T);

dt = 5
% dt = 0;
processedVertNormals = (-mesh.cotlaplacian*dt + speye(mesh.numVertices))\mesh.vertNormals;
% figure; clf; hold all; axis equal; rotate3d on; fa = 1; 
% patch('faces',T,'vertices',X,'facecolor','red','facealpha',fa,'edgecolor','none');
% quiver3(mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3),processedVertNormals(:,1),processedVertNormals(:,2),processedVertNormals(:,3))

ws = getBarycentricSamplingWeights(n);
% n, nX, 3
triverts = reshape(getSamplePointsFromBarycentricWeights(ws, X, T),[],3);
trivertNormals = reshape(getSamplePointsFromBarycentricWeights(ws, processedVertNormals, T),[],3);
trivertNormals = trivertNormals./vecnorm(trivertNormals,2,2);

faceDirections = [eye(3); -eye(3);];
cornerDirections = ((dec2bin([0:7]))-48.5)*2; cornerDirections = cornerDirections./vecnorm(cornerDirections,2,2);
edgeDirections = dec2base(0:26,3) - 49; edgeDirections = edgeDirections(sum(abs(edgeDirections)')==2,:); edgeDirections  = edgeDirections ./vecnorm(edgeDirections ,2,2);
primaryDirections = [faceDirections*scl(1); cornerDirections*scl(2); edgeDirections*scl(3)];

[~, labels] = max((trivertNormals * primaryDirections')');

dcolors = distinguishable_colors(size(primaryDirections,1));
figure; clf; hold all; axis equal; rotate3d on; fa = 1; 
for i=1:size(primaryDirections,1)
    bucketinds = labels==i;
    scatter3(triverts(bucketinds,1),triverts(bucketinds,2),triverts(bucketinds,3), 20, dcolors(i,:), 'filled');
end
patch('faces',T,'vertices',X,'facecolor','white','facealpha',1,'edgecolor','none');


% metaData = getMetaGraph(mesh, labels);
% computePolycubeEmbedding(metaData)

%% visualize
% figure; clf; hold all; axis equal; rotate3d on; fa = 1; 
% quiver3(triverts(:,1),triverts(:,2),tiverts(:,3),trivertNormals(:,1),trivertNormals(:,2),trivertNormals(:,3));
% patch('faces',T,'vertices',X,'facecolor','red','facealpha',fa,'edgecolor','none');
% patch('faces',T(labels==1,:),'vertices',X,'facecolor','red','facealpha',fa,'edgecolor','none');
% patch('faces',T(labels==2,:),'vertices',X,'facecolor','green','facealpha',fa,'edgecolor','none');
% patch('faces',T(labels==3,:),'vertices',X,'facecolor','blue','facealpha',fa,'edgecolor','none');
% patch('faces',T(labels==4,:),'vertices',X,'facecolor','red','facealpha',fa,'edgecolor','none');
% patch('faces',T(labels==5,:),'vertices',X,'facecolor','green','facealpha',fa,'edgecolor','none');
% patch('faces',T(labels==6,:),'vertices',X,'facecolor','blue','facealpha',fa,'edgecolor','none');
% quiver3(mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3),mesh.vertNormals(:,1),mesh.vertNormals(:,2),mesh.vertNormals(:,3))
