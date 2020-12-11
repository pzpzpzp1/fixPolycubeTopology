clear all; close all;

%% load mesh
fname = 'testMeshes/ex1.obj';
[X, T] = readOBJ(fname);

%% process mesh
mesh = processXTtri(X,T);
labels = getLabelsFromXT(X,T);

metaData = getMetaGraph(mesh, labels);
computePolycubeEmbedding(metaData)

%% visualize
figure; clf; hold all; axis equal; rotate3d on; fa = 1; 
patch('faces',T(labels==1,:),'vertices',X,'facecolor','red','facealpha',fa,'edgecolor','none');
patch('faces',T(labels==2,:),'vertices',X,'facecolor','green','facealpha',fa,'edgecolor','none');
patch('faces',T(labels==3,:),'vertices',X,'facecolor','blue','facealpha',fa,'edgecolor','none');
patch('faces',T(labels==4,:),'vertices',X,'facecolor','red','facealpha',fa,'edgecolor','none');
patch('faces',T(labels==5,:),'vertices',X,'facecolor','green','facealpha',fa,'edgecolor','none');
patch('faces',T(labels==6,:),'vertices',X,'facecolor','blue','facealpha',fa,'edgecolor','none');
quiver3(mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3),mesh.vertNormals(:,1),mesh.vertNormals(:,2),mesh.vertNormals(:,3))
