% triangle meshes
function mesh = processXTTri(X,T)

lite = 1;
name = 'mesh';
numEigs = 0;

if nargin == 0
    [X,T]=readOBJ("../../..\RFGParametrizationBenchmark/dragonstand_recon100K.obj"); name = 'dragonstand';
%     [X,T]=readOBJ("../../..\RFGParametrizationBenchmark/amphora.obj"); name = 'amphora';
end

mesh = [];
mesh.vertices = X;
mesh.triangles = double(T);
mesh.name = name;

% [mesh.cotLaplacian,mesh.areaWeights] = cotLaplacian(X,T);

% Change to negative cot Laplacian and rescale to area = 1
% mesh.areaWeights = mesh.areaWeights / sum(mesh.areaWeights);
% mesh.cotLaplacian = -1*mesh.cotLaplacian;

mesh.numVertices = size(X,1);
mesh.numTriangles = size(T,1);

% evec = [mesh.triangles(:,1) mesh.triangles(:,2)];
% evec = [evec; mesh.triangles(:,2) mesh.triangles(:,1)];
% evec = [evec; mesh.triangles(:,1) mesh.triangles(:,3)];
% evec = [evec; mesh.triangles(:,3) mesh.triangles(:,1)];
% evec = [evec; mesh.triangles(:,2) mesh.triangles(:,3)];
% evec = [evec; mesh.triangles(:,3) mesh.triangles(:,2)];
% evec = unique(evec,'rows');

% [~, idx] = unique(cotweights(:,[1 2]),'rows');
% cotweights=cotweights(idx,3);
% cotweights=cotweights(evec(:,1) < evec(:,2));

% mesh.edges = evec(evec(:,1) < evec(:,2),:);
% mesh.numEdges = size(mesh.edges,1);

% Compute LB eigenstuff
% areaMatrix = sparse(1:mesh.numVertices,1:mesh.numVertices,mesh.areaWeights);

if numEigs > 0
    mesh.cotLaplacian = cotmatrix(X,T);
    [evecs, evals] = eigs(mesh.cotLaplacian, speye(size(X,1)), max(numEigs,1), -1e-5);
    evals = diag(evals);
    mesh.laplaceBasis = evecs;
    mesh.eigenvalues = evals;
end



normalf = cross( mesh.vertices(mesh.triangles(:,2),:)'-mesh.vertices(mesh.triangles(:,1),:)', ...
                 mesh.vertices(mesh.triangles(:,3),:)'-mesh.vertices(mesh.triangles(:,1),:)' );
d = sqrt( sum(normalf.^2,1) ); d(d<eps)=1;
mesh.faceNormals = (normalf ./ repmat( d, 3,1 ))';

% ------------------------------------------------- paul adds stuff...
data = mesh;

%% extract surface mesh data.
rawedges = data.triangles(:,[1 2 2 3 3 1])';
edges = reshape(rawedges,2,[])';
data.edges = unique(sort(edges,2),'rows');
data.numEdges = size(data.edges,1);
T2V = sparse(repmat(1:data.numTriangles,3,1),data.triangles',ones(numel(data.triangles),1),data.numTriangles,data.numVertices);
V2E = sparse(repmat(1:data.numEdges,2,1),data.edges',ones(numel(data.edges),1),data.numEdges,data.numVertices)';
T2E = T2V*V2E; T2E = T2E == 2;
if any(sum(T2E)>2)
    nonManifoldEdges = find(sum(T2E)>2);
    error('Non Manifold Triangles Found!!!');
end
isBoundaryEdge = sum(T2E)==1;
[~, ia] = sort(isBoundaryEdge);
data.edges = data.edges(ia,:);
T2V = sparse(repmat(1:data.numTriangles,3,1),data.triangles',ones(numel(data.triangles),1),data.numTriangles,data.numVertices);
V2E = sparse(repmat(1:data.numEdges,2,1),data.edges',ones(numel(data.edges),1),data.numEdges,data.numVertices)';
T2E = T2V*V2E; T2E = T2E == 2;
isBoundaryEdge = sum(T2E)==1;
data.numInteriorEdges = sum(~isBoundaryEdge);
data.numBoundaryEdges = sum(isBoundaryEdge);
[ii1, jj1] = find(T2E');
[ii2, jj2] = find(T2E(:,1:data.numInteriorEdges));
data.triangles2edges = sort(reshape(ii1,3,[])',2);
data.edges2triangles = sort(reshape(ii2,2,[])',2); % data.numInteriorEdges x 2
data.isBoundaryEdge = isBoundaryEdge;

n=data.numTriangles;
adj = sparse(repmat((1:n)',3,1),data.triangles2edges(:),ones(3*n,1)); % tris x edges
[ii,jj] = find(adj(:,(data.numInteriorEdges+1):end));
data.edges2triangles = [data.edges2triangles; ii ii*0];
%data.boundaryEdges2triangle = ii;

%{
% deprecated edge loading code. not explicit boundary/nonmanifold-ness detection.
data.edges2 = [data.triangles(:,[1 2]) (1:data.numTriangles)'; data.triangles(:,[2 3]) (1:data.numTriangles)'; data.triangles(:,[3 1]) (1:data.numTriangles)'];
lr = data.edges2(data.edges2(:,1) > data.edges2(:,2),:);
rl = data.edges2(data.edges2(:,1) < data.edges2(:,2),[2 1 3]);
lr = sortrows(lr); rl = sortrows(rl);
data.edges2 = lr(:,1:2); 
if size(lr,1)~=size(rl,1); error('Mesh has boundary!!! Not handled still :('); end;
data.edges2triangles2 = sort([lr(:,3) rl(:,3)],2);
if data.numEdges ~= size(data.edges,1); error('Non Manifold Triangles Found!!!'); end; % run assert code below to find.
res = sortrows([lr(:,3) (1:data.numEdges)'; rl(:,3) (1:data.numEdges)']);
data.triangles2edges2 = sort(reshape(res(:,2),3,[])',2);
%}

% assert code to find non manifold triangles.
%{
for i = 1:mesh.numTriangles
%     assert(numel(find(sum(ismember(mesh.triangles, mesh.triangles(i,[1 2 3])),2)==3))==1)
    
    assert(numel(find(sum(ismember(mesh.triangles, mesh.triangles(i,[1 2])),2)==2))==2)
    assert(numel(find(sum(ismember(mesh.triangles, mesh.triangles(i,[3 2])),2)==2))==2)
    assert(numel(find(sum(ismember(mesh.triangles, mesh.triangles(i,[1 3])),2)==2))==2)
end
T=mesh.triangles;hold all;
ptc = patch('Faces',T,'Vertices',X,'FaceColor','green'); rotate3d on; axis equal;
alpha(ptc,.1)
patch('Faces',T(i,:),'Vertices',X,'FaceColor','red'); 
%}

% need edges to boundary triangles
% data.tt2e = sparse(data.edges2triangles(:,1),data.edges2triangles(:,2),0*(1:data.numEdges)+1,data.numTriangles,data.numTriangles);
% data.tt2e = data.tt2e + data.tt2e';

data.triangleBarycenters = (data.vertices(data.triangles(:,1),:)+data.vertices(data.triangles(:,2),:)+data.vertices(data.triangles(:,3),:))/3;

x1 = data.vertices(data.triangles(:,1),1);
y1 = data.vertices(data.triangles(:,1),2);
z1 = data.vertices(data.triangles(:,1),3);
x2 = data.vertices(data.triangles(:,2),1);
y2 = data.vertices(data.triangles(:,2),2);
z2 = data.vertices(data.triangles(:,2),3);
x3 = data.vertices(data.triangles(:,3),1);
y3 = data.vertices(data.triangles(:,3),2);
z3 = data.vertices(data.triangles(:,3),3);

A = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2);
B = sqrt((x2-x3).^2+(y2-y3).^2+(z2-z3).^2);
C = sqrt((x3-x1).^2+(y3-y1).^2+(z3-z1).^2);
S = (A+B+C)/2;
data.triangleAreas = sqrt(S.*(S-A).*(S-B).*(S-C));
data.triangleEdgeBLength = B; % an arbitrary edge per triangle

data.edgeLengths = sqrt(sum((data.vertices(data.edges(:,1),:) - data.vertices(data.edges(:,2),:)).^2,2));
data.dualEdgeLengths = sqrt(sum((data.triangleBarycenters(data.edges2triangles(1:data.numInteriorEdges,1),:) - data.triangleBarycenters(data.edges2triangles(1:data.numInteriorEdges,2),:)).^2,2));
data.primalOverDualWeight = data.edgeLengths(1:data.numInteriorEdges)./data.dualEdgeLengths;

data.primalIncidence = sparse(repmat(1:data.numEdges,1,2), data.edges(:), [ones(data.numEdges,1); -ones(data.numEdges,1)], data.numEdges, data.numVertices);

% cotweights=1./cotweights;
% (cotweights.*data.primalIncidence)'*(cotweights.*data.primalIncidence)

% dual incidence. IExT
data.incidenceMatrix = sparse(1:data.numInteriorEdges, data.edges2triangles(1:data.numInteriorEdges,1), ones(data.numInteriorEdges,1), data.numInteriorEdges, data.numTriangles);
data.incidenceMatrix = data.incidenceMatrix - sparse(1:data.numInteriorEdges, data.edges2triangles(1:data.numInteriorEdges,2), ones(data.numInteriorEdges,1), data.numInteriorEdges, data.numTriangles);
data.dualGraphL = data.incidenceMatrix'*data.incidenceMatrix;

% data.weightedIncidenceMatrix = data.incidenceMatrix./sqrt(data.primalOverDualWeight);
% data.dualL = data.weightedIncidenceMatrix'*data.weightedIncidenceMatrix;

data.triangle2verts = T2V; %sparse(repmat(1:data.numTriangles,1,3), data.triangles(:), ones(data.numTriangles*3,1), data.numTriangles, data.numVertices);

data.vertNormals = data.triangle2verts'*(data.faceNormals.*data.triangleAreas);
data.vertNormals = data.vertNormals ./ vecnorm(data.vertNormals,2,2);

data.cotlaplacian = cotmatrix(X,T);

if(~lite)
    disp('YOU PROBABLY DONT WANT OR NEED THE FULL LOAD.');
    
    %% load edge cycles and triangle cycles
    edgeCycles = cell(data.numVertices,1);
    triCycles = cell(data.numVertices,1);
    data.VV2E = sparse(data.edges(:,1),data.edges(:,2),1:data.numEdges,data.numVertices,data.numVertices); data.VV2E = data.VV2E+data.VV2E';
    for i = 1:data.numVertices
        adjEdges = data.VV2E(i,find(data.VV2E(i,:)));
        assert(numel(adjEdges) >= 2);
        edgeCycles{i} = adjEdges(1);
        adjEdges = adjEdges(2:end);
        triangles = data.edges2triangles(edgeCycles{i},:);
        nextEdge = intersect(ARemoveB(data.triangles2edges(triangles(1),:),edgeCycles{i}(end)),adjEdges);
        prevTri = triangles(1); triCycles{i} = prevTri;
        while(numel(adjEdges)~=0)
            edgeCycles{i} = [edgeCycles{i} nextEdge];
            adjEdges = ARemoveB(adjEdges,nextEdge);
            triangles = data.edges2triangles(nextEdge,:);
            prevTri = ARemoveB(data.edges2triangles(nextEdge,:), prevTri);
            triCycles{i} = [triCycles{i} prevTri];
            nextEdge = intersect(ARemoveB(data.triangles2edges(prevTri,:),edgeCycles{i}(end)),adjEdges);
        end
        %% orient cycle correctly.
        v1 = data.vertices(i,:);
        v2 = data.vertices(ARemoveB(data.edges(edgeCycles{i}(1),:),i),:);
        v3 = data.vertices(ARemoveB(data.edges(edgeCycles{i}(2),:),i),:);
        normal1 = cross(v2-v1,v3-v1); normal1 = normal1/norm(normal1);
        normal2 = data.faceNormals(triCycles{i}(1),:);
        if(dot(normal1,normal2)<0)
            data.triCycles{i} = fliplr(triCycles{i});
            data.edgeCycles{i} = fliplr(edgeCycles{i});
        end
    end
end



%% compute cotan weights.
%{
tic
e26v = [data.triangles(data.edges2triangles(:,1),:) data.triangles(data.edges2triangles(:,2),:)]
[ii,jj]=find(~([sum(e26v(:,1)==e26v(:,4:6),2) sum(e26v(:,2)==e26v(:,4:6),2) sum(e26v(:,3)==e26v(:,4:6),2)])); [~,perm] = sort(ii); jj = jj(perm); [J] = sub2ind(size(e26v),sort(ii),jj);
uniqueVertTri1 = e26v(J);
[ii,jj]=find(~([sum(e26v(:,4)==e26v(:,1:3),2) sum(e26v(:,5)==e26v(:,1:3),2) sum(e26v(:,6)==e26v(:,1:3),2)])); [~,perm] = sort(ii); jj = jj(perm); [J] = sub2ind(size(e26v),sort(ii),jj+3);
uniqueVertTri2 = e26v(J);
[ii,jj]=find(([sum(e26v(:,4)==e26v(:,1:3),2) sum(e26v(:,5)==e26v(:,1:3),2) sum(e26v(:,6)==e26v(:,1:3),2)])); [~,perm] = sort(ii); jj = jj(perm); [J] = sub2ind(size(e26v),sort(ii),jj+3);
sharedVerts = reshape(e26v(J),2,[])';
%e24v = [uniqueVertTri1 sharedVerts uniqueVertTri1];
cota = cot(abs(acos(dot(data.vertices(sharedVerts(:,1),:)-data.vertices(uniqueVertTri1,:),data.vertices(sharedVerts(:,2),:)-data.vertices(uniqueVertTri1,:),2))));
cotb = cot(abs(acos(dot(data.vertices(sharedVerts(:,1),:)-data.vertices(uniqueVertTri2,:),data.vertices(sharedVerts(:,2),:)-data.vertices(uniqueVertTri2,:),2))));
e2w = (cota+cotb)/2;
toc
%}

%% problem is there's many right angle equilateral back to back triangles. those have the same circumcenter, meaning dual edge length is 0.
% will use barycenters instead. :(. 
% faster cotan weights computation
% TR = triangulation(data.triangles,data.vertices);
% circents = circumcenter(TR);
% data.triangleCircumcenters = circents;
data.edgeWeights =  1./data.primalOverDualWeight;
assert(all(data.edgeWeights>=0));

% data.tri2vert = sparse(repmat([1:data.numTriangles],1,3),data.triangles(:),repmat(data.triangleAreas,3,1),data.numTriangles,data.numVertices);
% data.vertAreas = sum(data.tri2vert)/3;

% data.symmetricCotanL = data.primalIncidence'*diag(sparse(data.edgeWeights))*data.primalIncidence;
% data.symmetricDualCotanL = data.incidenceMatrix'*diag(sparse(1./data.edgeWeights))*data.incidenceMatrix;
% data.symmetricDualCotanL = data.incidenceMatrix'*data.incidenceMatrix;
% data.asymmetricDualCotanL = (1./data.triangleAreas).* data.symmetricDualCotanL;
% data.asymmetricCotanL = (1./data.vertAreas).* data.symmetricCotanL;

% tri2AdjTris2 = [data.edges2triangles(data.triangles2edges(:,1),:) data.edges2triangles(data.triangles2edges(:,2),:) data.edges2triangles(data.triangles2edges(:,3),:)];
% tri2AdjTris = tri2AdjTris2; assert(all(sum(tri2AdjTris == (1:data.numTriangles)',2)==3));
% ii = find(tri2AdjTris' == (1:data.numTriangles));
% tri2AdjTris=tri2AdjTris'; tri2AdjTris(ii)=[];
% data.triangle2AdjacentTriangles = reshape(tri2AdjTris,3,[])';


%% compute triXtri2edge
% adj = sparse(repmat((1:n)',3,1),data.triangles2edges(:),ones(3*n,1)); % tris x edges
%adj_edgelabeled = sparse(repmat((1:n)',3,1),data.triangles2edges(:),data.triangles2edges(:)); % tris x edges
adj_trilabeled = sparse(repmat((1:n)',3,1),data.triangles2edges(:),repmat([1:data.numTriangles]',3,1)); % tris x edges
% compute tri x tri 2 tet
edgeXedge2tri = (adj'*adj_trilabeled);
edgeXedge2tri(1:data.numEdges+1:end)=0;
data.edgeXedge2tri = edgeXedge2tri;
%{
for i=1:100000
    rtri = randsample(data.numTriangles,1);
    redges = randsample(data.triangles2edges(rtri,:),2);
    assert(edgexedge2tri(redges(1),redges(2))==rtri);
end
%}
% compute edges to triangles. one tri will be 0 for boundary edges
data.triXtri2edge = sparse(data.edges2triangles(1:data.numInteriorEdges,1),data.edges2triangles(1:data.numInteriorEdges,2),1:data.numInteriorEdges,data.numTriangles,data.numTriangles);
data.triXtri2edge = data.triXtri2edge + data.triXtri2edge';
data.triXedge2tri = sparse(data.edges2triangles(1:data.numInteriorEdges,1),1:data.numInteriorEdges,data.edges2triangles(1:data.numInteriorEdges,2),data.numTriangles,data.numEdges) +...
    sparse(data.edges2triangles(1:data.numInteriorEdges,2),1:data.numInteriorEdges,data.edges2triangles(1:data.numInteriorEdges,1),data.numTriangles,data.numEdges);
%{
% checked data.triXedge2tri
for i=1:100000
    redge = randsample(data.numEdges,1);
    if redge <= data.numInteriorEdges
        tris = data.edges2triangles(redge,:);
        assert(data.triXedge2tri(tris(1), redge)==tris(2));
        assert(data.triXedge2tri(tris(2), redge)==tris(1));
        assert(norm(data.triXedge2tri(ARemoveB(1:data.numTriangles,tris), redge))==0);
    else
        tri = data.edges2triangles(redge,1);
        assert(data.triXedge2tri(tri, redge)==0);
    end
end
%}

%{
bind = sum(adj);
[bind1,perm]=sort(bind);
adjp = adj(:,perm);
boundarypart = adjp(:,1:(find(bind1==2,1)-1));
intpart = adjp(:,find(bind1==2,1):end);
[first_boundaryTri2Tet, ~]=find(boundarypart); 
first_boundaryTri2Tet_padded = [first_boundaryTri2Tet zeros(numel(first_boundaryTri2Tet),1)];
[ii,~] = find(intpart); first_intTri2Tet = reshape(ii,2,[])';
first_tri2tet = [first_boundaryTri2Tet_padded; first_intTri2Tet];
triangles2tets = first_tri2tet*0; 
triangles2tets(perm,:) = first_tri2tet;
data.triangles2tets = triangles2tets; % [tet1 0;...] for boundary tris. [tet1 tet2] for interior tris.

tetxtet_trilabeled = adj_trilabeled*adj';
tetxtet_trilabeled(1:data.numTetrahedra+1:end)=0;
[ii,jj,kk]=find(tetxtet_trilabeled);
data.tetXtri2tet = sparse(ii,kk,jj,data.numTetrahedra,data.numTriangles);
data.tetXtri2tet(:,data.isBoundaryTriangle)=0;
%}




v1 = data.vertices(data.triangles(:,1),:);
v2 = data.vertices(data.triangles(:,2),:);
v3 = data.vertices(data.triangles(:,3),:);

function altitudes = alt(v1,v2,v3)
    e12 = v1-v2;
    e32 = v3-v2;
    e32n = e32./vecnorm(e32,2,2);
    
    altitudes = e12 - dot(e12,e32n,2).*e32n;
end

data.TriVertAltitude1 = alt(v1,v2,v3);
data.TriVertAltitude2 = alt(v2,v3,v1);
data.TriVertAltitude3 = alt(v3,v1,v2);


























mesh = data;

end

