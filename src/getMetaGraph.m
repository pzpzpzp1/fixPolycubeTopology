function [metaData] = getMetaGraph(mesh, labels)
    
    sameLabels = labels(mesh.edges2triangles(:,1))==labels(mesh.edges2triangles(:,2));
    
    g = addedge(graph, mesh.edges2triangles(sameLabels,1), mesh.edges2triangles(sameLabels,2));
    [faces2metaFaces, binSizes] = conncomp(g); faces2metaFaces = faces2metaFaces';
    metaData.faces2metaFaces = faces2metaFaces;
    metaData.faces2labels = labels(:);
    
    nmFaces = size(binSizes,2);
    metaData.nmFaces = nmFaces;
    for i=1:nmFaces
        metaFaces2Faces{i} = find(faces2metaFaces == i);
        metaFaces2labels(i,1) = labels(metaFaces2Faces{i}(1));
    end
    metaData.metaFaces2Faces = metaFaces2Faces;
    metaData.metaFaces2labels = metaFaces2labels;
    
    cgraph = addedge(graph, mesh.edges(~sameLabels,1), mesh.edges(~sameLabels,2));
    metaVertices2Vertices = find(cgraph.degree ~=2 & cgraph.degree ~=0);
    nmVertices = size(metaVertices2Vertices,1);
    metaData.nmVertices = nmVertices;
    metaData.metaVertices2Vertices = metaVertices2Vertices;
    vertices2metaVertices = zeros(mesh.numVertices,1);
    vertices2metaVertices(metaVertices2Vertices) = 1:numel(metaVertices2Vertices);
    metaData.vertices2metaVertices = vertices2metaVertices;
    metaData.cgraph = cgraph;
    
    % build oriented metaFaces
    metaEdgesNU = [];
    metaEdgeNUEdges = {};
    metaEdgeNUEdgeLength = {};
    for i=1:numel(metaFaces2Faces)
        tris = mesh.triangles(metaFaces2Faces{i},:);
        metaFaces2Vertices{i} = intersect(unique(tris), metaVertices2Vertices);
        
        orientedMetaFaceEdges = [tris(:,[1 2]); tris(:,[2 3]); tris(:,[3 1]);];
        [~, ia, ic] = unique(sort(orientedMetaFaceEdges,2),'rows');
        filt2 = orientedMetaFaceEdges(ia,:);
        metaFaceBoundaryEdgeIndices = find(accumarray(ic,ones(size(ic)))==1);
        
        orientedMetaFaceBoundaryEdges = filt2(metaFaceBoundaryEdgeIndices,:);
        g = addedge(digraph, orientedMetaFaceBoundaryEdges(:,1), orientedMetaFaceBoundaryEdges(:,2));
        src = orientedMetaFaceBoundaryEdges(end,2);
        snk = orientedMetaFaceBoundaryEdges(end,1);
        boundaryloop = shortestpath(g,src,snk);
        
        boundaryloop = circshift(boundaryloop, -find(ismember(boundaryloop,metaFaces2Vertices{i}),1)+1);
        
        metaFaces2Vertices{i} = boundaryloop(ismember(boundaryloop, metaFaces2Vertices{i}));
        metaFaces{i} = vertices2metaVertices(metaFaces2Vertices{i});
        
        % build metaEdges per face as well. will have duplicates.
        mF = metaVertices2Vertices(reshape(metaFaces{i},[],1));
        mE = [mF circshift(mF,1)];
        for j=1:size(mE,1)
            vertpath = shortestpath(g,mE(j,2),mE(j,1));
            metaEdgeNUEdges{end+1,1} = [vertpath(1:end-1); vertpath(2:end);]';
            
            metaEdgeNUEdgeLength{end+1,1} = sum(vecnorm(mesh.vertices(metaEdgeNUEdges{end}(:,1),:) - mesh.vertices(metaEdgeNUEdges{end}(:,2),:),2,2))
            
        end        
        metaEdgesNU = [metaEdgesNU; [mE repmat(i,size(mE,1),1)]];
        
    end    
    metaData.metaFaces = metaFaces;
    
    % build metaEdges
    [~, perm] = sortrows(sort(metaEdgesNU(:,1:2),2)); selectInds = perm(1:2:end);
    metaEdgesNU(selectInds,:)
    
    % build metaEdges
    metaEdges = zeros(0,2);
    for i = 1:numel(metaFaces)
        mF = reshape(metaFaces{i},[],1);
        mE = [mF circshift(mF,1)];
        metaEdges = [metaEdges; mE];
    end
    metaEdges = sortrows(sort(metaEdges,2)); assert(all(metaEdges(1:2:end,:)==metaEdges(2:2:end,:),'all'));
    metaEdges = metaEdges(1:2:end,:);
    metaData.metaEdges = metaEdges;
    nmEdges = size(metaEdges,1);
    metaData.nmEdges = nmEdges;
    
    % build metaEdges2Edges
    edges = sort(mesh.edges,2);
    edges2MetaEdges = zeros(mesh.numEdges,1);
    mcgraph = cgraph;
    for i=1:size(metaEdges,1)
        metaEdgeVertices = metaVertices2Vertices(sort(metaEdges(i,:)));
        path = shortestpath(mcgraph, metaEdgeVertices(1), metaEdgeVertices(2));
        mcgraph = rmedge(mcgraph, path(1), path(2));
        
        metaEdgeEdges = sort([path(1:end-1); path(2:end)]',2);
        [flags, metaEdges2Edges{i}] = ismember(metaEdgeEdges, edges, 'rows'); assert(all(flags)); 
        edges2MetaEdges(metaEdges2Edges{i}) = i;
        
        metaEdgeLength(i,1) = sum(mesh.edgeLengths(metaEdges2Edges{i}));
    end
    metaData.metaEdges2Edges = metaEdges2Edges;
    metaData.edges2MetaEdges = edges2MetaEdges;
    metaData.metaEdgeLength = metaEdgeLength;
    
    % build metaEdges2metaFaces
    metaEdges2metaFaces = 0
    
    if nargin==0
        figure; clf; hold all; axis equal; rotate3d on; fa = 1; 
        patch('faces',T(labels==1,:),'vertices',X,'facecolor','red','facealpha',fa,'edgecolor','none');
        patch('faces',T(labels==2,:),'vertices',X,'facecolor','green','facealpha',fa,'edgecolor','none');
        patch('faces',T(labels==3,:),'vertices',X,'facecolor','blue','facealpha',fa,'edgecolor','none');
        patch('faces',T(labels==4,:),'vertices',X,'facecolor','red','facealpha',fa,'edgecolor','none');
        patch('faces',T(labels==5,:),'vertices',X,'facecolor','green','facealpha',fa,'edgecolor','none');
        patch('faces',T(labels==6,:),'vertices',X,'facecolor','blue','facealpha',fa,'edgecolor','none');

        for i=1:nmVertices
            vx = mesh.vertices(metaVertices2Vertices(i),:);
            scatter3(vx(1),vx(2),vx(3),'k','filled')
        end

        patch('faces',mesh.edges(~sameLabels,[1 2 1]),'vertices',mesh.vertices,'edgecolor','cyan','linewidth',2)
        svinds = unique(mesh.edges(~sameLabels,:));
        scatter3(mesh.vertices(svinds,1),mesh.vertices(svinds,2),mesh.vertices(svinds,3),20,'filled')
    end
    
end