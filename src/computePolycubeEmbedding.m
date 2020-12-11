function mVertexX = computePolycubeEmbedding(metaData)
    load labelDirs.mat;
    
    cvx_begin
        cvx_solver mosek
        variable mVertexX(metaData.nmVertices,3);
%         maximize sum(sum(mVertexX(2:end,:)));
        subject to
        mVertexX(1,:) == 0;
        for i=1:size(metaData.nmFaces,1)
            label = metaData.metaFaces2labels(i);
            fixedcoord = find(labelDirs(:,label));
            
            vsi{i} = mVertexX(metaData.metaFaces{i}, fixedcoord);
            for j=2:numel(vsi{i})
                vsi{i}(1) == vsi{i}(j);
            end
        end
        norm(mVertexX,'fro') <= 1;        
%         for i=1:size(metaData.nmEdges,1)
%             el(i) = norm(mVertexX(metaData.metaEdges(i,:),:),'fro');
%         end        
    cvx_end
    
    figure; hold all; rotate3d on; axis equal;
    scatter3(mVertexX(:,1),mVertexX(:,2),mVertexX(:,3),20,'r','filled')
    
    
end