% orients curves so the outside is the left.
function [res_digraph, cycles] = orientCycles(cycles, E, U, mesh, J)
    X = mesh.vertices;
    T = mesh.triangles;
    res_digraph = digraph();
    for j=1:numel(cycles)
%         eind = edgeIndices{j}(1);
        einds = cycles{j}(1:2);
        [foundflag, eind] = ismember(sort(einds,2), sort(E,2), 'rows');
        assert(foundflag);
        triInd = J(eind);
        
        edgeVerts = U(einds, :);
        triVerts = X(T(triInd,:),:);
        
%         plot3(edgeVerts(:,1),edgeVerts(:,2),edgeVerts(:,3),'r')
%         plot3(triVerts(:,1),triVerts(:,2),triVerts(:,3),'g')
        
        normalDir = mesh.faceNormals(triInd,:);
        eDir = edgeVerts(2,:)-edgeVerts(1,:);
        eDirRot = [-eDir(2) eDir(1) 0]; % rotate 90 deg ccw
        if dot(eDirRot(1:2), normalDir(1:2)) < 0 
            cycles{j} = fliplr(cycles{j});
        end
        
        res_digraph = addedge(res_digraph, cycles{j}, circshift(cycles{j},-1));
        % [foundflags, edgeIndex] = ismember(sort([cycles{j}; circshift(cycles{j},-1)]',2), sort(E,2), 'rows');
        % assert(all(foundflags));
        % edgeIndices{j} = edgeIndex;
    end


end