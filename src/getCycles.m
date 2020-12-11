% Assumes each connected component of edge graph is a single loop.
% no guarantees if that's not the case.
function [cycles, edgeIndices] = getCycles(E)
    curvesG = graph();
    curvesG = addedge(curvesG, E(:,1), E(:,2));
    
    labels = conncomp(curvesG);
    nloops = max(labels);
    
    % res_digraph = digraph();
    for i=1:nloops
        seed(i) = find(labels==i,1);
        nbrs = neighbors(curvesG, seed(i));
        target(i) = nbrs(1);
        curvesG = rmedge(curvesG, seed(i), target(i));
        
        cycles{i} = shortestpath(curvesG, seed(i), target(i));
        % res_digraph = addedge(res_digraph, cycles{i}, circshift(cycles{i},1));
        
        [foundflags, edgeIndex] = ismember(sort([cycles{i}; circshift(cycles{i},-1)]',2), sort(E,2), 'rows');
        assert(all(foundflags));
        edgeIndices{i} = edgeIndex;
    end
    
    % figure; plot(res_digraph,'XData',U(:,1),'YData',U(:,2),'ZData',U(:,3));
end