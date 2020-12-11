% Given X, T of a tet mesh and vertex weights, computes Lw, the weighted laplacian. 
% weights are piecewise linearly interpolated, and integrated against by grad(f) of a hypothetical function f.
% L is the unweighted laplacian. M is its mass matrix. Lij are the 1form metric (hodge star).
function [Lw, L, M, Lij] = ComputeWeightedLaplacian(X, T, weights, data, tetra)
    assert(numel(weights)==size(X,1))
    
    if ~exist('tetra','var')
        tetra = triangulation(T, X);
    end
    if ~exist('data','var')
        data = processXTtetra(T, X, 1);
    end
    
    [L, M, Lij] = GeometricPrimalLM(tetra);
    
    Lw = data.linearVertsToConstantTetsOp'* diag(sparse(repelem(data.tetVolumes.*mean(weights(data.tetrahedra),2),3,1))) * data.linearVertsToConstantTetsOp;
    
    %% debugging only
    if nargin == 0
        f = randn(data.numVertices,1);
        D1 = vecnorm(reshape(data.linearVertsToConstantTetsOp*f,3,[])',2,2).^2'*(data.tetVolumes.*mean(weights(data.tetrahedra),2));
        D2 = f'*Lw*f;
        zzz = data.primalIncidenceMatrix' * diag(sparse(Lij)) * data.primalIncidenceMatrix - L;
        assert(abs(D1-D2) < .00001);
        assert(all(zzz(:)==0));
    end
end