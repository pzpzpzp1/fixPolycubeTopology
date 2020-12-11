function [L, M, Lij] = GeometricPrimalLM(tetra)

nv = length(tetra.Points);

%% Compute L

E = edges(tetra);
nE = size(E, 1);
edgeStars = edgeAttachments(tetra, E);
edgeStarSizes = cellfun(@length, edgeStars);
edgeStarIdx = repelem(1:nE, edgeStarSizes)';

neighTetIdx = cell2mat(edgeStars')';
nnt = length(neighTetIdx);
neighTetVerts = tetra.ConnectivityList(neighTetIdx, :);

selfEdgeVerts = repelem(E, edgeStarSizes, 1);
assert(length(selfEdgeVerts) == nnt);

oppEdgeMask = all(bsxfun(@ne, neighTetVerts, permute(selfEdgeVerts, [1 3 2])), 3);
neighTetVertsT = neighTetVerts';
oppEdgeVerts = reshape(neighTetVertsT(oppEdgeMask'), [2 nnt])';

v0 = tetra.Points(selfEdgeVerts(:, 1), :);
v1 = tetra.Points(selfEdgeVerts(:, 2), :);
w0 = tetra.Points(oppEdgeVerts(:, 1), :);
w1 = tetra.Points(oppEdgeVerts(:, 2), :);
oppEdgeLengths = vecnorm(w1 - w0, 2, 2);
oppEdgeVecs = (w1 - w0) ./ oppEdgeLengths;
n0 = normr(cross(oppEdgeVecs, v0 - w1, 2));
n1 = normr(cross(oppEdgeVecs, v1 - w1, 2));
t1 = cross(oppEdgeVecs, n1, 2);

% Make sure we are consistently using the positively oriented angle
alpha = abs(atan2(dot(n0, t1, 2), dot(n0, n1, 2)));
oppositeCotans = cot(alpha);

Lij = oppEdgeLengths .* oppositeCotans / 6;
Lij = accumarray(edgeStarIdx, Lij);

A = sparse(E(:, 1), E(:, 2), Lij, nv, nv);
A = A + A';
L = diag(sum(A, 1)) - A;

% incidence=sparse([1:size(E,1) 1:size(E,1)],E(:),[E(:,1)*0+1;E(:,1)*0-1], nE, nv);
% L == incidence'*diag(sparse(Lij))*incidence

%% Compute Mass Matrix

Mij = abs(dot(v1 - v0, cross(w0 - v0, w1 - v0, 2), 2)) / 120;
Mij = accumarray(edgeStarIdx, Mij);

vertexStars = vertexAttachments(tetra, (1:nv)');
vertexStarSizes = cellfun(@length, vertexStars);
vertexStarIdx = repelem(1:nv, vertexStarSizes)';

vNeighTetIdx = cell2mat(vertexStars')';
vNeighTetVerts = tetra.ConnectivityList(vNeighTetIdx, :);
v0 = tetra.Points(vNeighTetVerts(:, 1), :);
v1 = tetra.Points(vNeighTetVerts(:, 2), :);
v2 = tetra.Points(vNeighTetVerts(:, 3), :);
v3 = tetra.Points(vNeighTetVerts(:, 4), :);

Mii = abs(dot(v1 - v0, cross(v2 - v0, v3 - v0, 2), 2)) / 60;
Mii = accumarray(vertexStarIdx, Mii);

M = sparse(E(:, 1), E(:, 2), Mij, nv, nv);
M = M + M' + spdiags(Mii, 0, nv, nv);

end

