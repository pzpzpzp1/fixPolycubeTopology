function cycleData = markCycles(cycle, U, Visualize, useInflections)
    if nargin==0
%         U = userinput_points('a.txt', zeros(10)); 
        U = loadCachedShape('a.txt');
        
        U = U(1:end-1,:);
        U(:,3)=0;
        cycle = [1:size(U,1)];        
        Visualize = 1;
    end
    
    cycleXY = U(cycle,:); N = numel(cycle);
    
    % HACK: make cycleXY generic. 
    %{
    eps = 1e-3;
    eforw = circshift(cycleXY,-1) - cycleXY;
    edgeNormalsForw = [-eforw(:,2) eforw(:,1) 0*eforw(:,1)]; edgeNormalsForw = edgeNormalsForw./vecnorm(edgeNormalsForw,2,2);
    incid = sparse(N, N);    incid(1:N+1:end) = 1;    incid(N+1:N+1:end) = -1;    incid(end,1) = -1; D = diag(sparse(vecnorm(eforw,2,2))); Lap = incid'*D*incid; evs = eigs(Lap,2,'sm'); evs(2); Lap = Lap / evs(2); dt = .0001;
    pert = (speye(N) - dt*Lap)\edgeNormalsForw;
    cycleXY = cycleXY + eps * pert + eps*eps*randn(size(cycleXY));
    %{
    clf; hold all; axis equal;
    plot3(cycleXY(:,1),cycleXY(:,2),cycleXY(:,3),'.-')
    quiver3(edgeCentersForw(:,1), edgeCentersForw(:,2), edgeCentersForw(:,3), pert(:,1), pert(:,2), pert(:,3));
    %}
    %}
     
    eforw = circshift(cycleXY,-1) - cycleXY;
    eback = cycleXY - circshift(cycleXY,1);
    
    normalEps = 0; %.017; %1e-6;
    ef = eforw./vecnorm(eforw,2,2);
    eb = eback./vecnorm(eback,2,2);
    isConvex = cross(ef, eb)*[0;0;1] > normalEps;
    isConcave = cross(ef, eb)*[0;0;1] < -normalEps;
    isFlat = ~isConcave & ~isConvex;
    
    vertexconvexity(isConvex)=1;
    vertexconvexity(isConcave)=-1;
    vertexconvexity(isFlat)=0;
    
    edgeCentersForw = (circshift(cycleXY,-1) + cycleXY)/2;
    edgeNormalsForw = [-eforw(:,2) eforw(:,1) 0*eforw(:,1)]; edgeNormalsForw = edgeNormalsForw./vecnorm(edgeNormalsForw,2,2);
    
    
    
    xc = cycleXY(:,1)
    xf = circshift(xc,-1);
    xb = circshift(xc,1);
    yc = cycleXY(:,2)
    yf = circshift(yc,-1);
    yb = circshift(yc,1);
    
    isXmax = xc > xf & xc > xb;
    isXmin = xc < xf & xc < xb;
    isYmax = yc > yf & yc > yb;
    isYmin = yc < yf & yc < yb;
    isXYminmax = isXmax | isXmin | isYmax | isYmin;
    
    % get inflection points
    prevVal = vertexconvexity(1);
    inflectionEdges = 0*vertexconvexity;
    inflectionEdgeTransition = zeros(numel(vertexconvexity),2);
    if useInflections
        for i=2:numel(vertexconvexity)+1
            if i==numel(vertexconvexity)+1
                val = vertexconvexity(1);
            else
                val = vertexconvexity(i);
            end

            if val == 0
                continue;
            end

            if val ~= prevVal
                inflectionEdges(i-1) = 1;
                inflectionEdgeTransition(i-1,:) = [prevVal val];
            end

            prevVal = val;
        end
    end
    
    if Visualize
        figure; hold all; rotate3d on; axis equal;
        sz = 50;
        scatter3(cycleXY(isXmax,1),cycleXY(isXmax,2),cycleXY(isXmax,3),sz,'b','filled')
        scatter3(cycleXY(isXmin,1),cycleXY(isXmin,2),cycleXY(isXmin,3),sz,'b','filled')
        scatter3(cycleXY(isYmax,1),cycleXY(isYmax,2),cycleXY(isYmax,3),sz,'b','filled')
        scatter3(cycleXY(isYmin,1),cycleXY(isYmin,2),cycleXY(isYmin,3),sz,'b','filled')
        scatter3(edgeCentersForw(inflectionEdges==1,1), edgeCentersForw(inflectionEdges==1,2), edgeCentersForw(inflectionEdges==1,3),sz,'c','filled');
        plot3(cycleXY(:,1),cycleXY(:,2),cycleXY(:,3),'r.-')
    end
    
    MetaVerts = {};
    for i=1:N
        if isXYminmax(i)
            MetaVerts{end+1}.pos = cycleXY(i,:);
            MetaVerts{end}.isVert = 1;
            MetaVerts{end}.convexity = vertexconvexity(i);
        end
        if inflectionEdges(i)==1
            MetaVerts{end+1}.pos = edgeCentersForw(i,:);
            MetaVerts{end}.isVert = 0;
            MetaVerts{end}.convexity = 0;
            MetaVerts{end}.transition = inflectionEdgeTransition(i,:);
        end
    end
    MN = numel(MetaVerts);
    MetaEdges = cell(MN,1);
    for i=1:MN
        v1 = MetaVerts{i};
        v2 = MetaVerts{wrapN(i+1,MN)};
        mE.v1 = v1;
        mE.v2 = v2;
        if v1.isVert && v2.isVert
            mE.convexity = v1.convexity;
            if useInflections
                assert(v1.convexity == v2.convexity);
                assert(v1.convexity ~= 0);
            end
        elseif ~v1.isVert
            mE.convexity = v1.transition(2);
        elseif ~v2.isVert
            mE.convexity = v2.transition(1);
        else
            assert('unhandled');
        end
        if ~useInflections
            mE.convexity = 1;
        end
        metaEdges{i} = mE;
    end
    
    % create final cycleData
    cycleData = zeros(MN*2,3);
    for i=1:MN
        v1 = MetaVerts{i};
        v2 = MetaVerts{wrapN(i+1,MN)};
        
        edirQuad = sign(v2.pos-v1.pos); edirQuad = edirQuad(1:2);
        if all(edirQuad == [1 -1]) || all(edirQuad == [-1 1])
            if metaEdges{i}.convexity == 1
                newVert = [v2.pos(1) v1.pos(2)];
            else
                newVert = [v1.pos(1) v2.pos(2)];
            end
        elseif all(edirQuad == [1 1]) || all(edirQuad == [-1 -1])
            if metaEdges{i}.convexity == 1
                newVert = [v1.pos(1) v2.pos(2)];
            else
                newVert = [v2.pos(1) v1.pos(2)];
            end
        else
            error('uinhandled');
        end
            
        cycleData(2*i-1,:) = metaEdges{i}.v1.pos;
        cycleData(2*i,:) = [newVert metaEdges{i}.v1.pos(3)];
    end
    
    if Visualize
        figure; hold all; rotate3d on; axis equal;
        plot3(cycleXY(:,1),cycleXY(:,2),cycleXY(:,3),'k.-')
        for i=1:MN
            mE = metaEdges{i};
            v1 = MetaVerts{i};
            v2 = MetaVerts{wrapN(i+1,MN)};
            vs = [v1.pos; v2.pos];
            if mE.convexity == 1
                plot3(vs(:,1),vs(:,2),vs(:,3),'g','linewidth',1.1)
            else
                plot3(vs(:,1),vs(:,2),vs(:,3),'r','linewidth',1.1)
            end
        end
        plot3(cycleData(:,1),cycleData(:,2),cycleData(:,3),'b.-')
    end
    
%     %quiver3(edgeCentersForw(:,1), edgeCentersForw(:,2), edgeCentersForw(:,3), edgeNormalsForw(:,1), edgeNormalsForw(:,2), edgeNormalsForw(:,3));
%     plot3(cycleXY(isConvex,1),cycleXY(isConvex,2),cycleXY(isConvex,3),'g.','markersize',10)
%     plot3(cycleXY(isConcave,1),cycleXY(isConcave,2),cycleXY(isConcave,3),'y.','markersize',10)
%     plot3(cycleXY(isFlat,1),cycleXY(isFlat,2),cycleXY(isFlat,3),'m.','markersize',10)
    
end