
function [ws, interiorInds] = getBarycentricSamplingWeights(n)
    assert(n >= 5);
    persistent cachedWs cachedInteriorInds
    if isempty(cachedWs) || n > numel(cachedWs) || numel(cachedWs{n})==0
        % n would need to be 1e6 for this to become floating point imprecise.
        eps = 1e-6;
        
        u = linspace(0,1,n);
        v = linspace(0,1,n);
        [U,V] = ndgrid(u,v);
        W = 1-U-V;
        
        keep = find(W>=-eps);
        assert(numel(keep)==n*(n+1)/2); % another thresholding sanity check
        cachedWs{n} = [U(keep) V(keep) W(keep)];
        layersdeep = 3;
        cachedInteriorInds{n} = ~any(abs(cachedWs{n}) < (layersdeep-1)/(n-1) + eps,2); % 1 layer interior
    end
    ws = cachedWs{n};
    interiorInds = cachedInteriorInds{n};
end