files = dir('*.vtk');
for i=1:numel(files)
    fname = files(i).name;
    [~,name,ext] = fileparts(fname);
    fid = fopen(fname);
    [X,T] = loadVTKTET(fid);
    
    data = getTetDataRT(T,X,1,0,0);
    BTs = data.triangles(data.boundaryTriangles,:);
    [BX,BT] = minimizeMesh(X,BTs);
    writeOBJ([name '.obj'],BX,BT(:,[1 3 2]));
end