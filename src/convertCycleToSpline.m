function convertCycleToSpline(cycleXY)
    padsize = 10;
    assert(size(cycleXY,1) >= padsize);
    
    figure; hold all; rotate3d on; axis equal;
    plot(cycleXY(:,1), cycleXY(:,2),'b.-','linewidth',2,'markersize',10)

    cycleXYPadded = [cycleXY(end-padsize+1:end,:); cycleXY; cycleXY(1:padsize,:)];
    
    ts = 1:size(cycleXYPadded,1);
    denseTs = linspace(1,max(ts),10000);
    xspline = spline(ts, cycleXYPadded(:,1), denseTs);
    yspline = spline(ts, cycleXYPadded(:,2), denseTs);

    plot(xspline, yspline, 'r')
    
    
end

n=100;
ts = 1:n;
cycleXYPadded = ts*0; cycleXYPadded(1) = 100;
denseTs = linspace(1,n,n*100);
pp = spline(ts, cycleXYPadded);
pp.coefs
xspline = spline(ts, cycleXYPadded, denseTs);



