function gap=path_has_void(c,grid)
% Traverse the grid cells crossed by the segment joining both core centers.
delta=c(2,:)-c(1,:); stops=[0;1]; dims=size(grid);
for axis=1:3
    if abs(delta(axis))>1e-12
        t=((.5:(dims(axis)+.5))-c(1,axis))/delta(axis);
        stops=[stops;t(t>0 & t<1)']; %#ok<AGROW>
    end
end
stops=unique(stops); mid=(stops(1:end-1)+stops(2:end))/2;
v=floor(c(1,:)+mid.*delta+.5);
gap=any(~grid(sub2ind(dims,v(:,1),v(:,2),v(:,3))));
end
