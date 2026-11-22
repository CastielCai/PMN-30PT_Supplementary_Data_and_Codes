function [qr,sites,bins,counts] = tube(lattice,centers,cell_A)
% Flat-ended cylinder. All occupied sites in the sample are eligible.
axis = centers(2,:)-centers(1,:);
L = norm(axis); axis = axis/L;
q = lattice-centers(1,:);
s = sum(q.*axis,2);
r = vecnorm(q-s.*axis,2,2);
sites = find(s>=0 & s<=L & r<=8/cell_A+1e-10);
edges = (0:2/cell_A:L)';
if edges(end)<L, edges(end+1,1)=L; end
bins = discretize(s(sites),edges);
qr = cell_A*(edges(1:end-1)+edges(2:end))/2;
if nargout>3
    counts = accumarray(bins,1,[numel(qr),1]);
end
end
