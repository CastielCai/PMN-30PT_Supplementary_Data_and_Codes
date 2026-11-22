function [W,count] = neighbors(xyz,maxSquared,denominator)
% Row-normalized lattice neighborhood. Finite Gaussian weights use the
% original shell-rank convention at incomplete boundary neighborhoods.
dims = max(xyz,[],1);
grid = zeros(dims,'uint32');
grid(sub2ind(dims,xyz(:,1),xyz(:,2),xyz(:,3))) = uint32(1:size(xyz,1));
r = floor(sqrt(maxSquared));
[dx,dy,dz] = ndgrid(-r:r);
offsets = [dx(:) dy(:) dz(:)];
d2 = sum(offsets.^2,2);
offsets = offsets(d2<=maxSquared,:);
d2 = d2(d2<=maxSquared);
levels = unique(d2);
n = size(xyz,1);
rows = cell(size(offsets,1),1); cols = rows; vals = rows;
rank = zeros(n,1);
for level = 1:numel(levels)
    idx = find(d2==levels(level));
    present = false(n,1);
    for j = idx'
        dest = xyz+offsets(j,:);
        inside = all(dest>=1 & dest<=dims,2);
        sites = find(inside);
        dest = dest(inside,:);
        target = grid(sub2ind(dims,dest(:,1),dest(:,2),dest(:,3)));
        sites = sites(target>0);
        rows{j} = sites;
        cols{j} = double(target(target>0));
        if isinf(denominator)
            vals{j} = ones(numel(sites),1);
        else
            vals{j} = exp(-levels(rank(sites)+1)/denominator);
        end
        present(sites) = true;
    end
    rank = rank+present;
end
W = sparse(vertcat(rows{:}),vertcat(cols{:}),vertcat(vals{:}),n,n);
count = full(sum(spones(W),2));
W = spdiags(1./full(sum(W,2)),0,n,n)*W;
end
