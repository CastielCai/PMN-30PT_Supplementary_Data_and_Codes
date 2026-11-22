function [labels,seedSites] = grow_domains(xyz,vectors,seedSites)
% Six-face, downhill growth used to propose domains before manual review.
% xyz: integer lattice sites (as returned by pmn.lattice). vectors: 08's
% smoothed three-component field. seedSites optionally supplies reviewed
% seed rows in growth-priority order. The function returns candidates in
% memory. Production analyses load the deposited, reviewed domain masks.
[edge18,coordination]=pmn.neighbors(xyz,2,Inf);
edge18=spones(edge18);
edge6=spones(pmn.neighbors(xyz,1,Inf));
edge26=spones(pmn.neighbors(xyz,3,Inf));
average=pmn.neighbors(xyz,6,Inf)*vectors;
score=sum(abs(average),2);
interior=coordination==19;
if nargin<3
    dims=max(xyz,[],1);
    site=sub2ind(dims,xyz(:,1),xyz(:,2),xyz(:,3));
    map=zeros(dims); map(site(interior))=score(interior);
    [x,y,z]=ndgrid(-2:2); rr=x.^2+y.^2+z.^2;
    localMax=map>imdilate(map,rr>0 & rr<=4);
    [ix,iy,iz]=ind2sub(dims,find(localMax & map>0));
    [~,order]=sort(vecnorm([ix iy iz]-mean(xyz,1),2,2));
    [~,seedSites]=ismember([ix(order) iy(order) iz(order)],xyz,'rows');
end
seedSites=seedSites(:);
assert(all(seedSites>=1 & seedSites<=size(xyz,1)) && numel(unique(seedSites))==numel(seedSites));
labels=zeros(size(xyz,1),1,'uint16');
n=numel(seedSites); front=cell(n,1); active=true(n,1);
for d=1:n
    sites=find(edge26(:,seedSites(d)) & interior & labels==0);
    labels(sites)=d; front{d}=sites;
end
while any(active)
    for d=find(active)'
        frontier=front{d};
        candidates=find(any(edge6(:,frontier),2) & interior & labels==0);
        accepted=false(size(candidates));
        for j=1:numel(candidates)
            k=candidates(j);
            neighbors=find(edge18(:,k)); members=neighbors(labels(neighbors)==d);
            accepted(j)=score(k)>=.1*score(seedSites(d)) && any(score(members)>score(k));
        end
        if ~any(accepted), active(d)=false; continue; end
        labels(candidates(accepted))=d;
        members=find(labels==d);
        near=vecnorm(xyz(members,:)-xyz(seedSites(d),:),2,2)<=3;
        supported=full(sum(edge26(members,members),2))>=9;
        front{d}=members(near|supported); active(d)=~isempty(front{d});
    end
end
end
