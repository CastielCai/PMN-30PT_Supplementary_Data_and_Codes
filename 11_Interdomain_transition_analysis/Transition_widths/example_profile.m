function e=example_profile(domainPair)
% Reconstruct the plotted quantities for an observed (method 1) example.
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
field=fullfile(root,'08_polar&domain&POP'); addpath(field);
m=load(fullfile(root,'06_final_coordinates','Model_output_02_Final_model.mat'),'Blat_final');
f=load(fullfile(field,'output','polar_field.mat'),'centers_A','cell_A');
a=load(fullfile(root,'10_Polar_order_and_spatial_correlations','Domain_adjacency','output','adjacent_pairs.mat'),'pairs');
p=load(fullfile(fileparts(here),'Core_to_core_profiles','output','profiles.mat'),'profiles');
w=load(fullfile(here,'output','transition_widths.mat'),'width_A','scalar_method');
[found,i]=ismember(domainPair,a.pairs,'rows');
assert(found && all(w.scalar_method(i,:)==1),'Choose a pair with two observed half-depth widths.');
row=p.profiles{i}; c=f.centers_A(domainPair,:)/f.cell_A; xyz=pmn.lattice(m.Blat_final);
[~,~,~,count]=pmn.tube(xyz,c,f.cell_A);
sel=count>0; x=row(sel,1); y=row(sel,2:5); n=count(sel);
window=valley_window(x,y(:,3:4),n,norm(diff(c))*f.cell_A,f.cell_A);
e.qr_A=x; e.observed=y; e.bounds_A=zeros(3,2); e.reduction=zeros(numel(x),2);
for j=1:2
    ww=count;
    if j==1 && isequal(domainPair,[5 23]), ww(find(sel,1,'last'))=0; end
    [e.bounds_A(j,:),mode]=measure_depression(row(:,1),row(:,j+1),ww,window,false,2);
    assert(mode==1);
    valid=ww(sel)>0; pos=x(valid); obs=y(valid,j);
    base=interp1(pos([1 end]),obs([1 end]),x,'linear','extrap');
    dep=1-y(:,j)./base; height=max(dep(x>=window(1) & x<=window(2) & valid));
    e.baseline(:,j)=base; e.half_level(:,j)=base*(1-height/2);
    e.reduction(:,j)=100*dep/height;
end
curve=fit_profile(x,y(:,3:4),n,true);
[~,~,~,e.bounds_A(3,:),cuts,mass]=angular_change(curve,x(1),x(end));
e.grid_A=linspace(x(1),x(end),601)';
e.angle_fit=[ppval(curve{1},e.grid_A),ppval(curve{2},e.grid_A)];
k=discretize(e.grid_A,cuts); start=cuts(k);
advance=abs(e.angle_fit-[ppval(curve{1},start),ppval(curve{2},start)]);
e.cumulative=100*(mass(k)+sum(advance,2))/mass(end);
assert(max(abs(diff(e.bounds_A,1,2)'-w.width_A(i,:)))<1e-7);
end
