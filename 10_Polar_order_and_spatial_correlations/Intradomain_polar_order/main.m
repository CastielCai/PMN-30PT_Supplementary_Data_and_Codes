% Use accepted domain/core labels (08), the common POP (08), and sites (06).
domainOrder(fileparts(mfilename('fullpath')));

function domainOrder(here)
root = fileparts(fileparts(here)); fieldDir=fullfile(root,'08_polar&domain&POP'); addpath(fieldDir);
m=load(fullfile(root,'06_final_coordinates','Model_output_02_Final_model.mat'),'Blat_final','B_tr_rep');
f=load(fullfile(fieldDir,'output','polar_field.mat'),'POP','cell_A','pixel_A');
a=load(fullfile(fieldDir,'input','reviewed_domains.mat'),'domain_id','core_id');
xyz=pmn.lattice(m.Blat_final); nDomains=double(max(a.domain_id));
size_cells=zeros(nDomains,1); radial_POP=cell(nDomains,1);
for d=1:nDomains
    members=find(a.domain_id==d); cores=find(a.core_id==d);
    size_cells(d)=numel(members);
    distance=min(pdist2(xyz(members,:),xyz(cores,:)),[],2);
    shell=ceil(distance);
    bins=unique(shell);
    av=arrayfun(@(b) mean(f.POP(members(shell==b))),bins);
    radial_POP{d}=[bins*f.cell_A,av];
end
% Fig. 2E uses the physical volume defined by the fitted lattice vectors.
position_A=(double(m.B_tr_rep)+2)*f.pixel_A/3;
affine=[xyz,ones(size(xyz,1),1)]\position_A;
cell_volume_A3=abs(det(affine(1:3,:)));
equivalent_diameter_nm=(6*size_cells*cell_volume_A3/pi).^(1/3)/10;
out=fullfile(here,'output'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,'domain_order.mat'),'size_cells','equivalent_diameter_nm','radial_POP');
fprintf('Intradomain order: %d domains.\n',nDomains);
end
