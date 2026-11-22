% Widths and angular changes use profiles from the preceding sibling.
% Counts are recovered from the shared cylinder geometry, not saved twice.
measure(fileparts(mfilename('fullpath')));

function measure(here)
root=fileparts(fileparts(here)); fieldDir=fullfile(root,'08_polar&domain&POP');
addpath(fieldDir); addpath(here);
m=load(fullfile(root,'06_final_coordinates','Model_output_02_Final_model.mat'),'Blat_final');
f=load(fullfile(fieldDir,'output','polar_field.mat'),'centers_A','cell_A');
a=load(fullfile(root,'10_Polar_order_and_spatial_correlations','Domain_adjacency','output','adjacent_pairs.mat'),'pairs');
p=load(fullfile(fileparts(here),'Core_to_core_profiles','output','profiles.mat'),'profiles');
xyz=pmn.lattice(m.Blat_final); n=size(a.pairs,1);
width_A=nan(n,3); signed_net_deg=nan(n,2); total_variation_deg=nan(n,2);
scalar_method=zeros(n,2,'uint8');
grid=false(max(xyz,[],1)); grid(sub2ind(size(grid),xyz(:,1),xyz(:,2),xyz(:,3)))=true;
for i=1:n
    row=p.profiles{i}; ids=a.pairs(i,:);
    [qr,~,~,count]=pmn.tube(xyz,f.centers_A(ids,:)/f.cell_A,f.cell_A);
    assert(isequal(qr,row(:,1)) && isequal(count>0,isfinite(row(:,2))));
    keep=count>0; x=qr(keep); y=row(keep,2:5);
    c=f.centers_A(ids,:); L=norm(diff(c,1,1));
    window=valley_window(x,y(:,3:4),count(keep),L,f.cell_A);
    hasVoid=path_has_void(c/f.cell_A,grid);
    for j=1:2
        weights=count;
        % Study-specific rejection of the final POP observation in D5-D23.
        if j==1 && isequal(ids,[5 23]), weights(find(keep,1,'last'))=0; end
        [b,scalar_method(i,j)]=measure_depression(qr,row(:,j+1),weights,window,hasVoid,2);
        width_A(i,j)=diff(b);
    end
    curves=fit_profile(x,y(:,3:4),count(keep),true);
    [width_A(i,3),signed_net_deg(i,:),total_variation_deg(i,:)]=angular_change(curves,x(1),x(end));
end
width_columns={'POP','polarization_magnitude','angular_reorientation'};
out=fullfile(here,'output'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,'transition_widths.mat'),'width_A','width_columns','scalar_method','signed_net_deg','total_variation_deg');
fprintf('Width computation completed for %d pairs.\n',n);
end
