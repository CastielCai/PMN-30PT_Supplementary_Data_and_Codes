% Read the unsmoothed and reference polar fields from 08; coordinates from 06.
compareKernels(fileparts(mfilename('fullpath')));

function compareKernels(here)
root=fileparts(fileparts(here));
fieldDir=fullfile(root,'08_polar&domain&POP'); addpath(fieldDir);
m=load(fullfile(root,'06_final_coordinates','Model_output_02_Final_model.mat'),'Blat_final','B_tr_rep');
f=load(fullfile(fieldDir,'output','polar_field.mat'),'vectors','raw_displacement_pixel','pixel_A');
xyz=pmn.lattice(m.Blat_final);
wide=pmn.neighbors(xyz,8,0.455*8.1)*f.raw_displacement_pixel;
POP_wide=sum(wide.*(pmn.neighbors(xyz,6,Inf)*wide),2);
pos_A=(double(m.B_tr_rep)+2)/3*f.pixel_A;
unit=[f.vectors./vecnorm(f.vectors,2,2),wide./vecnorm(wide,2,2)];
edges=(0:600)'*(200/999); N=size(xyz,1); nb=numel(edges)-1;
F=zeros(nb,2);
% Center-weighted ordered pairs. Empty shells contribute zero at that center.
for first=1:32:N
    sites=first:min(first+31,N);
    dist=pdist2(pos_A,pos_A(sites,:),'euclidean');
    bins=ceil(dist/(200/999));
    keep=dist>0 & bins<=nb;
    group=bins+(0:numel(sites)-1)*nb;
    group=group(keep);
    counts=accumarray(group,1,[nb*numel(sites),1]);
    for k=1:2
        cols=(k-1)*3+(1:3);
        dotProducts=unit(:,cols)*unit(sites,cols)';
        totals=accumarray(group,dotProducts(keep),size(counts));
        perCenter=reshape(totals./max(counts,1),nb,[]);
        F(:,k)=F(:,k)+sum(perCenter,2);
    end
end
F=F/N; r_A=edges(1:end-1);
xi_A=zeros(1,2);
for k=1:2
    window=10:599;
    [~,loc]=findpeaks(F(window,k)); idx=window(loc);
    fit=lsqcurvefit(@(p,x)p(1)*exp(-x/p(2))+p(3),[1 30 0],r_A(idx),...
        F(idx,k),[0 1e-8 -1],[Inf Inf Inf],optimoptions('lsqcurvefit','Display','off','FunctionTolerance',1e-13));
    xi_A(k)=fit(2);
end
kernel_diameter_nm=[1.20 (2*sqrt(8)+1)*0.4];
out=fullfile(here,'output'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,'gaussian_comparison.mat'),'kernel_diameter_nm','r_A','F','xi_A','POP_wide');
fprintf('Correlation lengths (Angstrom): %.5f, %.5f.\n',xi_A);
end
