%% Comparison of experimental and multislice positions and displacements
here=fileparts(mfilename('fullpath'));
simulation=fileparts(fileparts(here)); rootDir=fileparts(simulation);
load(fullfile(simulation,'2_tracing','part_b','PMNPT_noisy_part3b_1226_b60_1000iter_defocus_gausstracing.mat'),'new_pos','maxXYZ_shift2','atom_types');
[~,D]=knnsearch(maxXYZ_shift2,new_pos');
res=.3079; D=D*res/3; rmsd=rms(D);
A_old_pos=maxXYZ_shift2(atom_types==82,:); A_new_pos=new_pos(:,atom_types==82)';
B_old_pos=maxXYZ_shift2(atom_types~=82,:); B_new_pos=new_pos(:,atom_types~=82)';
% Reference coordinates are deposited once in 06. Labels are only used for optional displays.
m=load(fullfile(rootDir,'06_final_coordinates','Model_output_02_Final_model.mat'),'B_tr_rep','Blat_final');
a=load(fullfile(rootDir,'08_polar&domain&POP','input','reviewed_domains.mat'),'domain_id');
B_pos=(double(m.B_tr_rep)+2)/3*res;
use_BLI=[double(m.Blat_final),double(a.domain_id)];
th1=107.7138100686118; th2l=67.35216416736165; th2h=121.1676920356952; th3=80.80604613444504;
part_Dind=B_pos(:,2)>th2l&B_pos(:,2)<th2h;
part_use_BLI=use_BLI(part_Dind,:);
% Each region has its own grid dimensions and number of B sites.
grid_polar_gau_old=[]; grid_polar_gau_new=[];
old_polar_gau=[]; new_polar_gau=[];
%%
tic
[I,B]=knnsearch(A_old_pos,B_old_pos,'K',8);
effind1=max(B,[],2)<50;
eff_I1=I(effind1,:);
old_A_cell=reshape(A_old_pos(eff_I1,:),size(eff_I1,1),size(eff_I1,2),[]);
old_A_cell_mean_pos=mean(old_A_cell,2);
old_polar=B_old_pos(effind1,:)-squeeze(old_A_cell_mean_pos);

% The experimental nearest-Pb indices are used for both coordinate sets.
% effind2=max(B,[],2)<50;
% eff_I2=I(effind2,:);
new_A_cell=reshape(A_new_pos(eff_I1,:),size(eff_I1,1),size(eff_I1,2),[]);
new_A_cell_mean_pos=mean(new_A_cell,2);
new_polar=B_new_pos(effind1,:)-squeeze(new_A_cell_mean_pos);
toc
%%
curr_polar=old_polar;

part_use_BLI(:,1:3)=part_use_BLI(:,1:3)-min(part_use_BLI(:,1:3),[],1)+[1 1 1 ];
eff_part_use_BLI=part_use_BLI(effind1,:);

[xq,yq,zq]=ndgrid(min(part_use_BLI(:,1)):max(part_use_BLI(:,1)),min(part_use_BLI(:,2)):max(part_use_BLI(:,2)),min(part_use_BLI(:,3)):max(part_use_BLI(:,3)));
Fx=scatteredInterpolant(eff_part_use_BLI(:,1),eff_part_use_BLI(:,2),eff_part_use_BLI(:,3),curr_polar(:,1));
Fy=scatteredInterpolant(eff_part_use_BLI(:,1),eff_part_use_BLI(:,2),eff_part_use_BLI(:,3),curr_polar(:,2));
Fz=scatteredInterpolant(eff_part_use_BLI(:,1),eff_part_use_BLI(:,2),eff_part_use_BLI(:,3),curr_polar(:,3));
grid_polar_x_old=Fx(xq(:),yq(:),zq(:));
grid_polar_y_old=Fy(xq(:),yq(:),zq(:));
grid_polar_z_old=Fz(xq(:),yq(:),zq(:));
grid_polar=reshape(cat(1,grid_polar_x_old,grid_polar_y_old,grid_polar_z_old),size(xq,1),size(xq,2),size(xq,3),[]);
gaukernel=fspecial3('gaussian',3,2);
for i=1:3
grid_polar_gau_old(:,:,:,i)=convn(grid_polar(:,:,:,i),gaukernel,'same');
end
%
curr_polar=new_polar;
eff_part_use_BLI=part_use_BLI(effind1,:);
% eff_part_use_BLI(:,1:3)=eff_part_use_BLI(:,1:3)-min(eff_part_use_BLI(:,1:3),[],1)+[1 1 1 ];
[xq,yq,zq]=ndgrid(min(part_use_BLI(:,1)):max(part_use_BLI(:,1)),min(part_use_BLI(:,2)):max(part_use_BLI(:,2)),min(part_use_BLI(:,3)):max(part_use_BLI(:,3)));
Fx=scatteredInterpolant(eff_part_use_BLI(:,1),eff_part_use_BLI(:,2),eff_part_use_BLI(:,3),curr_polar(:,1));
Fy=scatteredInterpolant(eff_part_use_BLI(:,1),eff_part_use_BLI(:,2),eff_part_use_BLI(:,3),curr_polar(:,2));
Fz=scatteredInterpolant(eff_part_use_BLI(:,1),eff_part_use_BLI(:,2),eff_part_use_BLI(:,3),curr_polar(:,3));
grid_polar_x_new=Fx(xq(:),yq(:),zq(:));
grid_polar_y_new=Fy(xq(:),yq(:),zq(:));
grid_polar_z_new=Fz(xq(:),yq(:),zq(:));
grid_polar=reshape(cat(1,grid_polar_x_new,grid_polar_y_new,grid_polar_z_new),size(xq,1),size(xq,2),size(xq,3),[]);
gaukernel=fspecial3('gaussian',3,2);
for i=1:3
grid_polar_gau_new(:,:,:,i)=convn(grid_polar(:,:,:,i),gaukernel,'same');
end
%% 比较前后的夹角区别和模长的相关程度
grid_polar_gau_old_r=reshape(grid_polar_gau_old,[],3);
grid_polar_gau_new_r=reshape(grid_polar_gau_new,[],3);
diff_angle=acosd(sum(grid_polar_gau_old_r.*grid_polar_gau_new_r,2)./vecnorm(grid_polar_gau_new_r,2,2)./vecnorm(grid_polar_gau_old_r,2,2));
figure
histogram(diff_angle)
title(sprintf('模长的相关系数corr %.3f 夹角的峰位16度 平均夹角%.2f°',corr(vecnorm(grid_polar_gau_old_r,2,2),vecnorm(grid_polar_gau_new_r,2,2)),mean(diff_angle)))
%% 看分畴的结果

ind=sub2ind(size(xq),part_use_BLI(:,1),part_use_BLI(:,2),part_use_BLI(:,3));
for i=1:3
gaux=grid_polar_gau_old(:,:,:,i);
old_polar_gau(:,i)=gaux(ind);
end

ind=sub2ind(size(xq),part_use_BLI(:,1),part_use_BLI(:,2),part_use_BLI(:,3));
for i=1:3
gaux=grid_polar_gau_new(:,:,:,i);
new_polar_gau(:,i)=gaux(ind);
end
save(fullfile(here,'part3b_pos&polar.mat'),'A_new_pos','B_new_pos',"A_old_pos","B_old_pos","new_polar_gau","old_polar_gau",'maxXYZ_shift2','new_pos','atom_types');
%%
curr_D_set=unique(eff_part_use_BLI(:,4));
for i=1:5:size(curr_D_set,1)
    curr_D=curr_D_set(i);
    if curr_D==0
        continue
    end

    ind=eff_part_use_BLI(:,4)==curr_D;
    new_D_polar=new_polar_gau(ind,:);
    D_BLI=eff_part_use_BLI(ind,:);
    old_D_polar=old_polar_gau(ind,:);
    figure
   % subplot(1,3,1)
   %  quiver3(D_BLI(:,1),D_BLI(:,2),D_BLI(:,3),eff_D_use_polar(:,1),eff_D_use_polar(:,2),eff_D_use_polar(:,3)); 
   subplot(1,2,1)
       quiver3(D_BLI(:,1),D_BLI(:,2),D_BLI(:,3),old_D_polar(:,1),old_D_polar(:,2),old_D_polar(:,3));
    title('experi')
   subplot(1,2,2)
quiver3(D_BLI(:,1),D_BLI(:,2),D_BLI(:,3),new_D_polar(:,1),new_D_polar(:,2),new_D_polar(:,3));
    title('multislice')

   sgtitle(sprintf('domian no. %.0f',curr_D))
end

%%
diff_angle=acosd(sum(new_polar.*old_polar,2)./vecnorm(new_polar,2,2)./vecnorm(old_polar,2,2));
figure
histogram(diff_angle)

%%
corr(vecnorm(new_polar,2,2),vecnorm(old_polar,2,2))
%%
figure
quiver3(B_old_pos(effind1,1),B_old_pos(effind1,2),B_old_pos(effind1,3),old_polar(:,1),old_polar(:,2),old_polar(:,3))
%%
function ind=get_sphere_ind(radius,RecSize,pos)
[xq,yq,zq]=ndgrid(-radius:radius,-radius:radius,-radius:radius);
pos=round(pos);

insphere=sqrt(xq(:).^2+yq(:).^2+zq(:).^2)<=radius; 

Xq=xq(insphere)+pos(1);
Yq=yq(insphere)+pos(2);
Zq=zq(insphere)+pos(3);

ind=sub2ind(RecSize,Xq,Yq,Zq);

end

function [Xq,Yq,Zq]=get_sphere_sub(radius)
[xq,yq,zq]=ndgrid(-radius:radius,-radius:radius,-radius:radius);


insphere=sqrt(xq(:).^2+yq(:).^2+zq(:).^2)<=radius; 

Xq=xq(insphere);
Yq=yq(insphere);
Zq=zq(insphere);

end


