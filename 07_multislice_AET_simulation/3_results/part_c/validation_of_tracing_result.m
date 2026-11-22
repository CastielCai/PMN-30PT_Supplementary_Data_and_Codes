%% Comparison of experimental and multislice positions and displacements
here=fileparts(mfilename('fullpath'));
simulation=fileparts(fileparts(here)); rootDir=fileparts(simulation);
load(fullfile(simulation,'2_tracing','part_c','PMNPT_part3c_1224_rec_500iter_bm3d_defocus_gausstracing0105.mat'),'new_pos','maxXYZ_shift2','atom_types');
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
% Optional interactive intensity and layer inspection of the reconstructed volume.
reviewReconstruction=false;
if reviewReconstruction
    reconfile=fullfile(simulation,'1_recon','part_c','PMNPT_part3c_1231_rec_gauss_bm3d_500iter_0102_defocus.mat');
    Recon_filename=matfile(reconfile);
    addpath(fullfile(simulation,'2_tracing'));
    slice=B_pos(:,2)>th1;
    B_layer_order=double(m.Blat_final(slice,2))-min(double(m.Blat_final(:,2)))+1;
% Th: intensity threshold for the local maxima pixel
% local maxima with intensity less than this value will not be traced
% because they are way too weak to become actual atoms
MaxIter = 14;   CritIter = 7;   Th = 0;

% numpeak: maxmium number of peaks to trace
numpeak=200000;
SearchRad = 3;  saveInterval = 5000; 

% if exist([inputDir '11_Model\atom_tracing'],'dir')==0
%     mkdir([inputDir '11_Model\atom_tracing'])
% end

%
% load reconstructed 3D volume
OBJ=Recon_filename.OBJ;
Recon=OBJ.reconstruction;
Dsetvol.reconstruction = Recon;
Dsetvol.Dim1 = size(Dsetvol.reconstruction,1);
Dsetvol.Dim2 = size(Dsetvol.reconstruction,2);
Dsetvol.Dim3 = size(Dsetvol.reconstruction,3);
%Dsetvol在原来的Recon尺寸上pad了更多的0
dimx=Dsetvol.Dim1;
dimy=Dsetvol.Dim2;
dimz=Dsetvol.Dim3;
Dsetvol=croppedOut(Dsetvol.reconstruction,[dimx,dimy,dimz]);
%无效的crop，crop的尺寸等于被crop的对象
% Support = importdata([inputDir 'support0408.mat']);
% Support = My_paddzero(Support,[dimx,dimy,dimz]);
% Dsetvol = Dsetvol(:,51:549,210:389);

% clear Support
%引入support，使用support，删除support

%
% BoxSize0=3; %box size used for average when sorting peaks
% BoxSize1=9; %box size used to find maxima
% BoxSize2=9; %box size used box for fitting off-center gauss
% BoxSize3=7; %used to compute visualization matrix
%
% upsampling the reconstruction matrix by 3*3*3 by linear interpolation
% better to run at super conputer since the size of the interpolated
% volume will be larger than 16G
xx = (1:size(Dsetvol,1)) - round((size(Dsetvol,1)+1)/2);
yy = (1:size(Dsetvol,2)) - round((size(Dsetvol,2)+1)/2);
zz = (1:size(Dsetvol,3)) - round((size(Dsetvol,3)+1)/2);

xxi = ((3*xx(1)):(xx(end)*3))/3;
yyi = ((3*yy(1)):(yy(end)*3))/3;
zzi = ((3*zz(1)):(zz(end)*3))/3;

xxi = xxi(3:end);
yyi = yyi(3:end);
zzi = zzi(3:end);

[Y,X,Z] = meshgrid(yy,xx,zz);
[Yi,Xi,Zi] = meshgrid(yyi,xxi,zzi);

Dsetvol = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'linear',0);
clear Xi Yi Zi
FinalVol = My_paddzero(Dsetvol,size(Dsetvol)+20);
%%
 loopTracing2(new_pos,maxXYZ_shift2',FinalVol,2,3) 
 %% RMSD
 options=optimset('MaxFunEvals',1e50,'Display','iter');

parasup=[0 0 0 0 0 0];
[paraf,fvalsup,flagsup]=fminsearch('register_two_model_PMNPT_moveOri',parasup,options,maxXYZ_shift2,(new_pos'),"rms")
%
 al=paraf(1)/180*pi;
be=paraf(2)/180*pi;
ga=paraf(3)/180*pi;
ori=paraf(4:6);
rotal=[cos(al) -sin(al) 0;
    sin(al) cos(al) 0
    0 0 1]';
rotbe=[cos(be) 0 -sin(be);
    0 1 0;
    sin(be) 0 cos(be)]';
rotga=[1 0 0;
    0 cos(ga) -sin(ga);
    0 sin(ga) cos(ga);]';
orimodel=maxXYZ_shift2;
mpos=mean(orimodel);
old_pos=(orimodel-mpos)*rotal*rotbe*rotga+ori+mpos;

[idx,D]=knnsearch(old_pos,new_pos');
res=0.3079;
D=D*res/3;
rmsd=rms(D);


%%  get local intensity of tracing result 



%%
tic
for radius=[1 5 10 15]   

[Xq,Yq,Zq]=get_sphere_sub(radius);
pos_round=round(new_pos);
X_pos_round=Xq+pos_round(1,:);
Y_pos_round=Yq+pos_round(2,:);
Z_pos_round=Zq+pos_round(3,:);
ind=sub2ind(size(FinalVol),X_pos_round(:),Y_pos_round(:),Z_pos_round(:));
ind=reshape(ind,size(Xq,1),[]);
subset=FinalVol(ind);
if radius==1
local_max_r1t=max(subset,[],1);
elseif radius==5
local_max_r5t=max(subset,[],1);  
elseif radius==10
local_max_r10t=max(subset,[],1);  
elseif radius==15
local_max_r15t=max(subset,[],1);  
end
end
toc
%%

radius=15;
col=-radius:radius;
% [xq,yq,zq]=ndgrid(-radius:radius,-radius:radius,-radius:radius);
subrec=FinalVol(col+round(pos(1)),col+round(pos(2)),col+round(pos(3)));
loop(subrec,2)

%%
% tic
% 
% local_max_r1=zeros(1,size(new_pos,2));
% for i=1:size(new_pos,2)
%     ind=get_sphere_ind(1,size(FinalVol),new_pos(:,i));
%     local_max_r1(i)=max(FinalVol(ind));
% end
% local_max_r5=zeros(1,size(new_pos,2));
% for i=1:size(new_pos,2)
%     ind=get_sphere_ind(5,size(FinalVol),new_pos(:,i));
%     local_max_r5(i)=max(FinalVol(ind));
% end
% %
% local_max_r15=zeros(1,size(new_pos,2));
% for i=1:size(new_pos,2)
%     ind=get_sphere_ind(15,size(FinalVol),new_pos(:,i));
%     local_max_r15(i)=max(FinalVol(ind));
% end
% %
% local_max_r10=zeros(1,size(new_pos,2));
% for i=1:size(new_pos,2)
%     ind=get_sphere_ind(10,size(FinalVol),new_pos(:,i));
%     local_max_r10(i)=max(FinalVol(ind));
% end
% toc
%%
pos=[ 1012.53 31 867.215];
ind=get_sphere_ind(5,size(FinalVol),pos);
max(FinalVol(ind))
ind=get_sphere_ind(15,size(FinalVol),pos);
max(FinalVol(ind))
%%
figure
hold on
histogram(local_max_r1t(atom_types==82),DisplayName='Pb')
histogram(local_max_r1t(atom_types==41),DisplayName='Nb')
histogram(local_max_r1t(atom_types==22),DisplayName='Ti')
histogram(local_max_r1t(atom_types==12),DisplayName='Mg')
legend
%%
figure
hold on
histogram(local_max_r15t(atom_types==82),DisplayName='Pb')
histogram(local_max_r15t(atom_types==41),DisplayName='Nb')
histogram(local_max_r15t(atom_types==22),DisplayName='Ti')
histogram(local_max_r15t(atom_types==12),DisplayName='Mg')
legend
%%
comp= local_max_r10t==local_max_r5t;
show_ind=(~comp) &atom_types==12;
binedge=linspace(0,100,100);
figure
hold on
for atty=[82 41 22 12]
    rms(D(comp&atom_types==atty & D'<0.6))
    histogram(D(comp&atom_types==atty& D'<0.6 )*100,DisplayName=sprintf('%d',atty),BinEdges=binedge)
sum(comp(atom_types==atty)& D(atom_types==atty)'<0.6 )/sum(atom_types==atty)
end
%%
figure
hold on

histogram(D(atom_types==82),DisplayName='Pb')
histogram(D(atom_types==41),DisplayName='Nb')
histogram(D(atom_types==22),DisplayName='Ti')
histogram(D(atom_types==12),DisplayName='Mg')
legend
 title(sprintf('Pb rms%.3fA B rms%.3fA Nb rms%.3fA Mg rms%.3fA Ti rms%.3fA ',rms(D(atom_types==82)),rms(D(atom_types~=82)),rms(D(atom_types==41)),rms(D(atom_types==12)),rms(D(atom_types==22))))
%% 大偏离
figure
hold on
th=0.6;
for atty=[82 41 22 12]
    rms(D(atom_types'==atty&D<th ))
    histogram(D(atom_types'==atty&D<th ),DisplayName=sprintf('%d',atty))
% sum(comp(atom_types==atty))/sum(atom_types==atty)
sum(atom_types'==atty&D<th)/sum(atom_types'==atty)
end
%%
rms(D(D<0.5))
figure
histogram(D(D<0.5))

%%

figure
hold on
for atty=[82 41 22 12]
    rms(D(atom_types==atty))
    histogram(D(atom_types==atty ),DisplayName=sprintf('%d',atty))
% sum(comp(atom_types==atty))/sum(atom_types==atty)
end
%%
rms(D(comp&atom_types~=82))
rms(D(atom_types~=82))
%% show single atom position
show_ind=D>th;
pps=find(show_ind);
for i=15%:size(pps)
    local_max_r5t(pps(i))
  local_max_r10t(pps(i))

pos=new_pos(:,pps(i));
radius=5;
col=-radius:radius;
% [xq,yq,zq]=ndgrid(-radius:radius,-radius:radius,-radius:radius);
subrec1=FinalVol(col+round(pos(1)),col+round(pos(2)),col+round(pos(3)));
ind=get_sphere_ind(radius,size(subrec1),[1 1 1]+radius);
maxposind1=find(subrec1==max(subrec1(ind)));
[x,y,z]=ind2sub(size(subrec1),maxposind1);
maxpos1=[x y z];

radius=10;
col=-radius:radius;
% [xq,yq,zq]=ndgrid(-radius:radius,-radius:radius,-radius:radius);

subrec2=FinalVol(col+round(pos(1)),col+round(pos(2)),col+round(pos(3)));
ind=get_sphere_ind(radius,size(subrec2),[1 1 1]+radius);
maxposind2=find(subrec2==max(subrec2(ind)));
[x,y,z]=ind2sub(size(subrec2),maxposind2);
maxpos2=[x y z];

loopTracing(maxpos1,subrec1,2,1)
loopTracing(maxpos2,subrec2,2,1)
% loop2(subrec1,subrec2,2)
colormap parula
% pause
end



%% show spatial distribute
%% tracing 失败的
show_ind=(~comp) &atom_types==12;
loopTracing2(new_pos(:,show_ind),maxXYZ_shift2(show_ind,:)',FinalVol,2,3) 
%%
% loopTracing2(new_pos(:,atom_types==22),new_pos(:,show_ind),FinalVol,2,3) 

% loopTracing2(new_pos(:,comp & atom_types==12 &D'>0.5),maxXYZ_shift2',FinalVol,2,3) 

% loopTracing2(new_pos,maxXYZ_shift2',FinalVol,2,3) 
%% show deviation by layer

A_old_pos=maxXYZ_shift2(atom_types==82,:);
A_new_pos=new_pos(:,atom_types==82)';
B_old_pos=maxXYZ_shift2(atom_types~=82,:);
B_new_pos=new_pos(:,atom_types~=82)';
B_atom_types=atom_types(atom_types~=82);
D_B=D(atom_types~=82);
%%
% ind2=comp(atom_types~=82) &B_atom_types==12
% ind2=comp(atom_types~=82) &B_atom_types==12 & D_B'>0.5 ;
ind2=comp(atom_types~=82) &B_atom_types==12 &local_max_r15t(atom_types~=82)<30;
for j=1:2
figure
tiledlayout(2,4)
for i=1:8
if i+(j-1)*8>13
    break
end
nexttile
hold on
%    

layer_ind=B_layer_order==i+(j-1)*8;
% old_pos_p=B_old_pos(ind,:);
% new_pos_p=B_new_pos(ind,:);

% [idx,D]=knnsearch(new_pos_p,old_pos_p);
% D=D*res/3;

% scatter(new_pos_p(idx,1),new_pos_p(idx,3),10,D,'filled');

% scatter(new_pos_p(idx,1),new_pos_p(idx,3),10,D,'filled');
% scatter(old_pos_p(:,1),old_pos_p(:,3),10,D,'filled');


  scatter(B_new_pos(layer_ind &B_atom_types==41,1),B_new_pos(layer_ind &B_atom_types==41,3),10,'filled',DisplayName='Nb');
    scatter(B_new_pos(layer_ind &B_atom_types==22,1),B_new_pos(layer_ind &B_atom_types==22,3),10,'filled',DisplayName='Ti');
  scatter(B_new_pos(layer_ind &B_atom_types==12,1),B_new_pos(layer_ind &B_atom_types==12,3),10,'filled',DisplayName='Mg');
 % scatter(B_new_pos(layer_ind &ind2,1),B_new_pos(layer_ind&ind2,3),30,'k','filled','^',DisplayName='not traced Nb');
  scatter(B_new_pos(layer_ind &ind2,1),B_new_pos(layer_ind&ind2,3),30,'k','filled','^',DisplayName='强度小 Mg');
% scatter(old_pos_p(idx,1),old_pos_p(idx,3),10,D,'filled');
axis equal tight
legend
% view(0,90)
colorbar
% clim([0.5  0.9])
% title(sprintf('B rms%.3fA Nb rms%.3fA Mg rms%.3fA Ti rms%.3fA  layer %d',rms(D),rms(D(B_atom_types(layer_ind)==41)),rms(D(B_atom_types(layer_ind)==12)),rms(D(B_atom_types(layer_ind)==22)),i+(j-1)*8))
 title(sprintf(' layer %d',i+(j-1)*8))

% nexttile
% hold on
% histogram(D(B_atom_types(ind)==41),DisplayName='Nb')
% histogram(D(B_atom_types(ind)==22),DisplayName='Ti')
% histogram(D(B_atom_types(ind)==12),DisplayName='Mg')
% legend


end
end
end
part_Dind=B_pos(:,2)>th1;
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
save(fullfile(here,'part3c_pos&polar.mat'),'A_new_pos','B_new_pos',"A_old_pos","B_old_pos","new_polar_gau","old_polar_gau",'maxXYZ_shift2','new_pos','atom_types');
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


