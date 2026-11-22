% Cross-check the whole model against the three overlapping simulated regions.
here=fileparts(mfilename('fullpath')); rootDir=fileparts(fileparts(here));
fieldDir=fullfile(rootDir,'08_polar&domain&POP'); addpath(fieldDir);
m=load(fullfile(rootDir,'06_final_coordinates','Model_output_02_Final_model.mat'),'A_trace','B_tr_rep','Blat_final','polar_set_f');
labels=load(fullfile(fieldDir,'input','reviewed_domains.mat'),'domain_id');
xyz=pmn.lattice(m.Blat_final); use_BLI=[xyz,double(labels.domain_id)];
res=.3079; A_pos=(double(m.A_trace)+2)/3*res; B_pos=(double(m.B_tr_rep)+2)/3*res;
Model=[A_pos;B_pos]'; all_atom_types=zeros(1,size(Model,2));
th1=1.077138100686118e+02;
th2h=1.211676920356952e+02;
th2l=67.352164167361650;
th3=80.806046134445040;

%%
old_polar_gau_3=zeros(size(B_pos,1),3,3);
new_polar_gau_3=zeros(size(B_pos,1),3,3);
D_3=zeros(size(Model,2),3);

for i=1:3

if i==1
    load(fullfile(here,'part_a','part3a_pos&polar.mat'))
    part_Bind=B_pos(:,2)<th3;
    part_ind=Model(2,:)<th3;
    id_name='part3a';
elseif i==2
 load(fullfile(here,'part_b','part3b_pos&polar.mat'))
part_Bind=B_pos(:,2)<th2h & B_pos(:,2)>th2l;
    part_ind=Model(2,:)<th2h & Model(2,:)>th2l;
 id_name='part3b';
elseif i==3
 load(fullfile(here,'part_c','part3c_pos&polar.mat'))
part_Bind=B_pos(:,2)>th1;
part_ind=Model(2,:)>th1;
 id_name='part3c';
end
[idx,D]=knnsearch(maxXYZ_shift2,new_pos');
res=0.3079;
D=D*res/3;
previous=all_atom_types(part_ind);
assert(all(previous==0 | previous==reshape(atom_types,1,[])),'Overlap species mismatch.');
all_atom_types(part_ind)=reshape(atom_types,1,[]);
D_3(part_ind,i)=D;
old_polar_gau_3(part_Bind,:,i)=old_polar_gau;
new_polar_gau_3(part_Bind,:,i)=new_polar_gau;

figure
hold on
histogram(D(atom_types==82),DisplayName='Pb')
histogram(D(atom_types==41),DisplayName='Nb')
histogram(D(atom_types==22),DisplayName='Ti')
histogram(D(atom_types==12),DisplayName='Mg')
legend
 title(sprintf([id_name 'Pb rms%.3fA B rms%.3fA Nb rms%.3fA Mg rms%.3fA Ti rms%.3fA '],rms(D(atom_types==82)),rms(D(atom_types~=82)),rms(D(atom_types==41)),rms(D(atom_types==12)),rms(D(atom_types==22))))

diff_angle=acosd(sum(old_polar_gau.*new_polar_gau,2)./vecnorm(old_polar_gau,2,2)./vecnorm(new_polar_gau,2,2));
figure
histogram(diff_angle)
title(sprintf([id_name '模长的相关系数corr %.3f 平均夹角%.2f°'],corr(vecnorm(new_polar_gau,2,2),vecnorm(old_polar_gau,2,2)),mean(diff_angle)))



end


%% 求总的RMSD，一致性

old_polar_gau_all=zeros(size(B_pos,1),3);
new_polar_gau_all=zeros(size(B_pos,1),3);

D_3(D_3==0)=inf;
[V,I]=min(D_3,[],2);
ind=sub2ind(size(D_3),1:size(D_3,1),I');
D_all=D_3(ind);

assert(all(all_atom_types~=0));
D_3_B=D_3(all_atom_types~=82,:);
[V2,I2]=min(D_3_B,[],2);
indB=sub2ind(size(D_3_B),1:size(D_3_B,1),I2');


for i=1:3
    gaui=old_polar_gau_3(:,i,:);
old_polar_gau_all(:,i)=gaui(indB);

 gaui=new_polar_gau_3(:,i,:);
new_polar_gau_all(:,i)=gaui(indB);
end
old_polar_gau_all=old_polar_gau_all/3*0.3079;
new_polar_gau_all=new_polar_gau_all/3*0.3079;
%%
% Exact historical traced-atom subset for the manuscript's multislice panels.
% Its membership is external to the archived package. If supplied, use that
% mask for the subset panels. The full-data validation below is independent.
maskFile=fullfile(here,'Rind0208.mat');
if isfile(maskFile)
    selection=load(maskFile,'Rind');
    Rind=logical(reshape(selection.Rind,1,[]));
    assert(numel(Rind)==numel(all_atom_types));
sum(Rind)/numel(Rind)
% RAind=Rind(all_atom_types==82);
RBind=Rind(all_atom_types~=82);
figure
hold on
binedge=linspace(0,100,100);
histogram(D_all(all_atom_types==82 &Rind)*100,DisplayName='Pb',BinEdges=binedge)
histogram(D_all(all_atom_types~=82&Rind)*100,DisplayName='B site atoms',BinEdges=binedge)
% histogram(D_all(all_atom_types==22&Rind)*100,DisplayName='Ti',BinEdges=binedge)
% histogram(D_all(all_atom_types==12&Rind)*100,DisplayName='Mg',BinEdges=binedge)
% rms(D_all(RAind))
% rms(D_all(RBind))
legend
 title(sprintf(['all%.3fA Pb%.3fA B%.3fA Nb%.3fA Mg%.3fA Ti%.3fA '],rms(D_all(Rind)),rms(D_all(all_atom_types==82&Rind)),rms(D_all(all_atom_types~=82&Rind)),rms(D_all(all_atom_types==41&Rind)),rms(D_all(all_atom_types==12&Rind)),rms(D_all(all_atom_types==22&Rind))))
xlim([0 80])
daspect([80 2500 1])
xlabel('Deviation(pm)')
ylabel('Counts')
 box on
end

%%
diff_angle=acosd(sum(old_polar_gau_all.*new_polar_gau_all,2)./vecnorm(old_polar_gau_all,2,2)./vecnorm(new_polar_gau_all,2,2));
figure
histogram(diff_angle)
title(sprintf(['all 模长的相关系数corr %.3f 平均夹角%.2f°'],corr(vecnorm(old_polar_gau_all,2,2),vecnorm(new_polar_gau_all,2,2)),mean(diff_angle)))
%% 畴区展示
% curr_D_set=unique(eff_part_use_BLI(:,4));
for curr_D=23
    % curr_D=curr_D_set(i);
    if curr_D==0
        continue
    end

    ind=use_BLI(:,4)==curr_D;
    % eff_D_use_polar=eff_part_use_polar(ind,:);
    new_D_polar=new_polar_gau_all(ind,:);
    D_BLI=use_BLI(ind,:);
    old_D_polar=old_polar_gau_all(ind,:);
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
%% 余弦相似度
cs=sum(new_polar_gau_all.*old_polar_gau_all,"all")/sqrt(sum(vecnorm(new_polar_gau_all,2,2).^2))/sqrt(sum(vecnorm(old_polar_gau_all,2,2).^2))

if isfile(maskFile)
%% XYZ 分别比较
xyzchar=['x','z','y'];
for i=1:3
figure
points=[new_polar_gau_all(RBind,i) old_polar_gau_all(RBind,i)]';
[X,Y,H,XMesh,YMesh,ZMesh]=getKsdensity(points);

scatter(points(1,:),points(2,:),6,H,'filled')
title(sprintf('corr %.3f',corr(points(1,:)',points(2,:)')))
% yline(0)
% xline(0)
axis equal
box on
ax=gca;
if i==1
xlim([-1 1])
ylim([-1 1])
ax.XTick=linspace(-1,1,5);
ax.YTick=linspace(-1,1,5);
elseif i==2
xlim([-0.5 0.5])
ylim([-0.5 0.5])
ax.XTick=linspace(-0.5,0.5,5);
ax.YTick=linspace(-0.5,0.5,5);
elseif i==3
xlim([-1 1])
ylim([-1 1])
ax.XTick=linspace(-1,1,5);
ax.YTick=linspace(-1,1,5);
end



xlabel(['simulation polar component ' xyzchar(i) '(Å)'])
ylabel(['experiment polar component ' xyzchar(i) '(Å)'])
end
end
%% 不同模长极化对应的角度差
diff_angle=acosd(sum(old_polar_gau_all.*new_polar_gau_all,2)./vecnorm(old_polar_gau_all,2,2)./vecnorm(new_polar_gau_all,2,2));
polar_len=vecnorm(old_polar_gau_all,2,2);
bins=50;
LenBinEdge=linspace(min(polar_len),max(polar_len),bins);
diff_angle_arr=zeros(bins-1,1);
num_arr=zeros(bins-1,1);
for i=1:bins-1
ind=polar_len>=LenBinEdge(i)&polar_len<LenBinEdge(i+1);
diff_angle_arr(i)=mean(diff_angle(ind));
num_arr(i)=sum(ind);
end
figure
hold on
yyaxis left
plot(LenBinEdge(2:end)*res/3*100,diff_angle_arr,'-o')
ylabel('夹角（°）')
yyaxis right
ylabel('不同bin中点数')
plot(LenBinEdge(2:end)*res/3*100,num_arr,'-o')
xlabel('pm')


%% 极化之差的模长统计
res=0.3079;
figure
hold on
histogram(vecnorm(old_polar_gau_all,2,2)/3*res,DisplayName='Old Polar')
diff_vec=new_polar_gau_all-old_polar_gau_all;
% figure
histogram(vecnorm(diff_vec,2,2)/3*res,DisplayName='diff')
xlabel('pm')
title(sprintf('polarMean %.2fA diffMean %.2fA',mean(vecnorm(old_polar_gau_all,2,2)/3*res),mean(vecnorm(diff_vec,2,2)/3*res)))
%%
figure
scatter(vecnorm(diff_vec,2,2)/3*res,vecnorm(old_polar_gau_all,2,2)/3*res,6,'filled')
axis equal tight

% The archived map used the L1 norm of the neighborhood-mean displacement.
W=pmn.neighbors(xyz,6,Inf);
dims=max(xyz,[],1); site=sub2ind(dims,xyz(:,1),xyz(:,2),xyz(:,3));
D_map02_new_polar=zeros(dims); D_map02_old_polar=zeros(dims);
D_map02_new_polar(site)=sum(abs(W*new_polar_gau_all),2);
D_map02_old_polar(site)=sum(abs(W*old_polar_gau_all),2);
cage=zeros(size(B_pos));
for j=2:9, cage=cage+double(m.A_trace(m.polar_set_f(:,j),:))/8; end
v=pmn.neighbors(xyz,3,25/18)*((double(m.B_tr_rep)-cage)/3);
D_map=zeros(dims);D_map(site)=sum(abs(W*v),2);
%%
nonzero=D_map02_old_polar~=0;
figure
points2=[D_map02_old_polar(nonzero) D_map02_new_polar(nonzero)]';
[X,Y,H,XMesh,YMesh,ZMesh]=getKsdensity(points2);
scatter(points2(1,:),points2(2,:),6,H,'filled')
axis tight equal
title(sprintf('corr %.2f',corr(D_map02_new_polar(nonzero),D_map02_old_polar(nonzero))))
xlabel(['simulation POP'])
ylabel(['experiment POP'])
xlim([0 14])
ylim([0 14])

%%
corr(D_map(:),D_map02_old_polar(:))
corr(D_map(:),D_map02_new_polar(:))
corr(D_map02_old_polar(:),D_map02_new_polar(:))
save(fullfile(here,'new_polar_gau&POPmap.mat'),'D_map02_new_polar',"new_polar_gau_all","old_polar_gau_all");
%%
function [X,Y,H,XMesh,YMesh,ZMesh]=getKsdensity(points)

X=points(1,:)';
Y=points(2,:)';
n=50;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);

end
