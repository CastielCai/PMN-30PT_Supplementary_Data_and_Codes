%% Main Atom Tracing
% code to perform polynomial tracing on reconstruction volume
% reconstruction volume should be upsampled by 3*3*3 linear interpolation
% each local maximum is fitted with 9*9*9 voxels (3*3*3 before interpolation)
% by 4th order polynomial equation to get the position

% This code trace each sub area reconstruction of Ta film, then patch them
% together, which follows the scheme of sAET.

%% Step 1: trace each sub area reconstruction
clear
% for nrecon=1%:65

here=fileparts(mfilename('fullpath'));
recon=fullfile(fileparts(fileparts(here)),'1_recon','part_b');
Recon_filename=matfile(fullfile(recon,'PMNPT_noisy_part3b_1226_b60_1000iter_defocus.mat'));

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
% FinalVol_single = single(FinalVol);
% FinalVol_u8=uint8(FinalVol);
% % 
% save([inputDir '\atom_tracing\FinalVol_single.mat'],"FinalVol_single")

%%
% get the local maxima from the reconstruction volume
FinalVol = double(FinalVol);
% FinalVol =gpuArray(FinalVol_u8);

se = strel3d(25);  %% voxel
% se= strel('sphere',7);
dilatedBW = imdilate(FinalVol,se);
maxPos = find(FinalVol==dilatedBW & FinalVol>Th); % 膨胀运算的特性
maxVals = FinalVol(maxPos);
[orderVals,sortInd] = sort(maxVals,'descend');
%%
effind=orderVals>0.0002;%这个截断需要根据直方图调整
maxPos = maxPos(sortInd);
maxPos=maxPos(effind);
maxNum = min(numpeak,length(maxPos));
%%

fprintf(1,'numpeak = %d \n',length(maxPos));
%
maxXYZ = zeros(length(maxPos),3);
for i=1:length(maxPos)
    [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(i));
    maxXYZ(i,:) = [xx yy zz];  
end% clear Dsetvol dilatedBW Xi Yi Zi
% save([inputDir '\atom_tracing\maxXYZ.mat'],"maxXYZ")
save('0904maxXYZ.mat','maxXYZ')
%%
atom_pos=maxXYZ';
rec=FinalVol;
% 去除边界的局部极值点
b1 = find(atom_pos(1,:)<3 | atom_pos(1,:)>size(rec,1)-3);
b2 = find(atom_pos(2,:)<3 | atom_pos(2,:)>size(rec,2)-3);
b3 = find(atom_pos(3,:)<3 | atom_pos(3,:)>size(rec,3)-3);
bT = union(union(b1,b2),b3);
atom_pos(:,bT) = [];
% tmp_atom_pos = atom_pos;


maxXYZ=atom_pos';
%% 如果有初始模型可以直接跳到这一步
res=0.3079;
maxXYZ=importdata(fullfile(here,'part3b_Model_1226.mat')).Model';
maxXYZ=double(maxXYZ);
maxXYZ=maxXYZ/res*3;
% maxXYZ=oriModel_t5*3+10-2-2;
% maxXYZ=round(maxXYZ);
maxXYZ(:,1)=maxXYZ(:,1)-min(maxXYZ(:,1))+81-24;
maxXYZ(:,2)=maxXYZ(:,2)-min(maxXYZ(:,2))+41-13;
maxXYZ(:,3)=maxXYZ(:,3)-min(maxXYZ(:,3))+63-32;
% maxXYZ(:,1)=maxXYZ(:,1)-3;
% maxXYZ(:,2)=maxXYZ(:,2)-10;
% maxXYZ(:,3)=maxXYZ(:,3)+3;
loopTracing(maxXYZ,FinalVol,2,3)
%% 整体平移对齐初始模型和重构
options=optimset('MaxFunEvals',1e50,'Display','iter');
% parasup=[1e-3 1e-3 1e-3 1e-3 1e-3 1e-3];
% parasup=[0 0 0 -0.131003268377540	0.164529854614848	0.0347833458867639];
% parasup=[   -0.2620   -0.1114    1.2378 ];
ffval=inf;
fpara=zeros(1,3);
for zdiva=-10:5:10
    for xdiva=-10:5:10
        for ydiva=-10:5:10
            parasup=[xdiva ydiva zdiva];
            % parasup=parasup+randn(size(parasup))*0.5;
            % parasup=parasup;
            % [parasupf,fvalsup,flagsup]=fminsearch('register_two_model',gather(gModel),options,gather(gnewModel),parasup)
            try
            [parasupf,fvalsup,flagsup]=fminsearch('match_orimodel_and_rec',parasup,options,FinalVol,(maxXYZ))
            if fvalsup< ffval
                ffval=fvalsup;
                fpara=parasupf;
            end
            catch
            end
        end
    end
end
fpara
ffval
%%
maxXYZ_shift=maxXYZ+fpara;
% maxXYZ_shift=round(maxXYZ_shift);
loopTracing(round(maxXYZ_shift),FinalVol,2,3)
%% 以初始模型的点为基础在其附近寻找极大值
NewmaxXYZ=round(maxXYZ_shift);
% mysearchradi=5;
cropHalfSize=5;
maxIterNum=5;

for i=1:size(maxXYZ,1)
% for i=1

%     disp(i)
    cnt=0;
    flag=true;

   
    tpos=maxXYZ(i,:);
    while cnt<maxIterNum & flag
         cropXind = NewmaxXYZ(i,1) + (-cropHalfSize:cropHalfSize);
    cropYind = NewmaxXYZ(i,2) + (-cropHalfSize:cropHalfSize);
    cropZind = NewmaxXYZ(i,3) + (-cropHalfSize:cropHalfSize);
    % try
    cropVol = FinalVol(cropXind,cropYind,cropZind);
        [M,I]=max(cropVol(:));
       
        [i1,i2,i3]=ind2sub(size(cropVol),I);
        if sum(([i1,i2,i3]-cropHalfSize-1)==0)==3
            flag=false;
        else
            NewmaxXYZ(i,:)=NewmaxXYZ(i,:)+[i1,i2,i3]-cropHalfSize-1;
%             NewmaxXYZ(i,:)
        end
    cnt=cnt+1;
    end
% loopTracing([i1,i2,i3]',cropVol,2,1)
end
loopTracing(NewmaxXYZ,FinalVol,2,3)
maxXYZ=NewmaxXYZ;
%% gauss tracing
DataMatrix=FinalVol;
start_pos=maxXYZ_shift';
BoxSize1=9;
% FinalVol = double(importdata([inputDir 'atom_tracing\FinalVol_single.mat']));
[box1coordinates, BoxCenter1, BoxRadius1, ~] = create_box(BoxSize1);
% [box2coordinates, BoxCenter2, BoxRadius2, sphere2] = create_box(BoxSize2);
DM = DataMatrix;
numpos = size(start_pos,2);
new_pos=zeros(3,numpos);
% pos_notmoved=zeros(1,numpos);
stats_fit          = cell(1,numpos);
stats_peak_bg      = zeros(1,numpos);
stats_peak_hight   = zeros(1,numpos);
stats_peak_FWHM    = zeros(3,numpos);
stats_localmax_pos = zeros(numpos,3);
stats_shiftdist01  = zeros(1,numpos);
stats_shiftdist02  = zeros(1,numpos);
stats_resnorm      = zeros(1,numpos);
stats_residual     = cell(1,numpos);
for ii = 1:numpos
    disp(ii)
    curr_pos = start_pos(:,ii);
    %Get Data ROI
    x=round(curr_pos(1)); y=round(curr_pos(2)); z=round(curr_pos(3));
    DataBox1 = get_data_box(DM,x,y,z,BoxRadius1);
    %Get local max 之前得到的已经是局部最大了，这里不再需要
    [~,mind] = max(DataBox1(:));
    [mx,my,mz] = ind2sub(size(DataBox1),mind);
    nx = x+mx-BoxCenter1; ny = y+my-BoxCenter1; nz = z+mz-BoxCenter1;
    %Start values and bounds for fitting       
    %(x1:Bg; x2:Hight; x3,4,5：position;  x6: x7: x8:    x9:psi; x10:theta; x11:phi)
    ymax_init       = DataBox1(BoxCenter1,BoxCenter1,BoxCenter1);
    fit_param_init  = [0,  ymax_init, 0,  0,  0,  0.2,  0.2,  0.8,  0.0,  0.0,  0.0];
    fixed           = [0,          0, 0,  0,  0,    0,    0,    0,    0,    0,    0];
    lb              = [0,    0   ,  -1*BoxRadius1, -1*BoxRadius1, -1*BoxRadius1,    0,    0,    0,  -pi,   0,  -pi];
    ub              = [inf,  inf ,     BoxRadius1,    BoxRadius1,    BoxRadius1,  inf,  inf,  inf,   pi,  pi,   pi];
    
    %perform the fit
    [fit_resultI,resnorm,residual] = fit_gauss3D_PD(fit_param_init, box1coordinates, DataBox1, fixed, lb, ub);
 
    %calc position of new atom and FWHM
    new_x = x+fit_resultI(3);
    new_y = y+fit_resultI(4);
    new_z = z+fit_resultI(5);
    dist01 = pdist2(curr_pos' ,[new_x,new_y,new_z]);
    dist02 = pdist2([nx,ny,nz],[new_x,new_y,new_z]);
    %{
    %integer new position
    rx=round(new_x);
    ry=round(new_y);
    rz=round(new_z);
    %sub pixel offset to integer position
    dx=new_x-rx;
    dy=new_y-ry;
    dz=new_z-rz;
    %substract fit from data
    fit_resultI(3:5)=[dx dy dz];                            %use subpixel offsets
    FitBox2 = calc_gauss3D_PD(fit_resultI, box2coordinates); %calculate gaussian
    FitBox2 = FitBox2-fit_resultI(1);                        %substract background

    DataBox2 = get_data_box(DM,rx,ry,rz,BoxRadius2); %get data
    DiffBox2 = DataBox2-FitBox2;                           %substact gaussian fit
    DiffBox2 = max(0,DiffBox2);                            %remove negative values

    %update data for next iteration
    DM(rx-BoxRadius2:rx+BoxRadius2,ry-BoxRadius2:ry+BoxRadius2,rz-BoxRadius2:rz+BoxRadius2)=DiffBox2;
    %}
    %save statistics
    stats_fit{ii}            = fit_resultI;
    stats_peak_bg(ii)        = fit_resultI(1);
    stats_peak_hight(ii)     = fit_resultI(2);
    stats_peak_FWHM(:,ii)    = 2.3548./(sqrt(2)*fit_resultI(6:8));
    stats_localmax_pos(ii,:) = [nx,ny,nz]';
    stats_shiftdist01(ii)    = dist01;
    stats_shiftdist02(ii)    = dist02;
    stats_resnorm(ii)        = resnorm;
    stats_residual{ii}       = residual(:);
    %save atom position
    c_newpos = [new_x;new_y;new_z];
    new_pos(:,ii) = c_newpos;
end
%--------------------------------------------------------------------------
stats.fit_result      = stats_fit;
stats.peak_bg         = stats_peak_bg;
stats.peak_hight      = stats_peak_hight;
stats.peak_FWHM       = stats_peak_FWHM;
stats.localmax_pos    = stats_localmax_pos;
stats.shiftdist_init  = stats_shiftdist01;     % distance with initial position
stats.shiftdist_local = stats_shiftdist02;     % distance with local max
stats.resnorm         = stats_resnorm;
stats.residual        = stats_residual;



%%
loopTracing(new_pos,FinalVol,2,3)
save('PMNPT_rec_defocus_partb2_1128_gauss_tracing_round1_shift.mat',"new_pos")
%% 对齐一次 再tracing
options=optimset('MaxFunEvals',1e50,'Display','iter');
% parasup=[1e-3 1e-3 1e-3 1e-3 1e-3 1e-3];
% parasup=[0 0 0 -0.131003268377540	0.164529854614848	0.0347833458867639];
resrms=100;
paraf=0;
for dx=[ 0 ]
    for dy=[ 0 ]
        for dz=[ 0 ]
parasup=[0 0 0 dx dy dz ];
% parasup=[0 0 0];
% parasup=parasup+randn(size(parasup))*0.5;
% parasup=parasup;
% [parasupf,fvalsup,flagsup]=fminsearch('register_two_model',gather(gModel),options,gather(gnewModel),parasup)
[parasupf,fvalsup,flagsup]=fminsearch('register_two_model_PMNPT_moveOri',parasup,options,maxXYZ_shift,(new_pos'),"rms")
if fvalsup<resrms
    resrms=fvalsup;
    paraf=parasupf;
end
        end
    end
end
%%
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
orimodel=maxXYZ_shift;
mpos=mean(orimodel);
maxXYZ_shift2=(orimodel-mpos)*rotal*rotbe*rotga+ori+mpos;
% maxXYZ_shift2=orimodel+paraf(4:6);
%%
% ind=sqrt((maxXYZ_shift2(:,1)-1022).^2+(maxXYZ_shift2(:,2)-1022).^2)<600;
% pd=pdist2(maxXYZ_shift2(ind,:),new_pos(:,ind)');
pd=pdist2(maxXYZ_shift2,new_pos');
 min5=mink(pd,1)*res/3;
 figure
 histogram(min5(1,:))
 rms(min5)
 %%
 atom_types=importdata(fullfile(here,'part3b_Model_1226.mat')).atom_types;
 save(fullfile(here,'PMNPT_noisy_part3b_1226_b60_1000iter_defocus_gausstracing.mat'),'maxXYZ_shift2','new_pos','atom_types')


%%

