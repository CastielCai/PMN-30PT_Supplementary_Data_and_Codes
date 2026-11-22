clc 
%% LOAD
here = fileparts(mfilename('fullpath'));
rootDir = fileparts(here);
addpath(fullfile(here,'src'),genpath(fullfile(rootDir,'02_RESIRE_package','src')));
Recvol0 = importdata(fullfile(rootDir,'04_final_reconstructoin_volume','PMN30PT_Recon_Volume.mat'));
Recvol0 = My_stripzero(Recvol0,[881 520 881]);
support = importdata(fullfile(here,'input','PMN30PT_support.mat'));
support(:,511:520,:)=0;

DM = Recvol0.*support;  % DataMatrix 
size0=size(DM);
 
%% find local maxima only for A
tic
SearchRad    = 8;
Th           = 200;
se           = strel3d(SearchRad);   %strel3d
dilatedBW    = imdilate(DM, se);
maxPos       = find( DM == dilatedBW & DM > Th );
maxXYZ = zeros(length(maxPos),3);
XZcen = round((size0([1,3]))/2);
clear pdist_arr
for ii=1:length(maxPos)
    [xx,yy,zz] = ind2sub(size(DM),maxPos(ii));
    maxXYZ(ii,:) = [xx yy zz];
    pdist_arr(ii) = pdist2([xx zz],XZcen);
end
toc
% loopTracing(maxXYZ,DM,2,2)
% loopTracing(maxXYZ,DM,3,2)
%--- filter out pos
Dth = max(pdist2([62 385; 360 75; 395 739; 742 378],XZcen));

%%
% maxPos_f = maxPos(pdist_arr'<Dth & ( maxXYZ(:,2)>75 | maxXYZ(:,2)<(size(maxXYZ,2)-SearchRad-1) ));
% maxPos_f = maxPos(pdist_arr'<Dth);% & ( maxXYZ(:,2)>75 | maxXYZ(:,2)<(size(maxXYZ,2)-SearchRad-1) ));
maxPos_f = maxPos;
maxXYZ_f = zeros(length(maxPos_f),3);
for ii=1:length(maxPos_f)
    [xx,yy,zz] = ind2sub(size(DM),maxPos_f(ii));
    maxXYZ_f(ii,:) = [xx yy zz];
end 
[maxXYZ_f,Ind]= clear_outer_pos3D(maxXYZ_f,support);
maxPos_f = maxPos_f(Ind);


loopTracing(maxXYZ_f,DM,2,2)
whos maxXYZ maxXYZ_f maxPos_f
% loopTracing(maxXYZ_f,DM,2,2)
% loopTracing(maxXYZ_f,DM,3,2)


%%
maxVals_f    = DM(maxPos_f);
[~,sortInd]  = sort(maxVals_f,'descend');
maxNum       = min( 100000, length(sortInd) );
maxXYZ_fs    = maxXYZ_f(sortInd(1:maxNum),:);
fprintf(1,'numpeak = %d \n',length(maxXYZ_fs));
clear dilatedBW
% loopTracing(maxXYZ,DM,2,2)
loopTracing(maxXYZ_fs,DM,2,2)
% save('local_max_pos','maxXYZ_fs', 'N')
min(maxXYZ_fs,[],1)
max(maxXYZ_fs,[],1)

%% manual delete out side A pos atom
Devide.X = 50:750;
Devide.Y = 90:530;
Devide.Z = 50:750;
fileN='No1A1';
%% Final version of Add & Delete
Pix_coordA = maxXYZ_fs;
if exist(['modify_' fileN '.mat'],'file')
    load(['modify_' fileN '.mat']);
    tempInfo.addpos = addpos;    tempInfo.deleteInd = deleteInd
else
    tempInfo.addpos =[];   tempInfo.deleteInd =[];
end

Model_manual_add_delete_v3(Pix_coordA', DM, Devide, fileN, tempInfo)
% this function creat an mat file of pos you modify manually

%%
load(['modify_' fileN '.mat'])
whos addpos deleteInd
modify_A = Pix_coordA';
modify_A(:,deleteInd) = [];
% modify_A = cat(2,modify_A,addpos);
%
max_final = modify_A';

 

%% extract lattice from local max
clc
coords = max_final;
start_round = round(mean(coords,1));
findN = 52;
[start_pos,pos_neighbor] = find_neighbour_pos_N_inset(start_round',coords',findN);
figure(1)
clf
scat3(pos_neighbor,15,'b');
hold on
scatter3(start_pos(1),start_pos(2),start_pos(3),25,'filled','MarkerEdgeColor','r')
start_idx = find_pos_ind_dim1(start_pos',coords);


%%

volind = cell(3,1);
volind{1} = round(start_pos(1)) + (-40:40);
volind{2} = round(start_pos(2)) + (-30:30);
volind{3} = round(start_pos(3)) + (-40:40);

cutvol = DM(volind{1},volind{2},volind{3});
volpos = [pos_neighbor,start_pos]-[volind{1}(1) volind{2}(1) volind{3}(1)]'+[1 1 1]';
loopTracing(volpos,cutvol,2,0.5,23)
fprintf('Check Lattice pos\n')
pause
% test 404   % 322   % 393
% find(coords(:,1)==404 & coords(:,2)==322 & coords(:,3)==393)


%% input centerpos and a,b,c vecotr pos manually
cubecen =  [41 31 41];
pospair1 = [50 31 51; 31 30 32];
pospair2 = [40 17 41; 40 43 42];
pospair3 = [51 30 32; 31 30 52];
SearchVecs_pixel_half(1,:) = mean([pospair1(1,:)-cubecen;cubecen-pospair1(2,:)]);
SearchVecs_pixel_half(2,:) = mean([pospair2(1,:)-cubecen;cubecen-pospair2(2,:)]);
SearchVecs_pixel_half(3,:) = mean([pospair3(1,:)-cubecen;cubecen-pospair3(2,:)])

% SearchVecs_pixel_half = SearchVecs_pixel(:,[1 3 5])
%%
% SearchVecs_pixel_half = fit_vec;
%%%%% update coords !!!!!
% coords = coords_extra;
SearchTh = 2;
pair_arr26 = Simple_cubic_search_0109(coords, SearchVecs_pixel_half, SearchTh);
% pair_arr26 = pair_arr;
pasum = sum(pair_arr26,2);
isoInd = find(pasum==0);
whos isoInd pair_arr26 

figure(1928); clf
scat3(coords); hold on
scat3(coords(isoInd,:)); hold on




    %% one by one check , [not neccesary]
close all
    for ii=1:size(pair_arr26,1)
        cind = pair_arr26(ii,:);
        cind(cind==0)=[];
            figure(24320);clf

        scat3(coords(cind,:));hold on
        scat3(coords(ii,:));
        scat3(coords(ii,:)+currSearchPos,'*');
        view([1 6 1])
        pause
    end


%% assignment Barvais Lattice Index
[ALI_arr,stats] = assign_lat_cubic26(pair_arr26,1);

%% extract effective lattice
close all
whos ALI_arr
Zind = find(ALI_arr(:,1)==0 & ALI_arr(:,2)==0 & ALI_arr(:,3)==0);
whos Zind
Aind = IndDiff(1:size(coords,1),Zind);
figure;scat3(coords(Zind,:));hold on
loopTracing(coords(Aind,:),DM,2,2);hold on
ALI_arr(1,:)

%% 
coords_extra = coords;%(Aind,:);
Alat_extra   = ALI_arr;%(Aind,:); 
%--- check repeat Lat index
rept_pair = [];
for ii=1:size(Alat_extra,1)
    if any(ii==Zind)
        continue 
    end
    clat = Alat_extra(ii,:);
    repInd = find( Alat_extra(ii+1:end,1)==clat(1) & Alat_extra(ii+1:end,2)==clat(2) & Alat_extra(ii+1:end,3)==clat(3) );
    if ~isempty(repInd)
        rept_pair(end+1,:) = [ii, repInd+ii];
    end
end
        if ~isempty(rept_pair)
            whos rept_pair
            % cat(2,Alat_extra(rept_pair(:,1),:),Alat_extra(rept_pair(:,2),:))
            % cat(2,pair_arr26(rept_pair(1,1),:)',pair_arr26(rept_pair(1,2),:)')
            % cat(2,coords(rept_pair(:,1),:),coords(rept_pair(:,2),:))
            figure;clf
            subplot(121)
            scat3(Alat_extra(rept_pair(:,1),:));axis equal
            subplot(122)
            scat3(coords_extra(rept_pair(:,1),:));hold on
            scat3(coords_extra(rept_pair(:,2),:),'.');hold on
        end
        
        %%
        for ckidx = 1:size(rept_pair,1)
            idxx = rept_pair(ckidx,:);
            % Alat_extra(idxx(1),:)
            c_round = pair_arr26(idxx(1),:)';
            c_round(c_round==0)=[];
            % Alat_extra(c_round,:)
            figure(3258);clf;
            subplot(121)
            scat3(Alat_extra(c_round,:));hold on
            scat3(Alat_extra(idxx(1),:),'r*');hold on
            c_round = pair_arr26(idxx(2),:)';
            c_round(c_round==0)=[];
            subplot(122)
            scat3(Alat_extra(c_round,:));hold on
            scat3(Alat_extra(idxx(2),:),'r*');hold on 
            pause
        end
 
%%  clear ALI_arr
ALI_arr(rept_pair(:,2),:) = repmat([0 0 0],[size(rept_pair,1),1]); 
Zind = find(ALI_arr(:,1)==0 & ALI_arr(:,2)==0 & ALI_arr(:,3)==0);
Zind = IndDiff(Zind,1); whos Zind

%% refine A lattice
    %% use LI_arr to fit lattice 
    Ydata = coords;       Xdata = ALI_arr;   
    % Ydata = fitpos_fill;       Xdata = Alat_fill;   
    Ydata(Zind,:) = [];   Xdata(Zind,:) = [];  
    svh = SearchVecs_pixel_half; 
    o1 = Ydata(1,1); o2 = Ydata(1,2); o3 = Ydata(1,3);
    a1 = svh(1,1); a2 = svh(1,2); a3 = svh(1,3); 
    b1 = svh(2,1); b2 = svh(2,2); b3 = svh(2,3); 
    c1 = svh(3,1); c2 = svh(3,2); c3 = svh(3,3);   
    %--------------------- 
    fun = @(x,xdata)x([1 2 3]) +x([4 5 6]).*xdata(:,1) +x([7 8 9]).*xdata(:,2) +x([10 11 12]).*xdata(:,3);
    x0 = [o1,o2,o3, a1,a2,a3, b1,b2,b3, c1,c2,c3];
    lb = [o1-4, o2-4, o3-4,  8, -2,  8, -2, -15, -2,    8, -2, -11];
    ub = [o1+4, o2+4, o3+4, 11,  2, 11,  2, -11,  2,   11,  2,  -8];
    
    x1 = lsqcurvefit(fun,x0,Xdata,Ydata,lb,ub);   % xdata: ALI_arr    ydata: c_A
    x = x1;
    % x2 = lsqcurvefit(fun,x0,Xdata,Ydata,lb,ub);   % xdata: ALI_arr    ydata: c_A
    % x = x2;
    fitA_pos = x([1 2 3])+x([4 5 6]).*Xdata(:,1)+x([7 8 9]).*Xdata(:,2)+x([10 11 12]).*Xdata(:,3);
    figure(2052);clf
    scat3(fitA_pos);hold on;pause
    scat3(Ydata,'.');hold on
%{
    fit_vec = [x([4 5 6]); 
             x([7 8 9]); 
             x([10 11 12])] ;
    save('fit_vec','fit_vec')
%}
    loopTracing2(fitA_pos',Ydata',DM,2,2) 

%% ----deviation
ds_arr = fitA_pos-Ydata;
for ii=1:size(Ydata,1)
    ds(ii) = norm(ds_arr(ii,:));
end
figure(1307);histogram(ds)

%% repair Lat & pos
repair_arr = zeros(size(fitA_pos,1),2);   % column 1:coords idx. col 2:distance
lat_r =  mean([norm(fit_vec(1,:)),norm(fit_vec(2,:)),norm(fit_vec(3,:))])/2;

for ii=1:size(fitA_pos,1)
    cpos = fitA_pos(ii,:);
    cds = pdist2(cpos,coords);
    [ds1,idx1] = min(cds);
    if ds1<lat_r 
        repair_arr(ii,:) = [idx1,ds1]; 
    end
end 
%--- cancel repeat
count = 0;
for ii=1:size(repair_arr,1)
    c_idx = repair_arr(ii,1);
    if c_idx==0
        continue 
    end 
    repInd = find(repair_arr(:,1)==c_idx);
    if length(repInd)>=2
        cdsarr = repair_arr(repInd,2);
        [~,minidx] = min(cdsarr);
        repInd(minidx)=[];
        repair_arr(repInd,:)=0;
        count = count+1;
    end  
end
Zidx = find(repair_arr(:,1)==0); whos Zidx
figure
scat3(fitA_pos(Zidx,:)); 

%%
% count = 0;
% for ii=1:size(Zidx,1)
%     cpos = fitA_pos(Zidx(ii),:);
%     cds = pdist2(cpos,coords);
%     [ds1,idx1] = min(cds);
%     if ds1<lat_r+2 && ~any(idx1==repair_arr(:,1))
%         repair_arr(Zidx(ii),:) = [idx1,ds1];  
%         count = count+1;
%     end
% end 

%% fill inner vacancy 
cIdx = find(repair_arr(:,1)~=0);
Lat_tofill    = Xdata(cIdx,:);
coords_tofill = coords(repair_arr(cIdx,1),:);
fit_tofill    = fitA_pos(cIdx,:);
whos Lat_tofill coords_tofill  fit_tofill 
% loopTracing(coords_tofill',DM,2,2)

[Alat_fill,coords_fill] = fill_Lattice_vacancy(Lat_tofill,coords_tofill);
[~,fit_fill] = fill_Lattice_vacancy(Lat_tofill,fit_tofill);
whos  Alat_fill coords_fill fit_fill
% loopTracing2(fit_fill',coords_fill',DM,2,2)
  


%%  Overall deviation
clc
Amax_final = coords_fill; 
ds_arr = fit_fill-Amax_final;
for ii=1:size(Amax_final,1)
    ds(ii) = norm(ds_arr(ii,:));
end
figure(1431);histogram(ds)  



%% Gen B pos 
vec_AB = sum(fit_vec,1)/2;
Alat_final = Alat_fill;
Blat_Init = Alat_final;

%% Polarset
    %% find polar set
polar_lat = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1];
pN = size(Blat_Init,1);
polar_set = zeros(pN,9);
polar_set(:,1) = 1:pN;
fullset = zeros(pN,1);
errorInd=[];
count=0;
for ii=1:pN
    curr_Blat = Blat_Init(ii,:);
    for jj=1:8
        c_lat = curr_Blat+polar_lat(jj,:);
        c_Aind = find(Alat_final(:,1)==c_lat(1) & Alat_final(:,2)==c_lat(2) & Alat_final(:,3)==c_lat(3));
        if ~isempty(c_Aind)
            if numel(c_Aind)>1
                count=count+1;
                errorInd(count,:) = c_Aind';
            end
            polar_set(ii,jj+1) = c_Aind(1);
        else
            polar_set(ii,jj+1) = 0;         
        end
    end
    fullset(ii) = ~any(polar_set(ii,:)==0);
end
% numel(unique(errorInd(:)))
goodUC_ind = find(fullset==1); whos goodUC_ind


%%
Blat_final = Blat_Init(goodUC_ind,:);
polar_set_f = polar_set(goodUC_ind,:);
whos Blat_final polar_set_f
Bpos_Init = zeros(length(goodUC_ind),3);
for bb=1:length(goodUC_ind)
    Bpos_Init(bb,:) = mean(Amax_final(polar_set_f(bb,2:9),:)); 
end 



%% 
save(fullfile(here,'output','Model_output_00_localMax.mat'), 'Alat_final', 'Amax_final', 'Blat_final', 'Bpos_Init', 'polar_set_f');




 
%% subfunction
function SearchVecs_pixel = mean_lattice_b(atom_pos,start_round,DataM)
    %% ----find neighbor atoms[start_round初始为空]
    findN = 52;
    [start_pos,pos_neighbor] = find_neighbour_pos_N_inset(start_round,atom_pos,findN);
    figure(1)
    clf
    scat3(pos_neighbor,15,'b')
    hold on
    scatter3(start_pos(1),start_pos(2),start_pos(3),25,'filled','MarkerEdgeColor','r')
    
    volind = cell(3,1);
    volind{1} = round(start_pos(1)) + (-40:40);
    volind{2} = round(start_pos(2)) + (-20:20);
    volind{3} = round(start_pos(3)) + (-40:40);
    
    cutvol = DataM(volind{1},volind{2},volind{3});
    volpos = [pos_neighbor,start_pos]-[volind{1}(1) volind{2}(1) volind{3}(1)]'+[1 1 1]';
    loopTracing(volpos,cutvol,2,2,23)
    fprintf('Check Lattice pos\n')
    pause
    %% ----extract primitive lat-vec
    neighborVec = pos_neighbor - repmat(start_pos',[findN 1])';
    Th = 2;
    nVec = abs(neighborVec);
    Lind = find((nVec(1,:)<Th & nVec(2,:)<Th ) | (nVec(2,:)<Th & nVec(3,:)<Th ) | (nVec(3,:)<Th & nVec(1,:)<Th ));
    SearchVecs_pixel_og = neighborVec(:,Lind)
    dt = delaunayTriangulation(SearchVecs_pixel_og');
    % calculate the volumn of the body
    [~, volume] = convexHull(dt);
    figure(2)
    clf
    scat3(SearchVecs_pixel_og,15,'b');
    hold on
    scatter3(0,0,0,25,'filled','MarkerEdgeColor','r')
    pause
end
%--------------------------------------------------------
function pair_arr = Simple_cubic_search_0109(pos_set, SearchVecs_pixel_half, SearchTh)
    t1 = clock;
Np = size(pos_set,1);
% Np=100;
pair_arr = zeros(Np,26);
[X,Y,Z] = meshgrid(-1:1);
cube_index = [X(:),Y(:),Z(:)];
cube_index(14,:)=[];
situation_overone = 0;
situation_non = 0;
for i1 = 1:Np
    if mod(i1,5000)==0
        fprintf('pos searched%d\n',i1)
    end
    currPos = pos_set(i1,:);
    for i2=1:26
        currSearchPos  = currPos+cube_index(i2,:)*SearchVecs_pixel_half;       %当前目标位置
         %-------DEBUG
         % for i3=1:26
         %         currSearchPos(i3,:)  = 0+cube_index(i3,:)*SearchVecs_pixel_half;
         % end
         % currSearchPos
         %-------DEBUG
%         pause    %-------DEBUG
        curr_dist      = pdist2(currSearchPos, pos_set);    %当前目标位置与各个atom_pos的距离
        ThInd = find(curr_dist < SearchTh);                    %距离阈值内的atom_pos的Ind
        %__________ 唯一 __________________________________________________
        if length( ThInd ) == 1
            pair_arr(i1,i2) = ThInd;
        %__________ 不唯一 ________________________________________________
        elseif length(ThInd)>1
            situation_overone = situation_overone+1;
            pair_arr(i1,i2) = -1;
        %__________ 距离阈值内没有待选原子 ________________________________
        elseif isempty(ThInd)
            situation_non = situation_non+1;
            pair_arr(i1,i2) = 0;
        end
     end
end
fprintf('empty:%d, overone:%d\n',situation_non, situation_overone)
    t2 = clock;
    consumeTime = round(etime(t2,t1)/60);
    fprintf('Time used: %dmin\n',consumeTime)

end

