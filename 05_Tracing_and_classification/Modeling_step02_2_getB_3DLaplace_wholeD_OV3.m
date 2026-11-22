%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a part of 
% 'Lattice_modeling_step02_nogroup_wholeOV_Btest3_wholeData.m'
% After A trace part
% with Gaussian method replaced By [project center sphere to 2D to locat atom]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
here = fileparts(mfilename('fullpath'));
load(fullfile(here,'output','Model_output_01_Atrace.mat')) 

%% Get polar set variable 
load(fullfile(here,'output','Model_output_00_localMax.mat'), 'polar_set_f')
setfull = polar_set_f;
setN = length(setfull);
polar_lat = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
edgepair  = [1 2; 2 3; 3 4; 4 1; 1 5; 2 6; 3 7; 4 8; 5 6; 6 7; 7 8; 8 5];
facepair  = [1 2 3 4; 5 6 7 8; 1 2 6 5; 4 3 7 8; 1 4 8 5; 2 3 7 6];
setinfo.polar_lat = polar_lat;
setinfo.edgepair  = edgepair;
setinfo.facepair  = facepair;
  

%% find B local max
%----- OTHER parameter 
distanceMode = 1; % 0 for edge, 1 for face 
cen_arr = zeros(setN,3);
emptyInd0 =[];
iter=0;
 
%% find B local max (project to 2D method )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  may need multiple time run  %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%---- Box parameter & Gaussian
boxD = 51;
[Boxcoordinates, BoxCenter, BoxR, ~] = create_box(boxD);
docheck=0;
% docheck=1;
  % iter = 0;
  % iter = 1;
%   iter = 2;
% close all
tic
Dth_ratio = 0.4;
emptyInd=[];

% count how many time whole iteration has been run
iter=iter+1; 
searched=[];
if iter>1
    for ii=1:iter-1
        searched = cat(2,searched,searchInd{ii});
    end 
end 
ts = IndDiff(1:setN,searched); % ind 2 be search
% searchInd{3} = moveto2;
ts = unique(cat(2,moveto2,searchInd{2}));

fprintf('current iter:%d\n',iter)
fprintf('number of UC to be search:%d\n\n',length(ts))
%-------------Laplace mask
mask = zeros(boxD,boxD,boxD);
mask_r = (boxD-1)/2;
[hx,hy,hz] = meshgrid(-mask_r:mask_r,-mask_r:mask_r,-mask_r:mask_r);
pooos = [hx(:),hy(:),hz(:)];
for ii=1:size(pooos,1)
    poolength(ii) = norm(pooos(ii,:)); 
end
%-------------
rtN = 9;
for rt=1:rtN
rad = 11+rt; 
mask_ind = find(poolength<=rad);
mask_arr{rt} = mask;
mask_arr{rt}(mask_ind) = true;
end
% trace_meth = zeros(length(ts),rtN); 
 
for ii=1:length(ts)
    %%
    %------ find local max  
    if iter>=2
    disp(ii)
    end
    c_8A = A_trace(setfull(ts(ii),2:9),:);
    cen_arr(ts(ii),:) = mean(c_8A);
    rc  = round(cen_arr(ts(ii),:)); 
    DataBox1 = get_data_box(interpVol,rc(1),rc(2),rc(3),BoxR); 
    %%%%%%-----------------------------------------------------------------
    if iter==1
        useBox = DataBox1; 
        useBox(useBox<0)=0; 
        %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
        %'for pp=1',just for fold
      for pp=1
        %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
        se           = strel3d(2); % whos se
        dilatedBW    = imdilate(useBox, se);
        maxInd       = find( useBox == dilatedBW & useBox > 10 );
        % maxVals      = useBox(maxInd);
        [mx,my,mz] = ind2sub(size(useBox),maxInd);
        maxincube = [mx,my,mz];
        incube_dist{ts(ii)} = pdist2([mx,my,mz],[BoxR+1,BoxR+1,BoxR+1]);
        % loopTracing([mx,my,mz]', useBox,3,3)
        maxpos = rc + [mx,my,mz] -(boxD+1)/2; 
        %------ distance & Intensity filter
        cendist_arr = pdist2(maxpos,cen_arr(ts(ii),:));
        ucdist_arr  = zeros(length(mx),1);  % distance relate to unit cell
        %----get distance to cube boder------------------------------------
        for jj=1:length(mx)
            cmpos = maxpos(jj,:); 
            if distanceMode==0
                dn=12;  c_ucd = zeros(dn,1);  UCD_cen = zeros(dn,1);
            elseif distanceMode==1
                dn=6;   c_ucd = zeros(dn,1);  UCD_cen = zeros(dn,1);
            end 
            for kk=1:dn 
                if distanceMode==0
                    c_ucd(kk) = p2ldist(c_8A(setinfo.edgepair(kk,:),:),cmpos);
                    UCD_cen(kk) = p2fdist(c_8A(setinfo.edgepair(kk,:),:),cen_arr(ts(ii),:),cen_arr(ts(ii),:));
                elseif distanceMode==1 
                    c_ucd(kk) = p2fdist(c_8A(setinfo.facepair(kk,:),:),cen_arr(ts(ii),:),cmpos);
                    UCD_cen(kk) = p2fdist(c_8A(setinfo.facepair(kk,:),:),cen_arr(ts(ii),:),cen_arr(ts(ii),:));
                end 
            end
            ucdist_arr(jj) = min(c_ucd);
        end

        D_arr = cat(2,cendist_arr,ucdist_arr);
        %----- dist validation
        UCDth(ts(ii)) = mean(UCD_cen);
        % asdasd = find(ucdist_arr>UCDth(ts(ii))*Dth_ratio & cendist_arr<15);
        asdasd = find(cendist_arr<15);
        % inner_N(ts(ii)) = length(asdasd);
        inner_intense = useBox(maxInd(asdasd));
        enough_intense_ind = find(inner_intense>max(inner_intense)/2); 
        if ~isempty(enough_intense_ind) 
            % figure(9876);clf;scat3(c_8A);hold on;scat3(maxpos);
            [~,ci] = min(cendist_arr(asdasd(enough_intense_ind)));
            B_start(ts(ii),:) = maxpos(asdasd(enough_intense_ind(ci)),:);
            B_dist(ts(ii),1) = cendist_arr(asdasd(enough_intense_ind(ci))); 
            B_dist(ts(ii),2) = ucdist_arr(asdasd(enough_intense_ind(ci))); 
        else
            emptyInd(end+1) = ts(ii);
        end  
        %%%%%%%%%%%%%%%%%%%%
        if docheck && any(ts(ii)==emptyInd0)
            % %
            disp(ts(ii))
            figure(9876);clf
            scat3([mx,my,mz]);hold on
            scat3(c_8A-rc+(boxD+1)/2);
            % scat3(maxpos(asdasd,:),'*r');  
            % figure(9876);clf;scat3(c_8A);hold on;
            scat3([(boxD+1)/2,(boxD+1)/2,(boxD+1)/2],'+r');
            view([1 8 0.7])
            % useBox(mx(ci),my(ci),mz(ci))
            % DataBox1(mx(ci),my(ci),mz(ci))
            % loopTracing([mx(ci),my(ci),mz(ci)],useBox,3,3)
            loopTracing([mx,my,mz],DataBox1,2,3); 
            loopTracing([mx,my,mz],useBox,2,3); 
            % maxincube(asdasd(enough_intense_ind(ci)),:)
            % loop(useBox)
            pause
            close all
        end
        %%%%%%%%%%%%%%%%%%%%%%%   
      end 
        %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    elseif iter>=2
        c_8A_incube = c_8A-rc+(boxD+1)/2;
        cen_incube  = mean(c_8A_incube); 
        if distanceMode==0
            dn=12;  c_ucd = zeros(dn,1);  UCD_cen = zeros(dn,1);
        elseif distanceMode==1
            dn=6;   c_ucd = zeros(dn,1);  UCD_cen = zeros(dn,1);
        end 
        parfor rt=1:rtN 
            c_mask = mask_arr{rt}; 
            J = regionfill3D(DataBox1,c_mask) ;
            useBox = DataBox1-J;
            useBox(useBox<0)=0; % only find pos in mask
        %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
        %apply local max search & tracing
        %'for pp=1',just for fold
        for pp=1
            se           = strel3d(2); % whos se
            dilatedBW    = imdilate(useBox, se);
            maxInd       = find( useBox == dilatedBW & useBox > 10 );
            % maxVals      = useBox(maxInd);
            [mx,my,mz] = ind2sub(size(useBox),maxInd);
            maxincube = [mx,my,mz];
            % incube_dist_iter2{ii} = pdist2([mx,my,mz],[BoxR+1,BoxR+1,BoxR+1]); 
            maxpos = maxincube; 
            % loopTracing(maxpos,useBox,3,2)
            %------ distance & Intensity filter
            cendist_arr = pdist2(maxpos,cen_incube);
            ucdist_arr  = zeros(length(mx),1);  % distance relate to unit cell 
            %========= spatial constrain
            asdasd = find(cendist_arr<15);
            % inner_N(ts(ii)) = length(asdasd);
            inner_intense = useBox(maxInd(asdasd));
            enough_intense_ind = find(inner_intense>max(inner_intense)/2);
            if ~isempty(enough_intense_ind) 
                % figure(9876);clf;scat3(c_8A);hold on;scat3(maxpos);
                [~,ci] = min(cendist_arr(asdasd(enough_intense_ind)));
                c_Bstart = maxpos(asdasd(enough_intense_ind(ci)),:); 
                %------------------apply trace:-------------------------
                B_locM(ii,:,rt) = c_Bstart;
                % loopTracing(c_Bstart,useBox,3,2)
                O_Ratio = 1; % interpolation done outside
                [c_Btrace,c_stats] = func_atom_retrace_Poly_parallel_2402(useBox, c_Bstart, 3, 4, 1); 
                B_Trace_poly(ii,:,rt) = c_Btrace; 
                B_Trace_flag(ii,rt) = c_stats.exitFlagArr;  
                if c_stats.exitFlagArr==0
                    trace_meth(ii,rt) = 1; 
                elseif c_stats.exitFlagArr~=0
                    [c_Btrace,c_stats] = func_atom_retrace_G_2402(useBox, c_Bstart, 3, 4); 
                    trace_meth(ii,rt) = 2;
                    B_Trace_Gau(ii,:,rt) = c_Btrace; 
                end
            end
        end
        %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
        end  
    end
 

end
emptyInd = unique(emptyInd);
% if iter==1
%     % searchInd{iter} = IndDiff(1:setN,emptyInd);
% else
%     % searchInd{iter} = IndDiff(emptyInd0,emptyInd); 
% end
emptyInd0 = emptyInd;
% close all
fprintf('%d empty UC\n',length(emptyInd0));
if iter==1
    inner_intrisic_lm = Indcommon(searchInd{1},innerInd);
    fprintf('Inner good UC:%d (%f). \n',[length(inner_intrisic_lm),length(inner_intrisic_lm)/length(innerInd)]);
end
% docheck=0;
% docheck=1;
% empInd = emptyInd0;
toc


%% Determine B pos (1) : trace iter1 local max
%% Determine B pos (1) : trace iter1 local max
%% Determine B pos (1) : trace iter1 local max 
cIdx = searchInd{1};
c_start = B_start(cIdx,:) ;%*3-2; 
% loopTracing(c_start,interpVol,3,3)
O_Ratio = 1; % interpolation done outside
[posB_cut7_P,statsB_cut7_P] = func_atom_retrace_Poly_parallel_2402(interpVol, c_start, 3, 4, O_Ratio); 
bad7p = find(posB_cut7_P(:,1)==0);     bad7f = find(statsB_cut7_P.exitFlagArr~=0);
whos  bad7f bad7p
good7 = IndDiff(1:size(c_start,1),bad7p);

%%
%----- gaussian 补全
[posB_cut7_G,statsB_cut7_G] = func_atom_retrace_G_2402(interpVol, c_start(bad7p,:), 7);
% figure;scat3(posB_cut7_G-c_start(bad7p,:)); 


% ----assemble
B_trace(cIdx(good7),:) = posB_cut7_P(good7,:);
B_trace(cIdx(bad7p),:) = posB_cut7_G;
  



%%
% B_trace(35361,3)=0;
figure(1541);clf
% scat3(B_trace(cIdx(good7),:)-c_start(good7,:));
% scat3(B_trace(cIdx02,:));
% scat3(A_trace,'r.');hold on
scat3(B_trace,'b.');



%% Determine B pos (2) : Laplacian
%% Determine B pos (2) : Laplacian
%% Determine B pos (2) : Laplacian 
%% statistic
    whos B_locM B_Trace trace_meth
    % find out all_zero UC(Laplace didnt worn)
    for c_r = 1:rtN
    % c_r = 1;  % real r=c_r+14
    c_method = trace_meth(:,c_r);
    fprintf('mask_R=%d, poly:%d, Gau:%d, empty:%d\n',[c_r+11 nnz(c_method==1) nnz(c_method==2) nnz(c_method==0)])
    goodtrace{c_r} = find(c_method~=0);
    end
    % goodtrace
    

%% main part use mask_R==15;

% searchInd{3} = moveto2;
% cIdx02 = unique(cat(2,searchInd{2},searchInd{3}));
cIdx02 = ts;

cpBN = length(cIdx02); % current part B number
c_part_B_final = zeros(cpBN,3);
clc
% R=15
c_method = trace_meth(:,1);
Round{1} = find(c_method>0);
left = IndDiff(1:cpBN,Round{1});
whos left
%--- 
% % look down(R=16  14  13)
% c_method = trace_meth(:,3);
% c_good =  find(c_method>0);
% round{2}  = Indcommon(left,c_good);
% left = IndDiff(left,round{2});
% whos left 

%% transition from in-cube to in-Volume
use_rt = 4;
c_method = trace_meth(:,use_rt);
c_Btr_poly = B_Trace_poly(:,:,use_rt);
c_Btr_Gau = B_Trace_Gau(:,:,use_rt);
c_part_B_final_incube = c_part_B_final;
c_part_B_final_incube(c_method==1,:) = c_Btr_poly(c_method==1,:);
c_part_B_final_incube(c_method==2,:) = c_Btr_Gau(c_method==2,:);
left = find(c_method==0); 
figure;scat3(c_part_B_final_incube);


%% assamble
    cIdxset = Round{1}; % belong to cIdx02 
    for ii=1:length(c_method)  
        c_8A = A_trace(setfull(cIdx02(ii),2:9),:);
        cen_arr = mean(c_8A);
        rc  = round(cen_arr); 
        if any(ii==cIdxset)
            c_part_B_final(ii,:) = ...
                c_part_B_final_incube(ii,:)+rc-(boxD+1)/2;
        elseif any(ii==left)
            c_part_B_final(ii,:) = rc; 
        end
    
    end 
    B_trace(cIdx02,:) = c_part_B_final;


    %%
    c_BIdx_nonempty = IndDiff(1:setN,cIdx02(left));
    figure;
    subplot(121); 
    scat3(Blat_final(c_BIdx_nonempty,:));
    subplot(122); 
    scat3(Blat_final(cIdx02(left),:));
    
    figure(2322);clf  
    scat3(B_trace,'b.'); hold on 
    scat3(B_trace(oddY,:),'ro');
    
    
    c_part_B_final(left,:);
    B_trace(B_leftIdx,:);
    
    
    B_leftIdx = cIdx02(left);


%% check with A center
for ii=1:size(setfull,1)
    Acen_arr(ii,:) = ...
            mean(A_trace(setfull(ii,2:9),:));
end

%% 
% figure(1729);
% scat3(B_trace-Acen_arr);
ABdif = B_trace-Acen_arr;
for ii=1:size(ABdif,1)
    p_arr(ii) = norm(ABdif(ii,:));
end
figure(2139);clf
%histogram(p_arr)
histogram(abs(ABdif(:,2)))

%%
% oddX = find(abs(ABdif(:,1))>8);
oddY = find(abs(ABdif(:,2))>8);
% oddZ = find(abs(ABdif(:,3))>8);
% whos oddX oddY oddZ
whos oddY 
moveto2 = Indcommon(oddY,searchInd{1});
size(Indcommon(oddY,searchInd{2}))
new_ts = unique(cat(2,searchInd{2},moveto2));
for ii=1:size(moveto2,2)
    rc = round(Acen_arr(moveto2(ii),:)); 
    DataBox1 = get_data_box(interpVol,rc(1),rc(2),rc(3),BoxR);  
    useBox = DataBox1; 
    useBox(useBox<0)=0; 
    loop(useBox)
    pause
end
  
% loopTracing((B_trace(oddY,:)+2)/3,DM,1,2)


%% replace oddY B pos
B_tr_rep = B_trace;
for rr=1:size(oddY)
clat =  Blat_final(oddY(rr),:);
[stpidx,croundidx] = find_neibposInd_Include_D(Blat_final,1.1,clat);
croundidx = IndDiff(croundidx,stpidx);
B_tr_rep(stpidx,:) = mean(B_trace(croundidx,:)); 
end

    figure(2322);clf  
    scat3(A_trace,'b.'); hold on 
    scat3(B_tr_rep,'r.');
    
 


%% SAVE Final model
%% SAVE Final model
    % pos
    % [B_trace]
    
    % pos properties
    cIdx01 = searchInd{1};
    B_prop = zeros(setN,2);
    % colom 01, Local max type: [1]-Intrinsic   [2]-Laplacian  [0]-empty
    B_prop(cIdx01,1) = 1;  % 
    B_prop(cIdx02(Round{1}),1) = 2;  % 
     
    % colom 02, trace type:     [1]-Poly   [2]-Gaussian  [0]-empty
    B_prop(cIdx01(good7),2) = 1;  % 
    B_prop(cIdx01(bad7p),2) = 2;  %  
    
    B_prop(cIdx02(c_method==1),2) = 1;  % 
    B_prop(cIdx02(c_method==2),2) = 2;  % 
    
    %---Check
    idxemp1 = find(B_prop(:,1)==0);
    idxemp2 = find(B_prop(:,2)==0);
    IndDiff(idxemp1,idxemp2)
    

    bp11 = nnz(B_prop(:,1)==1 & B_prop(:,2)==1);
    bp12 = nnz(B_prop(:,1)==1 & B_prop(:,2)==2);
    bp21 = nnz(B_prop(:,1)==2 & B_prop(:,2)==1);
    bp22 = nnz(B_prop(:,1)==2 & B_prop(:,2)==2);
    clc
    
    fprintf('Intrinsic ratio:%f%% \n',nnz(B_prop(:,1)==1)/setN*100)
    fprintf('Laplacian ratio:%f%% \n',nnz(B_prop(:,1)==2)/setN*100)
    fprintf('empty ratio:%f%% \n\n',length(B_leftIdx)/setN*100)
    
    fprintf('Intrinsic poly ratio:%f%% \n',bp11/(bp11+bp12)*100)
    fprintf('laplacian poly ratio:%f%% \n',bp21/length(cIdx02)*100)
    fprintf('Total poly ratio:%f%% \n',(bp11+bp21)/setN*100)
    fprintf('Total Gaussian ratio:%f%% \n',(bp12+bp22)/setN*100)
    
   
    %% 
    figure; hold on
    p1=scat3(Blat_final(c_BIdx_nonempty,:)); 
    p1.MarkerEdgeAlpha = 0.1;
    scat3(Blat_final(cIdx02(left),:));
    
    %% save final model
    figure;hold on; scat3(A_trace,'.');scat3(B_tr_rep,'.');
 
    % [Alat_final] [Blat_final] [A_trace] [B_trace] [B_prop]
    % [polar_set]  [setinfo]
    % discription
    discription{1} = ' Blat_final & B_trace are connected';
    discription{2} = ' B oddY is filled by mean position of its neighbour pos';

save(fullfile(fileparts(here),'06_final_coordinates','Model_output_02_Final_model.mat'), 'Alat_final', 'Blat_final',...
    'A_trace', 'B_tr_rep','B_prop', 'polar_set_f', 'setinfo', 'discription')
 

    
















