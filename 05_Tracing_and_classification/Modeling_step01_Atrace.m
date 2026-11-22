%% Lattice Modeling Step02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1.Do A atom tracing [do interp for whole Vol] %%%
%%%                                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
here = fileparts(mfilename('fullpath'));
load(fullfile(here,'output','Model_output_00_localMax.mat'));
% DM is the masked reconstructed volume prepared in step 00.
% [Alat_final] [Amax_final] [Blat_final] [Bpos_Init] [polar_set_f]
% [RecFilename] [useA_idx] 

%% Do interpolation 
tic;interpVol = oversampling_n(DM,3,3); toc
whos interpVol
% interpVol_single = single(interpVol);
% tic;save('interpVol_single','interpVol_single','-v7.3');toc
% tic;load('interpVol_single');toc
% DataMatrix = interpVol;

%% ov3_pos 
Apos_ov3 = 3*Amax_final-2;
Bpos_Init_ov3 = 3*Bpos_Init-2; 

%% relocate A local max
aN = size(Apos_ov3,1);
Amax_new = zeros(size(Apos_ov3));
r=2;
for ii=1:aN
rc  = round(Apos_ov3(ii,:)); 
DBox = get_data_box(interpVol,rc(1),rc(2),rc(3),r); 
% loop(DBox)
% rc  = round(Amax_final(ii,:)); 
% DBox = get_data_box(DM,rc(1),rc(2),rc(3),2); 
% loop(DBox)
[~,midx] = max(DBox(:));
[nx,ny,nz] = ind2sub(size(DBox),midx);
Amax_new(ii,:) = rc+[nx,ny,nz]-r-1; 

% DBox = get_data_box(interpVol,Amax_new(ii,1),Amax_new(ii,2),Amax_new(ii,3),r); 
% loop(DBox)
% pause

end

%% Trace_A 
useVol = interpVol;
c_start = Amax_new;
[Atr_cut7,stats_Acut7] = func_atom_retrace_Poly_parallel_2402(useVol, c_start, 3, 4, 1); 
bad7p = find(Atr_cut7(:,1)==0);    bad7f = find(stats_Acut7.exitFlagArr~=0);
whos  bad7p bad7f 

    %%
    good7  = IndDiff(1:size(Amax_new,1),bad7p);
    figure(1318);subplot(121);
    scat3(Atr_cut7(good7,:));hold on
    scat3(c_start(bad7f,:),'r*')

    [out_rmsd, distarr] = prmsd2(Atr_cut7(good7,:),c_start(good7,:)) ;
    out_rmsd/3*0.3079
    subplot(122);histogram(distarr/3*0.3079);  

%% gaussian 补全
[Atr_cut7_G,stats_Acut7_G] = func_atom_retrace_G_2402(useVol, c_start(bad7p,:), 9); 
[out_rmsd, distarr] = prmsd2(Atr_cut7_G,c_start(bad7p,:));
% [out_rmsd, distarr] = prmsd2(B_cutest7_G,cut_start(bad7p_ogRec,:));
out_rmsd/3*0.3079
figure(13189);histogram(distarr/3*0.3079);

%% ---assemble 
%---assemble 
A_trace(good7,:) = Atr_cut7(good7,:);
A_trace(bad7p,:) = Atr_cut7_G;
figure(131123);scat3(A_trace);
% loopTracing((A_trace+2)/3,DM,2,2)

save(fullfile(here,'output','Model_output_01_Atrace.mat'), 'A_trace')




%% end



%% subfunction
%% subfunction
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
    fprintf('Check Lattice pos/n')
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
        fprintf('pos searched%d/n',i1)
    end
    currPos = pos_set(i1,:);
    for i2=1:26
        currSearchPos  = currPos+cube_index(i2,:)*SearchVecs_pixel_half;       %当前目标位置
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
fprintf('empty:%d, overone:%d/n',situation_non, situation_overone)
    t2 = clock;
    consumeTime = round(etime(t2,t1)/60);
    fprintf('Time used: %dmin/n',consumeTime)

end

%--------------------------------------------------------
function outputdist = p2ldist(twopos,p0)
    p1 = twopos(1,:);
    p2 = twopos(2,:);
    % v1=p0-p1; v2=p0-p2; 
    % norm(v1)*norm(v2)*sind(vang2(v1,v2))
    outputdist = norm(cross(p0-p1, p0-p2))./norm(p2-p1); 
end


%----------------------------------------------------------
function outputdist = p2fdist(facepos,uc_center,p0)
    vec1 = mean(facepos)-uc_center;
    vec1_n = vec1/norm(vec1);
    vec2 = p0-uc_center;
    outputdist = norm(vec1)-dot(vec1_n,vec2); 
end















