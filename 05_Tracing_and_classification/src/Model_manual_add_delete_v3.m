function Model_manual_add_delete_v3(model, DataMatrix, Devide, savefileName,Tempinfo)

% crop DataMatrix
crop_shift = [Devide.X(1),Devide.Y(1),Devide.Z(1)]'-[1,1,1]';
cdM = DataMatrix(Devide.X,Devide.Y,Devide.Z); % [crop_DataM];
% get inner atompos
cropMask3D = zeros(size(DataMatrix));
cropMask3D(Devide.X,Devide.Y,Devide.Z)=1;
inner_Ind = find_inner_pos_ind(model,cropMask3D);
inner_pos = model(:,inner_Ind);
IPS = inner_pos-crop_shift;          % [inner_pos_shifted], figure(333);scat3(IPS)
aType = ones(1,length(model));   % Assume pos all Valid at begaining
if isempty(Tempinfo.addpos)
    addposTP=[];
else
    addposTP=Tempinfo.addpos-crop_shift;
end
if ~isempty(Tempinfo.deleteInd)
    aType(Tempinfo.deleteInd)=0;
end
IT = aType(inner_Ind);               % [inner_type], Get initial inner pos type
% figure(332);loopTracing(IPS,cdM,2);pause

%--------------------------------------------------------------------------
% 1st:
% [inner_pos]=[model](inner_Ind); %for [inner_pos_shifted]
% [inner_type]=[aType](inner_Ind);
% 2nd:
% [Layer_pos] = [inner_pos](Layer_ind);
% [Layer_type] = [inner_type](Layer_ind);
%--------------------------------------------------------------------------

thresh = 3;
dotS=  3;
color = ['r'; 'y'];
quit = 0;
click = 0;
doloop = 1;
%--------------------------------------------------------------------------
Direction = input('输入loop方向：1、2、3 == X、Y、Z\n');
switch Direction
    case 1
    cdMn = permute(cdM,[2,3,1]);
    IPSn = IPS([2,3,1],:);
    case 2
    cdMn = permute(cdM,[1,3,2]);
    IPSn = IPS([1,3,2],:);
    case 3
    cdMn = cdM;
    IPSn = IPS;
end
num = size(cdMn,3);
figure(211);clf
k=1; jj=1;
L_ind = abs(IPSn(3,:)-jj)<thresh;  % In-Layer pos ind
L_posP = IPSn([1,2],L_ind);        % In-Layer pos to plot
L_t = IT(L_ind);                   % In-Layer pos type
set(gcf,'KeyPressFcn',@callback)
    
function callback(~,event)
    switch event.Key
    case 'a'
        k = k-1; doloop=1;
    case 'd'
        k = k+1; doloop=1;
    case 'space'
        click=1; doloop=0;
    case 'escape'  % 即时保存
        click=0; doloop=0;
%         quit=1;
        if ~isempty(addposTP)
            addpos = addposTP+crop_shift;
        else 
            addpos = [];
        end
        aType(inner_Ind) = IT;
        deleteInd = find(aType==0);
        save(['modify_' savefileName], 'deleteInd', 'addpos')
%         load(['modify_' savefileName '.mat'])
        fprintf('modify saved\n')
%         close
    end
    %---- apply modify ---------------------------------------------
    if click==1 && doloop==0
    p=round(ginput(1));
    if ~isempty(p)
        [p_dist01,curr_ind01] = find_pos(p,L_posP,5); % 检查输入点集
        [p_pos,p_dist02,curr_ind02] = find_addpos(p,addposTP,5,Direction); % 检查add点集
        if  ~isempty(p_dist01) && isempty(p_dist02)
            situ = 1;
        elseif ~isempty(p_dist01)  && ~isempty(p_dist02)
            if p_dist01<=p_dist02
                situ = 1;
            else
                situ = 2;
            end
        elseif isempty(p_dist01) && ~isempty(p_dist02)
            situ = 2;
        elseif isempty(p_dist01) && isempty(p_dist02)
            situ = 3;
        end
        switch situ
        case 1    % 删除/保留 输入model的原子
            L_t(curr_ind01) = mod(L_t(curr_ind01)+1,2);
            IT(L_ind(curr_ind01)) = L_t(curr_ind01);
            plot(L_posP(2,curr_ind01),L_posP(1,curr_ind01),'o','MarkerSize',dotS,...
            'MarkerEdgeColor','r','MarkerFaceColor',color(L_t(curr_ind01)+1));
            fprintf('Curr model modifyed\n')
        case 2    % 从add集删去
            addposTP(:,curr_ind02)=[];fprintf('added pos deleted\n')
            plot(p_pos(2),p_pos(1),'*','MarkerSize',dotS+4,...
            'MarkerEdgeColor','r','MarkerFaceColor','r');
        case 3    % 添加至add集
            switch Direction
            case 1
            search_pos = [jj,   p(2), p(1)];
            case 2
            search_pos = [p(2), jj,   p(1)];
            case 3
            search_pos = [p(2), p(1), jj  ];
            end
            Th = 3;
            curr_addpos = find_maxima(search_pos,cdM,Th);
            addposTP(:,end+1) = curr_addpos;
            fprintf('pos added\n')
            switch Direction
            case 1
            p_pos = curr_addpos([2,3]);
            case 2
            p_pos = curr_addpos([1,3]);
            case 3
            p_pos = curr_addpos([1,2]);
            end
            plot(p_pos(2),p_pos(1),'*','MarkerSize',dotS+4,...
                 'MarkerEdgeColor','k','MarkerFaceColor','k');
        end
    end
    doloop = 0;
    click  = 0;
    end
    %------- loop to observe --------
    if doloop==1
        jj=mod(k,num);
        if jj==0
          jj=num; 
        end
        % get index & position in current Layer
        L_ind = find(abs(IPSn(3,:)-jj)<thresh);  % In-Layer pos ind
        L_posP = IPSn([1,2],L_ind);        % In-Layer pos to plot
        L_t = IT(L_ind);                   % In-Layer pos type
        %------- plot model pos--------
        clf
        imagesc(squeeze(cdMn(:,:,jj)));hold on; 
        colormap gray
        if ~isempty(L_ind)
            fprintf('Layer:%d, vali:%d, unvali:%d\n',jj, sum(L_t==1),sum(L_t==0))
        end
        plot(L_posP(2,L_t==1),L_posP(1,L_t==1),'o','MarkerSize',dotS,'MarkerEdgeColor','r','MarkerFaceColor','y'); %
        plot(L_posP(2,L_t==0),L_posP(1,L_t==0),'o','MarkerSize',dotS,'MarkerEdgeColor','r','MarkerFaceColor','r'); %
        %------- plot addposTP--------
        if ~isempty(addposTP)
            plotaind = find(abs(addposTP(Direction,:)-jj)<thresh);
            if ~isempty(plotaind)
                switch Direction
                case 1
                    plotapos = addposTP([2,3],plotaind);
                case 2
                    plotapos = addposTP([1,3],plotaind);
                case 3
                    plotapos = addposTP([1,2],plotaind);
                end
                plot(plotapos(2,:),plotapos(1,:),'*','MarkerSize',dotS,...
                     'MarkerEdgeColor','k'); %
            end
        end
        axis equal; title(sprintf('%d/%d',jj,num))
    end
end









%% SubfunctionS
%--------------------------------------------------------------------------
    function Ind = find_inner_pos_ind(pos,curr_support)
        if  any(curr_support(:)~=0 & curr_support(:)~=1)
            error('support matrix should only contain 0s and 1s!!!!')
        end
        if size(pos,1)~=3
            error('first demention of pos matrix must be 3!!!!')
        end
        N = size(pos,2);
        r_pos = round(pos);
        Ind = [];
        for ii=1:N
            try
                value = curr_support(r_pos(1,ii),r_pos(2,ii),r_pos(3,ii));
            catch
                value=0;
            end
            if value==1
                Ind(end+1)=ii;
            end
        end
    end
%--------------------------------------------------------------------------
    function [p_dist,ind] = find_pos(p,pset,Th)
        ind=[];
        p_dist=[];
        if ~isempty(pset)
            stp = [p(2),p(1)];   % start pos
            posdist = pdist2(stp,pset');
            posdist_sort = sort(posdist,'ascend');
            if ~isempty(posdist_sort) && posdist_sort(1)<=Th
                p_dist=posdist_sort(1);
                ind = find(posdist==posdist_sort(1));
            end
        end
    end
    
%--------------------------------------------------------------------------
    function [p_pos,p_dist,ind] = find_addpos(p,pset,Th,Direction)
        p_pos=[];
        p_dist=[];
        ind=[];
        if ~isempty(pset)
    	    cind = find(abs(pset(Direction,:)-jj)<thresh);
            switch Direction
            case 1
            pset_2 = pset([3,2],cind);
            case 2
            pset_2 = pset([1,3],cind);
            case 3
            pset_2 = pset([1,2],cind);
            end
            if ~isempty(pset_2)
                stp = [p(2),p(1)];   % start pos
                posdist = pdist2(stp,pset_2');
                posdist_sort = sort(posdist,'ascend');
                if ~isempty(posdist_sort) && posdist_sort(1)<=Th
                    p2ind = find(posdist==posdist_sort(1));
                    p_pos = (pset_2(:,p2ind));
                    p_dist = posdist_sort(1);
                    ind = cind(p2ind);
                else
                    ind=[];
                end
            end
        end
    end

%--------------------------------------------------------------------------
    function maxpos = find_maxima(search_pos,Datamatrix,Th)
        sp = search_pos;
        boxInd = -Th:Th;
        box = Datamatrix(sp(1)+boxInd, sp(2)+boxInd, sp(3)+boxInd);
        [~,maxind] = max(box(:));
        [x,y,z]= ind2sub(size(box),maxind);
        maxpos = [sp(1)+boxInd(x);...
                  sp(2)+boxInd(y);...
                  sp(3)+boxInd(z)];
    end

end