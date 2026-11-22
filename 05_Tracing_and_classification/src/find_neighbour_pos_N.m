function [start_pos,pos_neighbor] = find_neighbour_pos_N(atom_pos,N,varargin)
% Find Neibourpos with asigned Numbers
if nargin==2
    StartMedian = [median(atom_pos(1,:));median(atom_pos(2,:));median(atom_pos(3,:))];
elseif nargin==3
    StartMedian = varargin{1};
end
[~,StartPosInd] = min( sqrt( sum( (atom_pos - repmat(StartMedian,[1 size(atom_pos,2)])).^2,1) ) ); 
start_pos = atom_pos(:,StartPosInd);
Dist = pdist2(atom_pos(:,StartPosInd)',atom_pos');
Dist_sort = sort(Dist,'ascend');
Nearnum = N; %ÕÒN¸ö×î½üÁÚ
ind_neighbor = zeros(1,Nearnum);
for i=1:Nearnum
    ind_neighbor(i) = find(Dist==Dist_sort(i+1));
end
pos_neighbor = atom_pos(:,ind_neighbor);
% deletind = find(pos_neighbor(3,:)<start_pos(3)+2 & pos_neighbor(3,:)>start_pos(3)-2); 
% pos_neighbor(:,deletind) = [];
% while length(pos_neighbor)>8
%     [~,deletind] = max(abs(pos_neighbor(3,:)-start_pos(3)));
%     pos_neighbor(:,deletind) = [];
% end
end

