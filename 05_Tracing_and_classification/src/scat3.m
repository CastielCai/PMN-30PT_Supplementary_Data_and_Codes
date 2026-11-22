function scat3(pos,varargin)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
%check input data
if ~any(size(pos)==3) || length(size(pos))~=2
   error('Input Data format ERROR!!') 
end
% permute if needed
if size(pos,2)~=3
   pos = permute(pos,[2 1]);
end
dotsize = 15;
color   = 'b';
if nargin==2
    dotsize = varargin{1};
    color = 'b';
elseif nargin==3
    dotsize = varargin{1};
    color   = varargin{2};
end
scatter3(pos(:,1),pos(:,2),pos(:,3),dotsize,color,'filled','MarkerEdgeColor','k')
axis equal
end

