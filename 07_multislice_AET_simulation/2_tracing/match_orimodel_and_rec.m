function [err] = match_orimodel_and_rec(para,FinalVol,model)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
mr=round(model+para);
ind=sub2ind(size(FinalVol),mr(:,1),mr(:,2),mr(:,3));
val=FinalVol(ind);
err=-sum(val);

end

