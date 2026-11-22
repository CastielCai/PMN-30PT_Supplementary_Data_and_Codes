function [err]= register_two_model_PMNPT_moveOri(para,orimodel,newmodel,str)

% newmodel=gpuArray(newmodel);
    al=para(1)/180*pi;
be=para(2)/180*pi;
ga=para(3)/180*pi;
ori=para(4:6);
rotal=[cos(al) -sin(al) 0;
    sin(al) cos(al) 0
    0 0 1]';
rotbe=[cos(be) 0 -sin(be);
    0 1 0;
    sin(be) 0 cos(be)]';
rotga=[1 0 0;
    0 cos(ga) -sin(ga);
    0 sin(ga) cos(ga);]';
mpos=mean(orimodel);
orimodel=(orimodel-mpos)*rotal*rotbe*rotga+ori+mpos;

% orimodel=gpuArray(orimodel);
% ori=para(1:3);
% orimodel=orimodel+ori;
% tmodel=newmodel-mpos;

% gModel=gpuArray(orimodel);
% gnewModel=gpuArray(tmodel);
% Ymax=max(gModel(:,2));
% Ymin=min(gModel(:,2));

% slicenum=5;
% Yslice=linspace(Ymin,Ymax,slicenum+1);
% dmin=zeros(size(gModel, 1),1);
% dmin=dmin+1000;
% % slicenum
% % tic
% for i=1:slicenum
%     ind=gModel(:,2)>=Yslice(i) & gModel(:,2)<=Yslice(i+1);
%     gModel_s=gModel(ind,:);
%      ind2=gnewModel(:,2)>Yslice(i)-10 & gnewModel(:,2)<Yslice(i+1)+10;
%     gnewModel_s=gnewModel(ind2,:);
%     difpd=pdist2(gModel_s,gnewModel_s);
%     [v,~]=min(difpd,[],2);
%     dmin(ind)=min(dmin(ind),v);
% %     difpd=pdist2(orimodel,tmodel)
% end
% toc

 [~,D]=knnsearch(newmodel,orimodel);
dmin=D;
% % difpd=pdist2(orimodel,tmodel);
% [k,dist]=dsearchn(orimodel,tmodel);
% [v,~]=min(difpd);

% err=sum(dmin);
% dmin(dmin==1000)=0;
if strcmp(str,"rms")
    err=rms(dmin);
elseif strcmp(str,"sum")
    err=sum(dmin);
end
end