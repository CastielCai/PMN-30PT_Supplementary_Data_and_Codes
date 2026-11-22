function interpVol = oversampling_n(inputVol,n,dvn)
% n:oversampkking ratio
% n = 3;
% dvn = 6; %% how many parts to do oversampling sperate
%%
ddim = 2; %% divide Dim: do not change（for now）
halflap = 2;
Voldim = size(inputVol);
dL = round(Voldim(ddim)/dvn);
sc = zeros(dvn,2);  % sc: slice coords
ac = zeros(dvn,2);  % ac: assemble coords
cutVol_ov1 = cell(1,3);
cutVol_ov3 = cell(1,3);
for ii=1:dvn
    if ii==1
        sc(ii,:) = [1 ii*dL+halflap];
        % ac(ii,:) = [1
    elseif ii==dvn
        sc(ii,:) = [(ii-1)*dL-halflap+1 Voldim(ddim)];
        % ac(ii,:) = [1
    else
        sc(ii,:) = [(ii-1)*dL-halflap+1 ii*dL+halflap];
        % ac(ii,:) = [1
    end
    cutVol_ov1{ii} = inputVol(:,sc(ii,1):sc(ii,2),:);
end

for jj=1:dvn
    c_vol = cutVol_ov1{jj};
 %----interp3
    xx = (1:size(c_vol,1)) - round((size(c_vol,1)+1)/2);
    yy = (1:size(c_vol,2)) - round((size(c_vol,2)+1)/2);
    zz = (1:size(c_vol,3)) - round((size(c_vol,3)+1)/2); 
    xxi = ((n*xx(1)):(xx(end)*n))/n;
    yyi = ((n*yy(1)):(yy(end)*n))/n;
    zzi = ((n*zz(1)):(zz(end)*n))/n; 
    [Y,X,Z] = meshgrid(yy,xx,zz);
    [Yi,Xi,Zi] = meshgrid(yyi,xxi,zzi); 
    cutVol_ov3{jj} = interp3(Y,X,Z,c_vol,Yi,Xi,Zi,'linear',0); 
end
% interpVol=cutVol_ov3; 
% pitching cutvols
interpVol = [];
for kk=1:dvn 
    if kk==1 
        interpVol = cat(ddim,interpVol,cutVol_ov3{kk}(:,1 : end-n*halflap,:)); 
    elseif kk==dvn 
        interpVol = cat(ddim,interpVol,cutVol_ov3{kk}(:,n*(halflap-1)+2 : end,:));  
    else 
        interpVol = cat(ddim,interpVol,cutVol_ov3{kk}(:,n*(halflap-1)+2 : end-n*halflap,:));  
    end 
%     interpVol(:,ac(kk,1):ac(kk,2),:) = cutVol_ov3{kk};
end
interpVol(interpVol<0)=0;



end


