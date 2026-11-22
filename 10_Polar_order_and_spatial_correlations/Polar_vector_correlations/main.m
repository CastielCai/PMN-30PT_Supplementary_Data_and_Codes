% The baseline F(r) and its correlation length are in 09/Gaussian_kernel_sensitivity.
% This step uses the shared field (08) and experimental B positions (06)
% to obtain the angular distributions and their Gaussian and gamma peaks.
angularCorrelation(fileparts(mfilename('fullpath')));

function angularCorrelation(here)
root=fileparts(fileparts(here));
f=load(fullfile(root,'08_polar&domain&POP','output','polar_field.mat'),'vectors','pixel_A');
m=load(fullfile(root,'06_final_coordinates','Model_output_02_Final_model.mat'),'B_tr_rep');
radial=load(fullfile(root,'09_Stability_analysis','Gaussian_kernel_sensitivity','output','gaussian_comparison.mat'),'xi_A');
pos=(double(m.B_tr_rep)+2)/3*f.pixel_A; u=f.vectors./vecnorm(f.vectors,2,2);
edges_A=linspace(0,120,601); angle_edges_deg=0:180;
nShell=600; bins=180; histCount=zeros(nShell,bins);
sumAngle=zeros(nShell,1); sumSquare=sumAngle; sumLog=sumAngle;
positiveCount=sumAngle; sumPositive=sumAngle;
for i=1:32:size(pos,1)
    ids=i:min(i+31,size(pos,1));
    dist=pdist2(pos,pos(ids,:));
    keep=dist>0 & dist<=120; shell=ceil(dist(keep)/.2);
    dotProd=u*u(ids,:)'; angle=acosd(max(-1,min(1,dotProd(keep))));
    angularBin=min(floor(angle)+1,bins);
    histCount=histCount+accumarray([shell angularBin],1,[nShell bins]);
    sumAngle=sumAngle+accumarray(shell,angle,[nShell 1]);
    sumSquare=sumSquare+accumarray(shell,angle.^2,[nShell 1]);
    posAngle=angle>0; keys=shell(posAngle); ang=angle(posAngle);
    sumPositive=sumPositive+accumarray(keys,ang,[nShell 1]);
    positiveCount=positiveCount+accumarray(keys,1,[nShell 1]);
    sumLog=sumLog+accumarray(keys,log(ang),[nShell 1]);
end
counts=sum(histCount,2); [~,modeBin]=max(histCount,[],2);
observed_peak_deg=modeBin-.5; observed_peak_deg(counts==0)=NaN;
normal_peak_deg=sumAngle./counts;
normal_sigma_deg=sqrt(max(0,(sumSquare-counts.*normal_peak_deg.^2)./(counts-1)));
positiveMean=sumPositive./positiveCount; target=log(positiveMean)-sumLog./positiveCount;
gamma_shape=nan(nShell,1); gamma_scale=gamma_shape;
for j=find(isfinite(target) & target>0 & positiveCount>1)'
    theta=target(j);
    equation=@(a) log(a)-psi(a)-theta;
    gamma_shape(j)=fzero(equation,[1e-9,max(10,2/theta)]);
    gamma_scale(j)=positiveMean(j)/gamma_shape(j);
end
gamma_peak_deg=max(0,gamma_shape-1).*gamma_scale;
out=fullfile(here,'output'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,'angle_distributions.mat'),'edges_A','angle_edges_deg','histCount',...
    'observed_peak_deg','normal_peak_deg','normal_sigma_deg','gamma_shape','gamma_scale','gamma_peak_deg');
fprintf('Angle distributions: %d shells; shared baseline correlation length %.3f A.\n',nShell,radial.xi_A(1));
end
