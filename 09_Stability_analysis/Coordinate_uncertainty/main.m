% Read shared experimental coordinates and Pb cages from 06, and the
% unperturbed displacement field from 08. Values reported in pm or Angstrom.
coordinateTrials(fileparts(mfilename('fullpath')));

function coordinateTrials(here)
root=fileparts(fileparts(here));
fieldDir=fullfile(root,'08_polar&domain&POP'); addpath(fieldDir); addpath(here);
m=load(fullfile(root,'06_final_coordinates','Model_output_02_Final_model.mat'),'A_trace','polar_set_f','Blat_final');
f=load(fullfile(fieldDir,'output','polar_field.mat'),'pixel_A','vectors');
stream=load(fullfile(here,'input','trial_streams.mat'));
xyz=pmn.lattice(m.Blat_final);
S=pmn.neighbors(xyz,3,25/18); [~,support]=pmn.neighbors(xyz,2,Inf);
interior_index=find(support==19);
nB=size(xyz,1); nLead=size(m.A_trace,1); cages=double(m.polar_set_f(:,2:9));
C=sparse(repmat((1:nB)',8,1),cages(:),1/8,nB,nLead);
reference_A=f.vectors(interior_index,:)*f.pixel_A;
S=S(interior_index,:); SC=S*C;
nTrial=size(stream.initial,1); trial_rmsd_pm=zeros(nTrial,1);
for k=1:nTrial
    [ea,eb]=trial_errors(stream,k,nLead,nB);
    delta=S*eb-SC*ea;
    trial_rmsd_pm(k)=sqrt(mean(sum(delta.^2,2)));
    if k==1
        trial0_A=reference_A+delta/100;
        rms_difference_A=sqrt(mean((delta/100).^2,1));
        pearson_r=diag(corr(reference_A,trial0_A))';
    end
    if mod(k,100)==0, fprintf('Coordinate perturbations %d/%d\n',k,nTrial); end
end
mean_rmsd_pm=mean(trial_rmsd_pm);
percentile_pm=interp1(1:nTrial,sort(trial_rmsd_pm),1+(nTrial-1)*[.025 .975]);
out=fullfile(here,'output'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,'coordinate_uncertainty.mat'),'trial_rmsd_pm','mean_rmsd_pm','percentile_pm',...
    'interior_index','trial0_A','rms_difference_A','pearson_r');
fprintf('RMSD %.5f pm [%.5f, %.5f].\n',mean_rmsd_pm,percentile_pm);
end
