% Upstream inputs: lattice sites (06), POP/core directions/vectors (08),
% and the unique adjacency list (10). Each cell contains five columns.
makeProfiles(fileparts(mfilename('fullpath')));

function makeProfiles(here)
root = fileparts(fileparts(here));
fieldDir = fullfile(root,'08_polar&domain&POP'); addpath(fieldDir);
m = load(fullfile(root,'06_final_coordinates','Model_output_02_Final_model.mat'),'Blat_final');
f = load(fullfile(fieldDir,'output','polar_field.mat'));
a = load(fullfile(root,'10_Polar_order_and_spatial_correlations','Domain_adjacency','output','adjacent_pairs.mat'),'pairs');
lattice = pmn.lattice(m.Blat_final);
centers = f.centers_A/f.cell_A;
magnitudes = vecnorm(f.vectors,2,2);
unit = f.vectors./magnitudes;
profiles = cell(size(a.pairs,1),1);
for k = 1:numel(profiles)
    ids = a.pairs(k,:);
    [qr,site,bin] = pmn.tube(lattice,centers(ids,:),f.cell_A);
    angle = acosd(max(-1,min(1,unit(site,:)*f.core_direction(ids,:)')));
    obs = [f.POP(site),magnitudes(site)*f.polarization_per_pixel,angle];
    means = nan(numel(qr),4);
    for j = 1:4
        means(:,j) = accumarray(bin,obs(:,j),[numel(qr) 1],@mean,NaN);
    end
    profiles{k} = [qr means];
end
columns = {'qr_A','POP','polar_magnitude_uC_cm2','angle_start_deg','angle_end_deg'};
out = fullfile(here,'output'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,'profiles.mat'),'profiles','columns');
fprintf('Five-column profiles: %d pairs.\n',numel(profiles));
end
