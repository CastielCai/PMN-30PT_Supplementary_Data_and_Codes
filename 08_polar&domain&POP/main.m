% Shared experimental model: 06_final_coordinates. Source labels are the
% accepted domain and core masks in the same B-site row order as that model.
computeField(fileparts(mfilename('fullpath')));

function computeField(here)
model = load(fullfile(fileparts(here),'06_final_coordinates','Model_output_02_Final_model.mat'));
mask = load(fullfile(here,'input','reviewed_domains.mat'));
assert(all(model.B_prop(:,1)~=0),'The released model must have complete B sites.');
lattice = pmn.lattice(model.Blat_final);
pixel_A = 0.3079;
cell_A = 4;
polarization_per_pixel = 79.4382; % microC/cm^2 per displacement pixel
cages = double(model.polar_set_f(:,2:9));
lead = double(model.A_trace);
cageMean = zeros(size(cages,1),3);
for j = 1:8
    cageMean = cageMean + lead(cages(:,j),:)/8;
end
raw_displacement_pixel = (double(model.B_tr_rep)-cageMean)/3;
smooth = pmn.neighbors(lattice,3,25/18);
vectors = smooth*raw_displacement_pixel;
localAverage = pmn.neighbors(lattice,6,Inf);
POP = sum(vectors.*(localAverage*vectors),2);
nDomains = double(max(mask.domain_id));
centers_A = zeros(nDomains,3);
core_direction = zeros(nDomains,3);
for d = 1:nDomains
    selected = mask.core_id==d;
    assert(any(selected),'Missing core %d.',d);
    centers_A(d,:) = cell_A*mean(lattice(selected,:),1);
    meanVector = mean(vectors(selected,:),1);
    core_direction(d,:) = meanVector/norm(meanVector);
end
out = fullfile(here,'output');
if ~isfolder(out), mkdir(out); end
save(fullfile(out,'polar_field.mat'),'raw_displacement_pixel','vectors','POP',...
    'centers_A','core_direction','cell_A','pixel_A','polarization_per_pixel');
fprintf('Polar field: %d cells, %d cores.\n',numel(POP),nDomains);
end
