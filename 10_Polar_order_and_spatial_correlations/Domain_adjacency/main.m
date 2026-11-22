% Read the accepted labels from 08 and B-site lattice coordinates from 06.
identifyPairs(fileparts(mfilename('fullpath')));

function identifyPairs(here)
root = fileparts(fileparts(here));
fieldDir = fullfile(root,'08_polar&domain&POP');
addpath(fieldDir);
m = load(fullfile(root,'06_final_coordinates','Model_output_02_Final_model.mat'),'Blat_final');
a = load(fullfile(fieldDir,'input','reviewed_domains.mat'),'domain_id');
xyz = pmn.lattice(m.Blat_final); dims = max(xyz,[],1)+1;
grid = zeros(dims,'uint8');
grid(sub2ind(dims,xyz(:,1),xyz(:,2),xyz(:,3))) = a.domain_id;
contacts = cell(3,1);
for axis = 1:3
    dest = xyz; dest(:,axis) = dest(:,axis)+1;
    other = grid(sub2ind(dims,dest(:,1),dest(:,2),dest(:,3)));
    take = a.domain_id~=0 & other~=0 & a.domain_id~=other;
    contacts{axis} = sort([a.domain_id(take),other(take)],2);
end
pairs = double(unique(vertcat(contacts{:}),'rows'));
out = fullfile(here,'output'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,'adjacent_pairs.mat'),'pairs');
fprintf('Unique face-adjacent pairs: %d.\n',size(pairs,1));
end
