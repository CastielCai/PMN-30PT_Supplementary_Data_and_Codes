function[pos_ind, valid_ind] = find_pos_ind_dim1(pos,full_pset)   %
%%
N = size(pos,1);
valid_ind = 1:N;
pos_ind = zeros(1,N);
posdim = size(pos,2);
for ii=1:N
    if posdim==2
        c_ind = find(pos(ii,1)==full_pset(:,1) & pos(ii,2)==full_pset(:,2));
    elseif posdim==3
        c_ind = find(pos(ii,1)==full_pset(:,1) & pos(ii,2)==full_pset(:,2) & pos(ii,3)==full_pset(:,3));
    end
    if isempty(c_ind)
        continue
    end
    pos_ind(ii) = c_ind;
end
if any(pos_ind==0)
    valid_ind(pos_ind==0)=[];
    pos_ind(pos_ind==0)=[];
end
end