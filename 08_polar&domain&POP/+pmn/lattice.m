function xyz = lattice(Blat_final)
% One-based cubic-lattice coordinates in the row order of the B sites.
xyz = double(Blat_final)-min(double(Blat_final),[],1)+1;
end
