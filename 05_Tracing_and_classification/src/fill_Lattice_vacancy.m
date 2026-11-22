
function [lat_arr_f,pset_f] = fill_Lattice_vacancy(Lat_index,input_pset)
lat_arr_f = Lat_index;
pset_f    = input_pset;

%% Fill in empty place
    ub = max(lat_arr_f,[],1); lb = min(lat_arr_f,[],1);
    iter_count = 0;
    %%
    % fill_count = 0;
    fill_flag = 1;
    while fill_flag==1
        %---- set para
        iter_count = iter_count+1;
        bli_x = lat_arr_f(:,1);
        bli_y = lat_arr_f(:,2);
        bli_z = lat_arr_f(:,3);
        fill_count = 0;
        lat_arr_fill =[];
        pos_fill     =[];
        for ix=lb(1):ub(1)
        for iy=lb(2):ub(2)
        for iz=lb(3):ub(3)
            %---- is current pos already in BLI？
            flag = any(bli_x==ix & bli_y==iy & bli_z==iz);
            %---- if is a 'non-pos'
            if flag==0
                fx1 = find( bli_x<ix  & bli_y==iy & bli_z==iz );
                fx2 = find( bli_x>ix  & bli_y==iy & bli_z==iz );
                fy1 = find( bli_x==ix & bli_y<iy  & bli_z==iz );
                fy2 = find( bli_x==ix & bli_y>iy  & bli_z==iz );
                fz1 = find( bli_x==ix & bli_y==iy & bli_z<iz  );
                fz2 = find( bli_x==ix & bli_y==iy & bli_z>iz  );
                if (~isempty(fx1) + ~isempty(fx2) + ~isempty(fy1) +...  vacany is inside the volume
                    ~isempty(fy2) + ~isempty(fz1) + ~isempty(fz2))>4
                   fill_count = fill_count+1;
                   lat_arr_fill(fill_count,:) = [ix,iy,iz];         % fill Lattice indice
                   %--- fill pos set
                   c_fillpos = [];
                   if ~isempty(fx1) && ~isempty(fx2)
                       [cx1,ci1] = max(bli_x(fx1));
                       cpos1 = pset_f(fx1(ci1),:);
                       [cx2,ci2] = min(bli_x(fx2));
                       cpos2 = pset_f(fx2(ci2),:);
                       c_fillpos = cat(1,c_fillpos,(cpos1*(cx2-ix)+cpos2*(ix-cx1))/(cx2-cx1));
                   end
                   if ~isempty(fy1) && ~isempty(fy2)
                       [cy1,ci1] = max(bli_y(fy1));
                       cpos1 = pset_f(fy1(ci1),:);
                       [cy2,ci2] = min(bli_y(fy2));
                       cpos2 = pset_f(fy2(ci2),:);
                       c_fillpos = cat(1,c_fillpos,(cpos1*(cy2-iy)+cpos2*(iy-cy1))/(cy2-cy1));
                   end
                   if  ~isempty(fz1) && ~isempty(fz2)
                       [cz1,ci1] = max(bli_z(fz1));
                       cpos1 = pset_f(fz1(ci1),:);
                       [cz2,ci2] = min(bli_z(fz2));
                       cpos2 = pset_f(fz2(ci2),:);
                       c_fillpos = cat(1,c_fillpos,(cpos1*(cz2-iz)+cpos2*(iz-cz1))/(cz2-cz1));
                   end
                   pos_fill(fill_count,:) = mean(c_fillpos);
                end
            end
        end 
        end
        end
        fprintf('iter:%d, fill_count:%d\n',[iter_count fill_count])
        if isempty(lat_arr_fill)
            fill_flag=0;
        end
        % whos pos_fill
        lat_arr_f = cat(1,lat_arr_f,lat_arr_fill);
        pset_f    = cat(1,pset_f,pos_fill);
    end
    % BLI_output = lat_arr_f;
end
