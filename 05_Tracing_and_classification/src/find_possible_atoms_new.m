function [atom_pos, close_pos, stats, statsI] = find_possible_atoms_new(DataMatrix, Dmin, MaxNumberAtoms, BoxSize0,BoxSize1,BoxSize2,Th)

[box1coordinates, BoxCenter1, BoxRadius1, sphere1] = create_box(BoxSize1);
[box2coordinates, BoxCenter2, BoxRadius2, sphere2] = create_box(BoxSize2);

atom_pos=zeros(3,MaxNumberAtoms);   %store position of peaks
close_pos=[];                       %store position of peaks that are too close
nclose=0;                           %counter for too close atoms
i_atom=0;                           
i_maxima=0;                         %index for current atom
band_ind=0;
closeN=0;
[Xsize, Ysize, Zsize] = size(DataMatrix);

CurrData=DataMatrix;                %Matrix used for Gaussian substraction
M=find_maxima3D(CurrData,sphere1);  %find all maxima
Peaks = DataMatrix; Peaks(M==0) = 0;
M(Peaks < Th) = 0;
fprintf('Number of maxima found: %g\n', sum(M(:)));
Nmax = sum(M(:));

% find coordinates of peaks
ind_find=find(M);
if isempty(ind_find)
    error('No local maxima was found!!')
end   
[xx,yy,zz] = ind2sub(size(M),ind_find);        

% find peaks and sort them
[peak_sort, ind_sort] = get_peak_values3D(CurrData,xx,yy,zz,ind_find);
noatomM = zeros(size(DataMatrix));
% 3D mask for delete traced sphere
r = round(Dmin);
[x,y,z] = meshgrid(-r:r,-r:r,-r:r);
sph = (x.^2+y.^2+z.^2)<Dmin^2;
sph_ind = find(sph);
[sph_posx, sph_posy, sph_posz]  = ind2sub(size(sph),sph_ind);
sph_pos = [sph_posx'; sph_posy';sph_posz']-[r+1;r+1;r+1];

%continue this loop basicaly forever
%"inf" in a for loop is ugly, therefore 1e8
for iteration = 1:Nmax
    i_maxima=i_maxima+1;
    %take highest peak
    hind=ind_sort(iteration); hpeak=peak_sort(iteration);
    if ~isempty(find(band_ind==hind,1))
        continue
    end
    i_atom=i_atom+1;
    fprintf('Atom number: %04i\n',i_atom)
    %peak location
    x=xx(hind);y=yy(hind);z=zz(hind);
    
    %Get Data ROI
    [x y z];
    DataBox1 = get_data_box(CurrData,x,y,z,BoxRadius1);
    
    %start values and bounds for fitting       (x1   x2  x3,4,5:center  x6: x7: x8:      x9:psi  x10:theta  x11:phi  )
    ymax_init       = DataBox1(BoxCenter1,BoxCenter1,BoxCenter1);
    fit_param_init  = [0,  ymax_init, 0,  0,  0,  0.2,  0.2,  0.8,  0.0,  0.0,  0.0];
    fixed           = [0,  0   ,  0,  0,  0,    0,    0,    0,    0,    0,    0];
    lb              = [0,  0    , -1*BoxRadius1, -1*BoxRadius1, -1*BoxRadius1,    0,    0,    0, -pi, 0, -pi];
    ub              = [inf,inf      ,  BoxRadius1,  BoxRadius1,  BoxRadius1,  inf,  inf,  inf, pi, pi,  pi];
    
    %perform the fit
    [fit_resultI,resnorm,residual] = fit_gauss3D_PD(fit_param_init, box1coordinates, DataBox1, fixed, lb, ub);
 
    %calc position of new atom and FWHM
    new_x = x+fit_resultI(3);
    new_y = y+fit_resultI(4);
    new_z = z+fit_resultI(5);
    
    %integer new position
    rx=round(new_x);
    ry=round(new_y);
    rz=round(new_z);
    
    %sub pixel offset to integer position
    dx=new_x-rx;
    dy=new_y-ry;
    dz=new_z-rz;
    
    %get minimum distance
    if i_atom>1
        d0 = get_minimum_distance(new_x,new_y,new_z,atom_pos(:,1:i_atom-1));
    else
        d0 = inf;
    end
    
    valid_atom=1;
    
    %local area used to calculate maxima
    roisize=25; 
    roihalf=(roisize-1)/2;
    xhi=new_x+roihalf;
    yhi=new_y+roihalf;
    zhi=new_z+roihalf;

    xlo=new_x-roihalf;
    ylo=new_y-roihalf;
    zlo=new_z-roihalf;
    
    if xhi > Xsize || yhi > Ysize || zhi > Zsize || xlo < 1 || ylo < 1 || zlo < 1
        fprintf('Atom number: %04i \t outside volume!\n',i_atom)
        valid_atom=0;
    end
    
    %check minimum distance constraint
    if d0<Dmin
        fprintf('Atom number: %04i \t Too close!\n',i_atom)
        closeN = closeN+1;
        valid_atom=0;
    end
    
    %check, if new fit position is out of bounce
    %should not happen, if fit bounds are set correctly
    if new_x>Xsize || new_y>Ysize || new_z>Zsize || new_x<1 || new_y<1 || new_z<1
        fprintf('Atom number: %04i: \t Out of Border!\n',i_atom)
        valid_atom=0;
    end
    
    if valid_atom==0     %atom is to close
        
        nclose=nclose+1;
        close_pos(:,nclose)=[new_x new_y new_z];
        
        M(x,y,z) = 0;         %remove this peak from the list of potential atoms
        i_atom = i_atom-1;  %decrease atom counter, since it is not an atom
        noatomM(x,y,z) = 1;
    
    else                  %yes, it is a valid peak
        
        %substract fit from data
        
        fit_resultI(3:5)=[dx dy dz];                           %use subpixel offsets
        FitBox2 = calc_gauss3D_PD(fit_resultI, box2coordinates);  %calculate gaussian
        FitBox2 = FitBox2-fit_resultI(1);                      %substract background
        
        DataBox2 = get_data_box(CurrData,rx,ry,rz,BoxRadius2); %get data
        DiffBox2 = DataBox2-FitBox2;                           %substact gaussian fit
        DiffBox2 = max(0,DiffBox2);                            %remove negative values
        
        %update data for next iteration
        CurrData(rx-BoxRadius2:rx+BoxRadius2,ry-BoxRadius2:ry+BoxRadius2,rz-BoxRadius2:rz+BoxRadius2)=DiffBox2;
        
        %save statistics
        statsI.resnorm(i_atom) = resnorm;
        statsI.residual{i_atom} = residual(:);
        statsI.fit{i_atom}=fit_resultI;
        statsI.peak_bg(i_atom)=fit_resultI(1);
        statsI.peak_hight(i_atom)=fit_resultI(2);
        statsI.peak_FWHM(:,i_atom)=2.3548./(sqrt(2)*fit_resultI(6:8));
        statsI.localmax_pos{i_atom} = [x y z];

        
        %save atom position
        atom_pos(1,i_atom) = new_x;
        atom_pos(2,i_atom) = new_y;
        atom_pos(3,i_atom) = new_z;
        
        %% 
        fit_result  = stats_cal_PD(DataMatrix,atom_pos(:,i_atom),atom_pos(:,i_atom),BoxSize1);
          
        stats.fit{i_atom}=fit_result;
        stats.peak_bg(i_atom)=fit_result(1);
        stats.peak_hight(i_atom)=fit_result(2);
        stats.peak_FWHM(:,i_atom)=2.3548./(sqrt(2)*fit_result(6:8)); 
        stats.localmax_pos{i_atom} = [x y z];
          
        %find maxima again, but restrict to the region next to the
        %currently substracted atom
        
        %% do not trace current atom again
        band_pos = sph_pos+[rx;ry;rz];
%         band_x = rx+sph_pos(1,:);
%         band_y = ry+sph_pos(2,:);
%         band_z = rz+sph_pos(3,:);
        new_band_ind = sub2ind(size(M),band_pos(1,:),band_pos(2,:),band_pos(3,:));
        band_ind = cat(2,band_ind,new_band_ind);
        noatomM(new_band_ind)=1;   %mask for all the found atoms
    end
    
    if i_atom==MaxNumberAtoms || i_maxima==Nmax
        break
    end
    
end
closeN
atom_pos(:,i_atom+1:end)=[];


