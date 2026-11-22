function J = regionfill3D(I,varargin) 
narginchk(2,3);
% Error checking
matlab.images.internal.errorIfgpuArray(I, varargin{:})
mask = logical(varargin{1});
%---- Find the outer mask boundary pixels
maskPerimeter = findBoundaryVoxels(mask);  
% Fill the region in I specified by mask using maskPerimeter as boundary
% condition for Laplace's equation
J = regionfillLaplace_3D(I,mask,maskPerimeter); 
end 

%% SubFunction
%% SubFunction
%% SubFunction 
%----------------------------------------
function  maskPerimeter = findBoundaryVoxels(mask)
    SE = strel("sphere",1);
    maskDilated  = imdilate(mask,SE); 
    maskPerimeter = maskDilated & ~mask; 
end

%--------------------------------------------------------------------------
% Fill the region in I specified by mask using maskPerimeter as boundary
% condition for Laplace's equation
function J = regionfillLaplace_3D(I,mask,maskPerimeter)
% If mask==ones(size(I))
if all(all(mask))
    % Warn the user that there are no boundary pixels
    % The returned image will be all black because the right side of the
    % equation is all zeros
    warning(message('images:regionfill:maskIsAllWhite'))
end

% Initialize output
J = I;

% Image size
[nRow,nCol,nLay] = size(I);

% Form the right side of the equation
rightSide = formRightSide_3D(I,maskPerimeter);

% Location of mask pixels
maskIdx = find(mask);

% Only keep values for pixels that are in the mask
rightSide = rightSide(maskIdx);

% Number the mask pixels in a grid matrix
grid = zeros(nRow,nCol,nLay);
grid(maskIdx) = 1:numel(maskIdx);
% Pad with zeros to avoid "index out of bounds" errors in the for loop
grid = padMatrix(grid);
gridIdx = find(grid);

%%------------------------------------------------------------------
% 构建三维的线性方程组并解之
% Form the connectivity matrix D=sparse(i,j,s)
% Connect each mask pixel to itself
i = (1:numel(maskIdx))';
j = (1:numel(maskIdx))';
% k = (1:numel(maskIdx))';
% The coefficient is the number of neighbors over which we average
numNeighbors = computeNumberOfNeighbors_3D(nRow,nCol,nLay);
s = numNeighbors(maskIdx);

% Now connect the N,E,S,W neighbors if they exist
%%% Old code %%% m = nRow+2; % number of rows in grid
%%% new: 3D should have 3 direction%%  
m = nRow+2; % number of rows in grid
n = (nRow+2)*(nCol+2);

for direction = [-1 m 1 -m n -n]
    % Possible neighbors in the current direction
    neighbors = grid(gridIdx+direction);
    % Connect mask points to neighbors with -1's
    i = [i; grid(gridIdx(neighbors~=0))]; %#ok<*AGROW>
    j = [j; nonzeros(neighbors)];
    % is there a equation for 'k'?  Ans:[no] 
    s = [s; -ones(nnz(neighbors),1)];
end

D = sparse(i,j,s);
% D = sparse(i,j,k,s); [wrong format]

% Solve the linear system
sol = D \ rightSide;    % ans is linear vecotr
%%-------------------------------------------------------------------
% Output
J(maskIdx) = sol;
end

%--------------------------------------------------------------------------
% 【【【Need to transfer to an 3D version】】】
% Form the right side of the linear equation to solve
% The right side is made of the known neighboring pixel values.
% In the general case the equation associated with a mask pixel (pixel to
% be filled) looks like:
%
%         4 x - x_n - x_s - x_e - x_w = 0
%
% 【3D version may look like】
% 【      6 x - x_n - x_s - x_e - x_w - x_t - x_b = 0         】

%
% where x is the unknown value of the current pixel and x_n, x_s, x_e, x_w
% are the values of the North, South, East, and West pixels.
% If the current pixel is on the edge of the mask, then at least one of
% these neighbor is known. This is what goes in "rightSide":
%
%         4 x - x_n - x_s - x_w = val_e
% 【3D version may look like】
% 【      6 x - x_n - x_s - x_e - x_w - x_b = val_t          】
%
% The right side of the equation is made of the sum of the known neighbor
% values.
function rightSide = formRightSide_3D(I, maskPerimeter)
% Image size
[nRow,nCol,nLay] = size(I);
% [nRow,nCol] = size(I);

% Get the values of the pixels on the mask perimeter
perimeterValues = zeros(nRow,nCol,nLay);
% perimeterValues = zeros(nRow,nCol);
perimeterValues(maskPerimeter) = I(maskPerimeter);

% Initialize
rightSide = zeros(nRow,nCol,nLay);
%%% IIIINTERIOR
% For pixels that are INTERIOR to the image, neighbors are N,S,W,E,T,B
rightSide(2:nRow-1,2:nCol-1,2:nLay-1) = ...
      perimeterValues(1:nRow-2,2:nCol-1,2:nLay-1) ... % N  dim1-1
    + perimeterValues(3:nRow  ,2:nCol-1,2:nLay-1) ... % S  dim1+1
    + perimeterValues(2:nRow-1,1:nCol-2,2:nLay-1) ... % W  dim2-1
    + perimeterValues(2:nRow-1,3:nCol  ,2:nLay-1) ... % E  dim2+1
    + perimeterValues(2:nRow-1,2:nCol-1,3:nLay  ) ... % T  dim3+1
    + perimeterValues(2:nRow-1,2:nCol-1,1:nLay-2);    % B  dim3-1

%%% FFFFcaes have 5 neighbors; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left face: N,S,E,T,B
rightSide(2:nRow-1,1,2:nLay-1) = ...
      perimeterValues(1:nRow-2,1,2:nLay-1) ... % N
    + perimeterValues(3:nRow  ,1,2:nLay-1) ... % S
    + perimeterValues(2:nRow-1,2,2:nLay-1) ... % E
    + perimeterValues(2:nRow-1,2,3:nLay  ) ... % T
    + perimeterValues(2:nRow-1,1,1:nLay-2);    % B      

% Right face: N,S,W,T,B
rightSide(2:nRow-1,nCol,2:nLay-1) = ... 
      perimeterValues(1:nRow-2,nCol  ,2:nLay-1) ... % N
    + perimeterValues(3:nRow  ,nCol  ,2:nLay-1) ... % S
    + perimeterValues(2:nRow-1,nCol-1,2:nLay-1) ... % W
    + perimeterValues(2:nRow-1,nCol  ,3:nLay  ) ... % T
    + perimeterValues(2:nRow-1,nCol  ,1:nLay-2);    % B     

% North face: S,W,E,T,B
rightSide(1,2:nCol-1,2:nLay-1) = ...
      perimeterValues(2,2:nCol-1,2:nLay-1) ... % S
    + perimeterValues(1,1:nCol-2,2:nLay-1) ... % W
    + perimeterValues(1,3:nCol  ,2:nLay-1) ... % E
    + perimeterValues(1,2:nCol-1,3:nLay  ) ... % T
    + perimeterValues(1,2:nCol-1,1:nLay-2) ;   % B

% South face: N,W,E,T,B
rightSide(nRow,2:nCol-1,2:nLay-1) = ...
      perimeterValues(nRow-1,2:nCol-1,2:nLay-1) ... % N
    + perimeterValues(nRow  ,1:nCol-2,2:nLay-1) ... % W
    + perimeterValues(nRow  ,3:nCol  ,2:nLay-1) ... % E
    + perimeterValues(nRow  ,2:nCol-1,3:nLay  ) ... % T
    + perimeterValues(nRow  ,2:nCol-1,1:nLay-2) ;   % B


% Z-top face: N,S,W,E,B
rightSide(2:nRow-1,2:nCol-1,nLay) = ...
      perimeterValues(1:nRow-2,2:nCol-1,nLay)   ... % N
    + perimeterValues(3:nRow  ,2:nCol-1,nLay)   ... % S
    + perimeterValues(2:nRow-1,1:nCol-2,nLay)   ... % W
    + perimeterValues(2:nRow-1,3:nCol  ,nLay)   ... % E
    + perimeterValues(2:nRow-1,2:nCol-1,nLay-1);    % B 

% Z-bottom face: N,S,W,E,T
rightSide(2:nRow-1,2:nCol-1,1) = ...
      perimeterValues(1:nRow-2,2:nCol-1,1)   ... % N
    + perimeterValues(3:nRow  ,2:nCol-1,1)   ... % S
    + perimeterValues(2:nRow-1,1:nCol-2,1)   ... % W
    + perimeterValues(2:nRow-1,3:nCol  ,1)   ... % E
    + perimeterValues(2:nRow-1,2:nCol-1,2);      % T

%%%%% EEEEdges have 4 neighbors;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% North&Top    edge: S,W,E,B
rightSide(1,2:nCol-1,nLay) = ... 
    + perimeterValues(2,2:nCol-1,nLay  )   ... % S
    + perimeterValues(1,1:nCol-2,nLay  )   ... % W
    + perimeterValues(1,3:nCol  ,nLay  )   ... % E
    + perimeterValues(1,2:nCol-1,nLay-1);      % B

% North&Bottom     edge: S,W,E,T
rightSide(1,2:nCol-1,1) = ... 
    + perimeterValues(2,2:nCol-1,1)   ... % S
    + perimeterValues(1,1:nCol-2,1)   ... % W
    + perimeterValues(1,3:nCol  ,1)   ... % E
    + perimeterValues(1,2:nCol-1,2);      % T

% North&Left     edge: S,E,T,B
rightSide(1,1,2:nLay-1) = ... 
    + perimeterValues(2,1,2:nLay-1)   ... % S
    + perimeterValues(1,2,2:nLay-1)   ... % W
    + perimeterValues(1,1,3:nLay)   ... % T
    + perimeterValues(1,1,1:nLay-2);      % B

% North&Right    edge: S,W,T,B
rightSide(1,nCol,2:nLay-1) = ... 
    + perimeterValues(2,nCol  ,2:nLay-1 )   ... % S
    + perimeterValues(1,nCol-1,2:nLay-1 )   ... % W
    + perimeterValues(1,nCol  ,3:nLay   )   ... % T
    + perimeterValues(1,nCol  ,2:nLay-1 );      % B

% South&Top  edge: N,W,E,B
rightSide(nRow,2:nCol-1,nLay) = ... 
    + perimeterValues(nRow-1,2:nCol-1,nLay  )   ... % N
    + perimeterValues(nRow  ,1:nCol-2,nLay  )   ... % W
    + perimeterValues(nRow  ,3:nCol  ,nLay  )   ... % E
    + perimeterValues(nRow  ,2:nCol-1,nLay-1);      % B

% South&Bottom   edge: N,W,E,T
rightSide(nRow,2:nCol-1,1) = ... 
    + perimeterValues(nRow-1,2:nCol-1,1)   ... % N
    + perimeterValues(nRow  ,1:nCol-2,1)   ... % W
    + perimeterValues(nRow  ,3:nCol  ,1)   ... % E
    + perimeterValues(nRow  ,2:nCol-1,2);      % T

% South&Left   edge: N,E,T,B
rightSide(nRow,1,2:nLay-1) = ... 
    + perimeterValues(nRow-1,1,2:nLay-1)   ... % N
    + perimeterValues(nRow  ,2,2:nLay-1)  ... % E
    + perimeterValues(nRow  ,1,3:nLay  )  ... % T
    + perimeterValues(nRow  ,1,1:nLay-2);      % B

% South&Right  edge: N,W,T,B
rightSide(nRow,nCol,2:nLay-1) = ... 
    + perimeterValues(nRow-1 ,nCol  ,2:nLay-1)   ... % N
    + perimeterValues(nRow   ,nCol-1,2:nLay-1)  ... % W
    + perimeterValues(nRow   ,nCol  ,3:nLay  )  ... % T
    + perimeterValues(nRow   ,nCol  ,1:nLay-2);      % B

% Left&Top  edge: N,S,E,B
rightSide(2:nRow-1,1,nLay) = ... 
    + perimeterValues(1:nRow-2,1,nLay )   ... % N
    + perimeterValues(3:nRow  ,1,nLay )   ... % S
    + perimeterValues(2:nRow-1,2,nLay )   ... % E
    + perimeterValues(2:nRow-1,1,nLay-1);     % B

% Right&Top edge: N,S,W,T
rightSide(2:nRow-1,nCol,nLay) = ... 
    + perimeterValues(1:nRow-2,nCol  ,nLay )   ... % N
    + perimeterValues(3:nRow  ,nCol  ,nLay )   ... % S
    + perimeterValues(2:nRow-1,nCol-1,nLay )   ... % W
    + perimeterValues(2:nRow-1,nCol  ,nLay-1);     % B

% Left&Bottom   edge: N,S,E,T
rightSide(2:nRow-1,1,1) = ... 
    + perimeterValues(1:nRow-2,1,1)   ... % N
    + perimeterValues(3:nRow  ,1,1)   ... % S
    + perimeterValues(2:nRow-1,2,1)   ... % E
    + perimeterValues(2:nRow-1,1,2);      % T

% Right&Bottom  edge: N,S,W,B
rightSide(2:nRow-1,nCol,1) = ... 
    + perimeterValues(1:nRow-2,nCol  ,1)   ... % N
    + perimeterValues(3:nRow  ,nCol  ,1)   ... % S
    + perimeterValues(2:nRow-1,nCol-1,1)   ... % W
    + perimeterValues(2:nRow-1,nCol  ,2);      % T

% Corners have 3 neighbors
rightSide(1,1,1)          = perimeterValues(1,2) + perimeterValues(2,1); % E+S+B
rightSide(1,nCol,1)       = perimeterValues(1,nCol-1) + perimeterValues(2,nCol); % W+S+B
rightSide(nRow,1,1)       = perimeterValues(nRow-1,1) + perimeterValues(nRow,2); % N+E+B
rightSide(nRow,nCol,1)    = perimeterValues(nRow-1,nCol) + perimeterValues(nRow,nCol-1); % N+W+B
rightSide(1,1,nLay)       = perimeterValues(1,2) + perimeterValues(2,1); % E+S+T
rightSide(1,nCol,nLay)    = perimeterValues(1,nCol-1) + perimeterValues(2,nCol); % W+S+T
rightSide(nRow,1,nLay)    = perimeterValues(nRow-1,1) + perimeterValues(nRow,2); % N+E+T
rightSide(nRow,nCol,nLay) = perimeterValues(nRow-1,nCol) + perimeterValues(nRow,nCol-1); % N+W+T
end

%--------------------------------------------------------------------------

% Return a matrix of size (nRow,nCol) that contains the number of neighbors
% for each pixel. Interior pixels have 4 neighbors. Border pixels have 3.
% Corner pixels have 2.
% This is used when forming the linear equations to be solved. 
function numNeighbors = computeNumberOfNeighbors_3D(nRow,nCol,nLay)
% Initialize
numNeighbors = zeros(nRow,nCol,nLay);
% Interior voxels have 6 neighbors
numNeighbors(2:nRow-1,2:nCol-1,2:nLay-1) = 6;
% Face voxels have 5 neighbors
numNeighbors([1 nRow],2:nCol-1,2:nLay-1) = 5;
numNeighbors(2:nRow-1,[1 nCol],2:nLay-1) = 5;
numNeighbors(2:nRow-1,2:nCol-1,[1 nLay]) = 5;
% Edge voxels have 4 neighbors
numNeighbors(2:nRow-1,[1 nCol],[1 nLay]) = 4;
numNeighbors([1 nRow],2:nCol-1,[1 nLay]) = 4;
numNeighbors([1 nRow],[1 nCol],2:nLay-1) = 4; 
% Corner voxels have 3 neighbors
numNeighbors([1 nRow],[1 nCol],[1 nLay]) = 3;
end

%--------------------------------------------------------------------------
% Pad the input matrix
% When the mask touches the border of the image, we must pad the grid to
% avoid index out of bounds errors.
function gridPadded = padMatrix(grid)
[nRow,nCol,nLay] = size(grid);
gridPadded = zeros(nRow+2,nCol+2,nLay+2);
gridPadded(2:nRow+1,2:nCol+1,2:nLay+1) = grid;
gridPadded = cast(gridPadded,'like',grid);
end

%%



