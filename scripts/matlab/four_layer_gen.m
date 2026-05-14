clear all
close all
clc

atreyu_folder_name = 'ERCEM321_Atreyu'; 

atreyu_path = fullfile(pwd, '..', atreyu_folder_name);
atreyu_examples_path = fullfile(atreyu_path, 'Examples', 'matlab_shared_lib_ThreeLayers');

addpath(fullfile(atreyu_examples_path, 'aux'));
addpath(fullfile(atreyu_examples_path, 'mesh'));

% Cross-Platform Library Loading
if ispc % Windows OS
    library_name = 'matlab_shared_lib_ThreeLayers';
    library_filename = 'matlab_shared_lib_ThreeLayers.dll';
    libpath = fullfile(atreyu_path, 'build', 'bin', 'Release', library_filename);
else % Linux/Mac OS
    library_name = 'libmatlab_shared_lib_ThreeLayers';
    library_filename = 'libmatlab_shared_lib_ThreeLayers.so';
    libpath = fullfile(atreyu_path, 'build', 'lib', library_filename);
end

header_name = regexprep(library_filename, '\.(so|dll)$', '.h');
header_name = regexprep(header_name, '^lib', '');
libheaderpath = fullfile(atreyu_examples_path, header_name);

%% Dataset
dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
mat_dir = fullfile(dataset_dir, 'mat_files');
msh_dir = fullfile(dataset_dir, 'msh_files');

%% Load matlab_shared_lib_example libraries
if libisloaded(library_name), unloadlibrary(library_name); end
[notfound,warnings]=loadlibrary(libpath, libheaderpath);
libfunctions(library_name, "-full")

%% Load meshes
brain_msh = fullfile(msh_dir, 'Surface_cortical_surface.msh');
brain_gmsh = fullfile(msh_dir, 'cortical_surface.gmsh');
brain_skull_msh = fullfile(msh_dir, 'Surface_brain_skull.msh');
brain_skull_gmsh = fullfile(msh_dir, 'Surface_brain_skull.gmsh');
skull_skin_msh = fullfile(msh_dir, 'Surface_skull_skin.msh');
skull_skin_gmsh = fullfile(msh_dir, 'Surface_skull_skin.gmsh');
skin_air_msh = fullfile(msh_dir, 'Surface_skin_air.msh');
skin_air_gmsh = fullfile(msh_dir, 'Surface_skin_air.gmsh');

struct_path = fullfile(mat_dir, 'Surfaces_structure.mat');
load(struct_path);

%% Create dipoles
dr = 0.0015;

norm_magnitude = vecnorm(surf_struct.cortical_surface.vertex_normals, 2, 1); 
momDipole_list = surf_struct.cortical_surface.vertex_normals ./ norm_magnitude;

posDipole_list = surf_struct.cortical_surface.vertices - (dr .* momDipole_list);

[~, dist_nearest] = knnsearch(surf_struct.cortical_surface.vertices', posDipole_list');

tollerance = dr * 0.05; 
bad_idx = find(dist_nearest < (dr - tollerance));
good_idx = setdiff(1:size(posDipole_list,2), bad_idx);

%% Create spaces

sigma = [1.0, 1.79, 0.01, 1.0];

[points1, cells1, normals1, dim1, nNode1] = readMshFile_Samuele(brain_msh);
nPoints1 = size(points1,1);
nCells1 = size(cells1,1);

[points2, cells2, normals2, dim2, nNode2] = readMshFile_Samuele(brain_skull_msh);
nPoints2 = size(points2,1);
nCells2 = size(cells2,1);

[points3, cells3, normals3, dim3, nNode3] = readMshFile_Samuele(skull_skin_msh);
nPoints3 = size(points3,1);
nCells3 = size(cells3,1);

[points4, cells4, normals4, dim4, nNode4] = readMshFile_Samuele(skin_air_msh);
nPoints4 = size(points4,1);
nCells4 = size(cells4,1);

if(exist(brain_gmsh,'file')~=2)
    gmsh_mesh3d_write(brain_gmsh, dim1, nPoints1, points1, nNode1, nCells1, cells1);
end

if(exist(brain_skull_gmsh,'file')~=2)
    gmsh_mesh3d_write(brain_skull_gmsh, dim2, nPoints2, points2, nNode2, nCells2, cells2);
end

if(exist(skull_skin_gmsh,'file')~=2)
    gmsh_mesh3d_write(skull_skin_gmsh, dim3, nPoints3, points3, nNode3, nCells3, cells3);
end

if(exist(skin_air_gmsh,'file')~=2)
    gmsh_mesh3d_write(skin_air_gmsh, dim4, nPoints4, points4, nNode4, nCells4, cells4);
end

space1 = calllib(library_name, 'createSpace', brain_gmsh);
space2 = calllib(library_name, 'createSpace', brain_skull_gmsh);
space3 = calllib(library_name, 'createSpace', skull_skin_gmsh);
space4 = calllib(library_name, 'createSpace', skin_air_gmsh);

%% Data from spaces

funcSpacePatch1 = calllib(library_name, 'createFunctionalSpacePatch', space1);
funcSpacePatch2 = calllib(library_name, 'createFunctionalSpacePatch', space2);
funcSpacePatch3 = calllib(library_name, 'createFunctionalSpacePatch', space3);
funcSpacePatch4 = calllib(library_name, 'createFunctionalSpacePatch', space4);


funcSpacePyramid1 = calllib(library_name, 'createFunctionalSpacePyramid', space1);
funcSpacePyramid2 = calllib(library_name, 'createFunctionalSpacePyramid', space2);
funcSpacePyramid3 = calllib(library_name, 'createFunctionalSpacePyramid', space3);
funcSpacePyramid4 = calllib(library_name, 'createFunctionalSpacePyramid', space4);


funcSpaceGradPyramid1 = calllib(library_name, 'createFunctionalSpaceGradPyramid', space1);
funcSpaceGradPyramid2 = calllib(library_name, 'createFunctionalSpaceGradPyramid', space2);
funcSpaceGradPyramid3 = calllib(library_name, 'createFunctionalSpaceGradPyramid', space3);
funcSpaceGradPyramid4 = calllib(library_name, 'createFunctionalSpaceGradPyramid', space4);

nc1 = calllib(library_name, 'getNbPatchFunctions', funcSpacePatch1);
nc2 = calllib(library_name, 'getNbPatchFunctions', funcSpacePatch2);
nc3 = calllib(library_name, 'getNbPatchFunctions', funcSpacePatch3);
nc4 = calllib(library_name, 'getNbPatchFunctions', funcSpacePatch4);

nv1 = calllib(library_name, 'getNbPyramidFunctions', funcSpacePyramid1);
nv2 = calllib(library_name, 'getNbPyramidFunctions', funcSpacePyramid2);
nv3 = calllib(library_name, 'getNbPyramidFunctions', funcSpacePyramid3);
nv4 = calllib(library_name, 'getNbPyramidFunctions', funcSpacePyramid4);


nctot = nc1 + nc2 + nc3 + nc4;
nvtot = nv1 + nv2 + nv3 + nv4;

%% Evaluate h parameter

x1 = zeros(nv1,1); y1 = zeros(nv1,1); z1 = zeros(nv1,1);
[x1, y1, z1] = calllib(library_name, 'getVertices',...
    x1, y1, z1, space1);
points1 = [x1, y1, z1];

x2 = zeros(nv2,1); y2 = zeros(nv2,1); z2 = zeros(nv2,1);
[x2, y2, z2] = calllib(library_name, 'getVertices',...
    x2, y2, z2, space2);
points2 = [x2, y2, z2];

x3 = zeros(nv3,1); y3 = zeros(nv3,1); z3 = zeros(nv3,1);
[x3, y3, z3] = calllib(library_name, 'getVertices',...
    x3, y3, z3, space3);
points3 = [x3, y3, z3];

x4 = zeros(nv4,1); y4 = zeros(nv4,1); z4 = zeros(nv4,1);
[x4, y4, z4] = calllib(library_name, 'getVertices',...
    x4, y4, z4, space4);
points4 = [x4, y4, z4];

v1_a = zeros(nc1,1);
v1_b = zeros(nc1,1);
v1_c = zeros(nc1,1);
[v1_a, v1_b, v1_c] = calllib(library_name, 'getCells',...
    v1_a, v1_b, v1_c, space1);
cells1 = [v1_a, v1_b, v1_c]+1;

v2_a = zeros(nc2,1);
v2_b = zeros(nc2,1);
v2_c = zeros(nc2,1);
[v2_a, v2_b, v2_c] = calllib(library_name, 'getCells',...
    v2_a, v2_b, v2_c, space2);
cells2 = [v2_a, v2_b, v2_c]+1;

v3_a = zeros(nc3,1);
v3_b = zeros(nc3,1);
v3_c = zeros(nc3,1);
[v3_a, v3_b, v3_c] = calllib(library_name, 'getCells',...
    v3_a, v3_b, v3_c, space3);
cells3 = [v3_a, v3_b, v3_c]+1;

v4_a = zeros(nc4,1);
v4_b = zeros(nc4,1);
v4_c = zeros(nc4,1);
[v4_a, v4_b, v4_c] = calllib(library_name, 'getCells',...
    v4_a, v4_b, v4_c, space4);
cells4 = [v4_a, v4_b, v4_c]+1;

%h1
edges1 = computeEdges(cells1);
lengthE1 = getLengthEdges(edges1,points1);
lengthSum1 = sum(lengthE1);
h1 = lengthSum1/length(lengthE1);

% h2
edges2 = computeEdges(cells2);
lengthE2 = getLengthEdges(edges2,points2);
lengthSum2 = sum(lengthE2);
h2 = lengthSum2/length(lengthE2);

% h3
edges3= computeEdges(cells3);
lengthE3 = getLengthEdges(edges3,points3);
lengthSum3 = sum(lengthE3);
h3 = lengthSum3/length(lengthE3);

% h4
edges4 = computeEdges(cells4);
lengthE4 = getLengthEdges(edges4,points4);
lengthSum4 = sum(lengthE4);
h4 = lengthSum4/length(lengthE4);

hlist = (h1+h2+h3+h4)/4;

N_list = nvtot+nctot; %number of unknowns

%% S matrices (Safe Allocation)
Z_R = zeros(nc1, nc1); Z_I = zeros(nc1, nc1); [Z_R, ~] = calllib(library_name, 'computeSPatch', Z_R, Z_I, funcSpacePatch1, funcSpacePatch1); S11 = Z_R;
Z_R = zeros(nc1, nc2); Z_I = zeros(nc1, nc2); [Z_R, ~] = calllib(library_name, 'computeSPatch', Z_R, Z_I, funcSpacePatch1, funcSpacePatch2); S12 = Z_R; S21 = S12';
Z_R = zeros(nc2, nc2); Z_I = zeros(nc2, nc2); [Z_R, ~] = calllib(library_name, 'computeSPatch', Z_R, Z_I, funcSpacePatch2, funcSpacePatch2); S22 = Z_R;
Z_R = zeros(nc3, nc3); Z_I = zeros(nc3, nc3); [Z_R, ~] = calllib(library_name, 'computeSPatch', Z_R, Z_I, funcSpacePatch3, funcSpacePatch3); S33 = Z_R;
Z_R = zeros(nc2, nc3); Z_I = zeros(nc2, nc3); [Z_R, ~] = calllib(library_name, 'computeSPatch', Z_R, Z_I, funcSpacePatch2, funcSpacePatch3); S23 = Z_R; S32 = S23';

%% N matrices 
Z_R = zeros(nv1, nv1); Z_I = zeros(nv1, nv1); [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid1, funcSpacePyramid1); N11 = Z_R;
Z_R = zeros(nv1, nv2); Z_I = zeros(nv1, nv2); [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid1, funcSpacePyramid2); N12 = Z_R; N21 = N12';
Z_R = zeros(nv2, nv2); Z_I = zeros(nv2, nv2); [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid2, funcSpacePyramid2); N22 = Z_R;
Z_R = zeros(nv3, nv3); Z_I = zeros(nv3, nv3); [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid3, funcSpacePyramid3); N33 = Z_R;
Z_R = zeros(nv3, nv2); Z_I = zeros(nv3, nv2); [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid3, funcSpacePyramid2); N32 = Z_R; N23 = N32';
Z_R = zeros(nv4, nv4); Z_I = zeros(nv4, nv4); [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid4, funcSpacePyramid4); N44 = Z_R;
Z_R = zeros(nv4, nv3); Z_I = zeros(nv4, nv3); [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid4, funcSpacePyramid3); N43 = Z_R; N34 = N43';

%% D matrices 
Z_R = zeros(nc1, nv1); Z_I = zeros(nc1, nv1); [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch1, funcSpacePyramid1); D11 = Z_R; D11_compl=D11';
Z_R = zeros(nc1, nv2); Z_I = zeros(nc1, nv2); [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch1, funcSpacePyramid2); D12 = Z_R; D21_compl=D12';
Z_R = zeros(nc2, nv1); Z_I = zeros(nc2, nv1); [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch2, funcSpacePyramid1); D21 = Z_R; D12_compl=D21';
Z_R = zeros(nc2, nv2); Z_I = zeros(nc2, nv2); [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch2, funcSpacePyramid2); D22 = Z_R; D22_compl=D22';
Z_R = zeros(nc2, nv3); Z_I = zeros(nc2, nv3); [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch2, funcSpacePyramid3); D23 = Z_R; D32_compl=D23';
Z_R = zeros(nc3, nv3); Z_I = zeros(nc3, nv3); [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch3, funcSpacePyramid3); D33 = Z_R; D33_compl=D33';
Z_R = zeros(nc3, nv2); Z_I = zeros(nc3, nv2); [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch3, funcSpacePyramid2); D32 = Z_R; D23_compl=D32';
Z_R = zeros(nc3, nv4); Z_I = zeros(nc3, nv4); [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch3, funcSpacePyramid4); D34 = Z_R; D43_compl=D34';

%% Z and forward

values = zeros(1, nv4*nv4); rowIndices = values; colIndices = values;
[values, rowIndices, colIndices] = calllib(library_name,...
    'computeGramMatrixPyramid', values, rowIndices, colIndices, funcSpacePyramid4);
occupiedIndices = find(rowIndices);
GPyramid4 = sparse(rowIndices(1:occupiedIndices(end)), colIndices(1:occupiedIndices(end)), values(1:occupiedIndices(end)), nv4, nv4);

sigma_12 = sigma(1) + sigma(2); sigma_12inv = 1/sigma(1) + 1/sigma(2);
sigma_23 = sigma(2) + sigma(3); sigma_23inv = 1/sigma(2) + 1/sigma(3);
sigma_34 = sigma(3) + sigma(4); sigma_34inv = 1/sigma(3) + 1/sigma(4);

Z = [sigma_12*N11,          -sigma(2)*N12,      zeros(nv1,nv3),          zeros(nv1,nv4),        -2*D11_compl,                D12_compl,                    zeros(nv1, nc3);
     -sigma(2)*N21,    sigma_23*N22,          -sigma(3)*N23,       zeros(nv2,nv4),          D21_compl,               -2*D22_compl,                       D23_compl;
     zeros(nv3, nv1),    -sigma(3)*N32,        sigma_34*N33,        -sigma(4)*N34,      zeros(nv3,nc1),               D32_compl,                     -2*D33_compl;
     zeros(nv4,nv1),       zeros(nv4,nv2),        -sigma(4)*N43,        sigma(4)*N44,        zeros(nv4,nc1),             zeros(nv4,nc2),                    D43_compl;
     -2*D11,                  D12,              zeros(nc1,nv3),          zeros(nc1,nv4),        sigma_12inv*S11,         (-1/sigma(2))*S12,             zeros(nc1, nc3);
     D21,                 -2*D22,                  D23,               zeros(nc2,nv4),      (-1/sigma(2))*S21,     sigma_23inv*S22,          (-1/sigma(3))*S23;
     zeros(nc3,nv1),                D32,                 -2*D33,                    D34,               zeros(nc3,nc1),       (-1/sigma(3))*S32,         sigma_34inv*S33];

nZGram = [zeros(nv1,1); zeros(nv2,1); zeros(nv3,1); GPyramid4'*ones(nv4,1); zeros(nc1,1); zeros(nc2,1); zeros(nc3,1)];
Z = Z + nZGram;

QDiag = [ones(nv1, 1) / sqrt(max(sigma(1),sigma(2)));
    ones(nv2, 1) / sqrt(max(sigma(2),sigma(3)));
    ones(nv3, 1) / sqrt(max(sigma(3),sigma(4)));
    ones(nv4, 1) / sqrt(sigma(4));
    ones(nc1, 1) * sqrt(min(sigma(1),sigma(2)));
    ones(nc2, 1) * sqrt(min(sigma(2),sigma(3)));
    ones(nc3, 1) * sqrt(min(sigma(3),sigma(4)))];
Q = spdiags(QDiag, 0, length(QDiag), length(QDiag));

Z_inv = pinv(full(Q * Z * Q));

%% Forward Model
G_four_layers = zeros(nv4, size(posDipole_list,2));

for i = 1:length(good_idx)
    n = good_idx(i);
    posDipole = posDipole_list(:,n)'; 
    momDipole = momDipole_list(:, n)'; 
    
    rhs_GradNPotPyramid_An1 = zeros(nv1, 1);
    rhs_PotPatch1 = zeros(nc1, 1);
    
    [rhs_GradNPotPyramid_An1, ~] = calllib(library_name, 'computeGradNPotPyramidRHSAnal', rhs_GradNPotPyramid_An1, zeros(nv1, 1), funcSpacePyramid1, posDipole, momDipole);
    [rhs_PotPatch1, ~] = calllib(library_name, 'computePotPatchRHSNum', rhs_PotPatch1, zeros(nc1, 1), funcSpacePatch1, posDipole, momDipole);
   
    RHS = [+rhs_GradNPotPyramid_An1; zeros(nv2+nv3+nv4, 1); -rhs_PotPatch1 / sigma(1); zeros(nc2, 1); zeros(nc3, 1)];
    
    y = Z_inv * (Q * RHS); 
    solcond = Q * y;
    
    Pot_scalp = solcond(nv1+nv2+nv3+1 : nvtot);
    
    G_four_layers(:, n) = Pot_scalp - mean(Pot_scalp);
    
    if mod(i, 1000) == 0
        fprintf('Sorgente %d/%d calcolata\n', n, size(posDipole_list,2)); 
    end
end

%% Reordering
TVB_scalp_vertices = surf_struct.skin_air.vertices'; 
[idx_mapping, dist] = knnsearch(points4, TVB_scalp_vertices);

if max(dist) > 1e-3
    warning('Attenzione: Mismatch geometrico > 1mm tra la mesh C++ e TVB!');
end
G_four_layers = G_four_layers(idx_mapping, :);

%% Save G
fileName = 'G_four_layers.mat';
fullPath = fullfile(mat_dir, fileName);

save(fullPath, 'G_four_layers', '-v7.3');