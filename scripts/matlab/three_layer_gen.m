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
brain_skull_msh = fullfile(msh_dir, 'Surface_brain_skull.msh');
brain_skull_gmsh = fullfile(msh_dir, 'Surface_brain_skull.gmsh');
skull_skin_msh = fullfile(msh_dir, 'Surface_skull_skin.msh');
skull_skin_gmsh = fullfile(msh_dir, 'Surface_skull_skin.gmsh');
skin_air_msh = fullfile(msh_dir, 'Surface_skin_air.msh');
skin_air_gmsh = fullfile(msh_dir, 'Surface_skin_air.gmsh');

struct_path = fullfile(mat_dir, 'Surfaces_structure.mat');
load(struct_path);

%% Create spaces

sigma = [1.0, 0.01, 1.0]; %same as projection matrix from TVB

[points1, cells1, normals1, dim1, nNode1] = readMshFile_Samuele(brain_skull_msh);
nPoints1 = size(points1,1);
nCells1 = size(cells1,1);

[points2, cells2, normals2, dim2, nNode2] = readMshFile_Samuele(skull_skin_msh);
nPoints2 = size(points2,1);
nCells2 = size(cells2,1);

[points3, cells3, normals3, dim3, nNode3] = readMshFile_Samuele(skin_air_msh);
nPoints3 = size(points3,1);
nCells3 = size(cells3,1);

if(exist(brain_skull_gmsh,'file')~=2)
    gmsh_mesh3d_write(brain_skull_gmsh, dim1, nPoints1, points1, nNode1, nCells1, cells1);
end

if(exist(skull_skin_gmsh,'file')~=2)
    gmsh_mesh3d_write(skull_skin_gmsh, dim2, nPoints2, points2, nNode2, nCells2, cells2);
end

if(exist(skin_air_gmsh,'file')~=2)
    gmsh_mesh3d_write(skin_air_gmsh, dim3, nPoints3, points3, nNode3, nCells3, cells3);
end

space1 = calllib(library_name, 'createSpace', brain_skull_gmsh);
space2 = calllib(library_name, 'createSpace', skull_skin_gmsh);
space3 = calllib(library_name, 'createSpace', skin_air_gmsh);

%% Data from spaces

funcSpacePatch1 = calllib(library_name, 'createFunctionalSpacePatch', space1);
funcSpacePatch2 = calllib(library_name, 'createFunctionalSpacePatch', space2);
funcSpacePatch3 = calllib(library_name, 'createFunctionalSpacePatch', space3);


funcSpacePyramid1 = calllib(library_name, 'createFunctionalSpacePyramid', space1);
funcSpacePyramid2 = calllib(library_name, 'createFunctionalSpacePyramid', space2);
funcSpacePyramid3 = calllib(library_name, 'createFunctionalSpacePyramid', space3);


funcSpaceGradPyramid1 = calllib(library_name, 'createFunctionalSpaceGradPyramid', space1);
funcSpaceGradPyramid2 = calllib(library_name, 'createFunctionalSpaceGradPyramid', space2);
funcSpaceGradPyramid3 = calllib(library_name, 'createFunctionalSpaceGradPyramid', space3);

nc1 = calllib(library_name, 'getNbPatchFunctions', funcSpacePatch1);
nc2 = calllib(library_name, 'getNbPatchFunctions', funcSpacePatch2);
nc3 = calllib(library_name, 'getNbPatchFunctions', funcSpacePatch3);

nv1 = calllib(library_name, 'getNbPyramidFunctions', funcSpacePyramid1);
nv2 = calllib(library_name, 'getNbPyramidFunctions', funcSpacePyramid2);
nv3 = calllib(library_name, 'getNbPyramidFunctions', funcSpacePyramid3);


nctot = nc1 + nc2 + nc3;
nvtot = nv1 + nv2 + nv3;

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

% h1
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
edges3 = computeEdges(cells3);
lengthE3 = getLengthEdges(edges3,points3);
lengthSum3 = sum(lengthE3);
h3 = lengthSum3/length(lengthE3);

hlist = (h1+h2+h3)/3;

N_list = nvtot+nctot; %number of unknowns

%% S matrices
if(exist('Ssph_real.mat','file')~=2)
    Z_R=zeros(nc1,nc1); Z_I=zeros(nc1,nc1);
    [Z_R, ~] = calllib(library_name, 'computeSPatch', Z_R, Z_I, funcSpacePatch1, funcSpacePatch1); S11 = Z_R;
    
    Z_R=zeros(nc1,nc2); Z_I=zeros(nc1,nc2);
    [Z_R, ~] = calllib(library_name, 'computeSPatch', Z_R, Z_I, funcSpacePatch1, funcSpacePatch2); S12 = Z_R;
    S21 = S12';
    
    Z_R=zeros(nc2,nc2); Z_I=zeros(nc2,nc2);
    [Z_R, ~] = calllib(library_name, 'computeSPatch', Z_R, Z_I, funcSpacePatch2, funcSpacePatch2); S22 = Z_R;
    save('Ssph_real.mat', 'S11','S12','S21','S22');
else
    load('Ssph_real.mat');
end

%% N matrices
if(exist('Nsph_real.mat','file')~=2)
    Z_R=zeros(nv1,nv1); Z_I=zeros(nv1,nv1);
    [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid1, funcSpacePyramid1); N11 = Z_R;
    
    Z_R=zeros(nv1,nv2); Z_I=zeros(nv1,nv2);
    [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid1, funcSpacePyramid2); N12 = Z_R;
    N21 = N12';
    
    Z_R=zeros(nv2,nv2); Z_I=zeros(nv2,nv2);
    [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid2, funcSpacePyramid2); N22 = Z_R;
    
    Z_R=zeros(nv2,nv3); Z_I=zeros(nv2,nv3);
    [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid2, funcSpacePyramid3); N23 = Z_R;
    N32 = N23';
    
    Z_R=zeros(nv3,nv3); Z_I=zeros(nv3,nv3);
    [Z_R, ~] = calllib(library_name, 'computeNPyramid', Z_R, Z_I, funcSpacePyramid3, funcSpacePyramid3); N33 = Z_R;
    save('Nsph_real.mat', 'N11','N12','N21','N22','N23','N32','N33');
else
    load('Nsph_real.mat');
end

%% D matrices
if(exist('Dsph_real.mat','file')~=2)
    Z_R=zeros(nc1,nv1); Z_I=zeros(nc1,nv1);
    [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch1, funcSpacePyramid1); D11 = Z_R; D11_compl=D11';
    
    Z_R=zeros(nc1,nv2); Z_I=zeros(nc1,nv2);
    [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch1, funcSpacePyramid2); D12 = Z_R; D21_compl=D12';
    
    Z_R=zeros(nc2,nv1); Z_I=zeros(nc2,nv1);
    [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch2, funcSpacePyramid1); D21 = Z_R; D12_compl=D21';
    
    Z_R=zeros(nc2,nv2); Z_I=zeros(nc2,nv2);
    [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch2, funcSpacePyramid2); D22 = Z_R; D22_compl=D22';
    
    Z_R=zeros(nc2,nv3); Z_I=zeros(nc2,nv3);
    [Z_R, ~] = calllib(library_name, 'computeDPatchPyramid', Z_R, Z_I, funcSpacePatch2, funcSpacePyramid3); D23 = Z_R; D32_compl=D23';
    save('Dsph_real.mat', 'D11','D11_compl','D12','D21_compl','D21','D12_compl','D22','D22_compl','D23','D32_compl');
else
    load('Dsph_real.mat');
end

%% Compute G

pDipole_coords = surf_struct.cortical_surface.vertices; 
n_cortical_sources = size(pDipole_coords, 2);
pDipole_norms = surf_struct.cortical_surface.vertex_normals;

norm_magnitude = vecnorm(pDipole_norms, 2, 1);
pDipole_norms = pDipole_norms ./ norm_magnitude;

G_three_layers = zeros(nv3, n_cortical_sources);

sigma_12 = sigma(1) + sigma(2);
sigma_12inv = 1/sigma(1) + 1/sigma(2);
sigma_23 = sigma(2) + sigma(3);
sigma_23inv = 1/sigma(2) + 1/sigma(3);

Z = [sigma_12*N11,        -sigma(2)*N12,    zeros(nv1, nv3),  -2*D11_compl,        D12_compl;
     -sigma(2)*N21,        sigma_23*N22,    -sigma(3)*N23,     D21_compl,       -2*D22_compl;
     zeros(nv3,nv1),      -sigma(3)*N32,     sigma(3)*N33,     zeros(nv3,nc1),     D32_compl;
     -2*D11,               D12,              zeros(nc1,nv3),   sigma_12inv*S11,   -S12/sigma(2);
     D21,                 -2*D22,            D23,             -S21/sigma(2),       sigma_23inv*S22];


Z_inv = pinv(full(Z));

% Forward Solution
for n = 1:n_cortical_sources
    
    posDipole = pDipole_coords(:, n)'; 
    momDipole = pDipole_norms(:, n)'; 
   
    rhs_GradNPotPyramid_An1 = zeros(nv1, 1);
    rhs_PotPatch1 = zeros(nc1, 1);
    imagPartv1 = zeros(nv1, 1);
    imagPartc1 = zeros(nc1, 1);

    [rhs_GradNPotPyramid_An1, ~] = calllib(library_name, ...
        'computeGradNPotPyramidRHSAnal', rhs_GradNPotPyramid_An1, ...
        imagPartv1, funcSpacePyramid1, posDipole, momDipole);
        
    [rhs_PotPatch1, ~] = calllib(library_name, ...
        'computePotPatchRHSNum', rhs_PotPatch1, ...
        imagPartc1, funcSpacePatch1, posDipole, momDipole);
        
    RHS = [+rhs_GradNPotPyramid_An1; zeros(nv2+nv3, 1); ...
           -rhs_PotPatch1/sigma(1); zeros(nc2, 1)];

    unknowns = Z_inv * RHS;
    
    Pot_omega3 = unknowns(1+nv1+nv2 : nvtot);
    
    G_three_layers(:, n) = Pot_omega3 - mean(Pot_omega3);
    
    if mod(n, 1000) == 0, fprintf('Sorgente %d/%d calcolata\n', n, n_cortical_sources); end
end

%% reorder vertices and save G

TVB_scalp_vertices = surf_struct.skin_air.vertices'; 
[idx_mapping, dist] = knnsearch(points3, TVB_scalp_vertices);

if max(dist) > 1e-3
    warning('Mismatch geometrico tra la mesh C++ e TVB!');
end

G_three_layers = G_three_layers(idx_mapping, :);

fileName = 'G_three_layers.mat';
fullPath = fullfile(mat_dir, fileName);

save(fullPath, 'G_three_layers', '-v7.3');