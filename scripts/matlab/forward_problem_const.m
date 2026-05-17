clear all
clc

%% Create the structure if it does not exist and import other files

dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
h5_dir  = fullfile(dataset_dir, 'h5_files');
mat_dir = fullfile(dataset_dir, 'mat_files');

results_dir = fullfile(pwd,'simulations','matlab');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

struct_name = "Surfaces_structure.mat";
struct_path = fullfile(mat_dir, struct_name);

if isfile(struct_path)
    load(struct_path);
else
    surf_struct = h5_load_surf(h5_dir);
    if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end
    save(struct_path, "surf_struct");
end

sensor_name = "Sensors_EEG.h5";
sensor_path = fullfile(h5_dir,sensor_name);

if isfile(sensor_path)
    sensors = h5read(sensor_path,"/locations");
else
    error("Sensors files doesn't exist or is not found %s",sensor_name);
end

%% Parameters

sig_brain = 1.0;
sig_skull = 0.01;
sig_scalp =  1.0;
sig_air = 0.0;

%in order to simplify the code we insert in order of depth (first the
%deeper)

sig_in = [sig_brain, sig_skull, sig_scalp];
sig_out = [sig_skull, sig_scalp, sig_air];

%% Faces centers and areas calculation

%We extract the vertices and triangles

% 1. brain_skull
bnd(1).v = surf_struct.brain_skull.vertices'; 
bnd(1).f = surf_struct.brain_skull.triangles' + 1;

% 2. skull_skin
bnd(2).v = surf_struct.skull_skin.vertices'; 
bnd(2).f = surf_struct.skull_skin.triangles' + 1;

% 3. scalp_air
bnd(3).v = surf_struct.skin_air.vertices'; 
bnd(3).f = surf_struct.skin_air.triangles' + 1;

n_surf = length(bnd); %number of surfaces
n_total = 0; %total number of triangles

for i=1:n_surf
    %we first find the triangles centers

    p1 = bnd(i).v(bnd(i).f(:, 1), :);
    p2 = bnd(i).v(bnd(i).f(:, 2), :);
    p3 = bnd(i).v(bnd(i).f(:, 3), :);

    bnd(i).centers = (p1 + p2 + p3)/3;

    %triangles counter

    bnd(i).n_tri = size(bnd(i).f, 1);
    n_total = n_total + bnd(i).n_tri;
end

%vectors that contains all vertex in order of depth togheter
all_centers = vertcat(bnd.centers);

%we create vectors with conductivity for every triangle
sigma_minus = zeros(n_total, 1);
sigma_plus  = zeros(n_total, 1);

idx = 1;
for i = 1:n_surf
    last_idx = idx + bnd(i).n_tri - 1;
    sigma_minus(idx:last_idx) = sig_in(i);
    sigma_plus(idx:last_idx)  = sig_out(i);
    idx = last_idx + 1;
end

all_v1 = zeros(n_total, 3);
all_v2 = zeros(n_total, 3);
all_v3 = zeros(n_total, 3);

idx_v = 1;
for i = 1:n_surf
    p1 = bnd(i).v(bnd(i).f(:, 1), :);
    p2 = bnd(i).v(bnd(i).f(:, 2), :);
    p3 = bnd(i).v(bnd(i).f(:, 3), :);
    
    n_t = bnd(i).n_tri;
    all_v1(idx_v:idx_v+n_t-1, :) = p1;
    all_v2(idx_v:idx_v+n_t-1, :) = p2;
    all_v3(idx_v:idx_v+n_t-1, :) = p3;
    
    idx_v = idx_v + n_t;
end

% assemble B matrix (solid angle)

B = zeros(n_total,n_total);

num_sigma = (sigma_minus - sigma_plus)';
den_sigma = sigma_minus + sigma_plus;

for k = 1:n_total
    r_k = all_centers(k, :);
    
    R1 = all_v1 - r_k;
    R2 = all_v2 - r_k;
    R3 = all_v3 - r_k;
    
    R1_norm = vecnorm(R1, 2, 2);
    R2_norm = vecnorm(R2, 2, 2);
    R3_norm = vecnorm(R3, 2, 2);

    num = sum(R1 .* cross(R2, R3, 2), 2);
    
    dot12 = sum(R1 .* R2, 2);
    dot23 = sum(R2 .* R3, 2);
    dot31 = sum(R3 .* R1, 2);
    den = R1_norm .* R2_norm .* R3_norm + dot12 .* R3_norm + dot23 .* R1_norm + dot31 .* R2_norm;

    solid_angle = 2 * atan2(num, den);
    
    solid_angle(k) = 0; 
    
    B(k, :) = (1 / (2*pi)) * (num_sigma ./ den_sigma(k)) .* solid_angle';
end

%deflation

C = B - ones(n_total)/n_total;

A = eye(n_total)-C;

%% K matrix calculation

src_pos = surf_struct.cortical_surface.vertices'; 
src_dir = surf_struct.cortical_surface.vertex_normals'; 

n_dipoles = size(src_pos, 1);
K = zeros(n_total, n_dipoles);

for k = 1:n_total
    r_k = all_centers(k, :);
    
    R_vec = r_k - src_pos; 
    R_norm = vecnorm(R_vec, 2, 2);
   
    R_norm(R_norm < 1e-10) = Inf; 
    
    dipole_potential = sum(src_dir .* R_vec, 2) ./ ((R_norm.^3));

    K(k, :) = (1 / (2*pi*den_sigma(k))) * dipole_potential';
end

%% S constant

n_scalp_faces = bnd(3).n_tri;
n_brain_skull_faces = bnd(1).n_tri + bnd(2).n_tri;

S_faces = [zeros(n_scalp_faces, n_brain_skull_faces), eye(n_scalp_faces)];

%%  W (from face to vertex)

n_scalp_vertices = size(bnd(3).v, 1);
W_scalp = sparse(n_scalp_vertices, n_scalp_faces);

for i = 1:n_scalp_faces
    v1 = bnd(3).f(i, 1);
    v2 = bnd(3).f(i, 2);
    v3 = bnd(3).f(i, 3);
    
    W_scalp(v1, i) = 1;
    W_scalp(v2, i) = 1;
    W_scalp(v3, i) = 1;
end

W_scalp = W_scalp ./ sum(W_scalp, 2);

%% G final pass

V_mesh = A \ K;

G_faces = S_faces * V_mesh;

G_first_BEM = W_scalp * G_faces;

G_first_BEM = G_first_BEM - mean(G_first_BEM, 1);

%% Save G matrix

G_name = "First_BEM.mat";
G_path = fullfile(mat_dir, G_name);
save(G_path, 'G_first_BEM');

%% A visual

fig1 = figure;
imagesc(abs(A));
colormap('parula');
cb = colorbar;
ylabel(cb, 'Absolute Value Magnitude');

% Add lines to highlight the 3 bnd (Brain, Skull, Scalp)
hold on;
b1_idx = size(bnd(1).f,1);
b2_idx = size(bnd(1).f,1) + size(bnd(2).f,1);

% Vertical separators
xline(b1_idx, 'w--', 'LineWidth', 1.5);
xline(b2_idx, 'w--', 'LineWidth', 1.5);
% Horizontal separators
yline(b1_idx, 'w--', 'LineWidth', 1.5);
yline(b2_idx, 'w--', 'LineWidth', 1.5);
hold off;

title('Structure of the BEM Matrix A (Constant Elements)');
xlabel('Triangle Index (Observation)');
ylabel('Triangle Index (Integration)');
axis square;

exportgraphics(fig1, fullfile(results_dir,'Const_BEM_A_visual.pdf'), 'ContentType', 'vector');