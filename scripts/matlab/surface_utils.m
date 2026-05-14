clear all
close all
clc

dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
h5_dir  = fullfile(dataset_dir, 'h5_files');
mat_dir = fullfile(dataset_dir, 'mat_files');
stl_dir = fullfile(dataset_dir, 'stl_files_generated');
obj_dir = fullfile(dataset_dir, 'obj_files_generated');

if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end
if ~exist(stl_dir, 'dir'), mkdir(stl_dir); end
if ~exist(obj_dir, 'dir'), mkdir(obj_dir); end

%% open h5 files and create surface struct with conversion to meters

struct_name = "Surfaces_structure.mat";
struct_path = fullfile(mat_dir, struct_name);

surf_struct = h5_load_surf(h5_dir);

h5_to_stl(surf_struct, stl_dir);
h5_to_obj(surf_struct, obj_dir)

sensor_name = "Sensors_EEG.h5";
sensor_path = fullfile(h5_dir, sensor_name);

sensors_raw = h5read(sensor_path, "/locations");
surf_struct.sensors = sensors_raw / 1000;

save(struct_path, "surf_struct", "-v7.3");

%% surf check
%to see geometry and normals direction

vertices = surf_struct.cortical_surface.vertices;
mom_vertices = surf_struct.cortical_surface.vertex_normals;

figure;

scatter3(vertices(1,:), vertices(2,:), vertices(3,:), 36, 'b', '.');
hold on;
axis equal; 
grid on;

step = 1;
idx = 1:step:size(vertices, 2);

scale_factor = 1.0; 
quiver3(vertices(1,idx), vertices(2,idx), vertices(3,idx), ...
        mom_vertices(1,idx), mom_vertices(2,idx), mom_vertices(3,idx), ...
        scale_factor, 'r', 'MaxHeadSize', 0.5);

title('Vertices and Normals');
legend('Vertices ', 'Vertices Normals');
view(3);

%% cortical surface source scaling
%do the scaling of the cortical surface vertices in order to shift the
%dipoles from the boundary, the result has only the good vertices

dr = 0.002;

pos_orig = surf_struct.cortical_surface.vertices;
normals = surf_struct.cortical_surface.vertex_normals;
norm_magnitude = vecnorm(normals, 2, 1);
normals = normals ./ norm_magnitude;

bnd_f = surf_struct.cortical_surface.triangles' + 1;

pos_new = pos_orig - (dr .* normals);

[~, dist_nearest] = knnsearch(pos_orig', pos_new');

tollerance = dr * 0.05; 
bad_idx = find(dist_nearest < (dr - tollerance));

good_idx = setdiff(1:size(pos_new,2), bad_idx);

figure;

trisurf(bnd_f, pos_orig(1,:), pos_orig(2,:), pos_orig(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
hold on; axis equal; grid on; rotate3d on; %cortical surface

scatter3(pos_new(1, good_idx), pos_new(2, good_idx), pos_new(3, good_idx), ...
    10, 'g', 'filled', 'MarkerEdgeColor', 'k');

% if ~isempty(bad_idx)
%     scatter3(pos_new(1, bad_idx), pos_new(2, bad_idx), pos_new(3, bad_idx), ...
%         30, 'r', 'filled', 'MarkerEdgeColor', 'k');
% end

title(sprintf('Internal Shift (dr = %g m)', dr));
legend('Original Surface', 'Ok Vertices', 'Bad Vertices', 'Location', 'best');
view(3);