%% Load Data
dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
h5_dir  = fullfile(dataset_dir, 'h5_files');
mat_dir = fullfile(dataset_dir, 'mat_files');
results_dir = fullfile(pwd,'simulations','matlab');

struct_name = "Surfaces_structure.mat";
struct_path = fullfile(mat_dir, struct_name);

My_G_name = "First_BEM.mat";
My_G_path = fullfile(mat_dir, My_G_name);

Three_layer_G_name = "G_three_layers.mat";
Three_layer_G_path = fullfile(mat_dir, Three_layer_G_name);

G_tvb_name = "ProjectionMatrix_EEG.h5";
G_tvb_path = fullfile(h5_dir, G_tvb_name);
G_tvb = h5read(G_tvb_path, "/projection_data");
G_tvb = G_tvb';

if isfile(struct_path)
    load(struct_path);
else
    error("Genera prima Surfaces_structure.mat");
end

load(My_G_path, 'G_first_BEM'); 
load(Three_layer_G_path, 'G_three_layers'); 

sensors = surf_struct.sensors;

%% Mapping to sensors

bnd_v = surf_struct.skin_air.vertices; % [3 x 1082]
bnd_f = surf_struct.skin_air.triangles' + 1;

n_sensors = size(sensors, 2);
n_vertices = size(bnd_v, 2);
W_map = zeros(n_sensors, n_vertices);

for i = 1:n_sensors
    dists = vecnorm(bnd_v - sensors(:, i), 2, 1);
    
    [sorted_dists, idx] = sort(dists);
    k_nearest = 3; %we take the k nearest vertices to the sensor
    nearest_idx = idx(1:k_nearest);
    
    weights = 1 ./ sorted_dists(1:k_nearest);
    weights = weights / sum(weights); % normalize
    
    W_map(i, nearest_idx) = weights;
end

G_my_sens  = W_map * G_first_BEM;
G_sym_sens = W_map * G_three_layers;

%% Source Setup
src_idx = 100; %index of activated source

% Scalp Geometry
bnd_v = surf_struct.skin_air.vertices;
bnd_f = surf_struct.skin_air.triangles' + 1;

% Source Position
x_src = surf_struct.cortical_surface.vertices(1, src_idx);
y_src = surf_struct.cortical_surface.vertices(2, src_idx);
z_src = surf_struct.cortical_surface.vertices(3, src_idx);

% Source Orientation (Normal Vector)
nx_src = surf_struct.cortical_surface.vertex_normals(1, src_idx);
ny_src = surf_struct.cortical_surface.vertex_normals(2, src_idx);
nz_src = surf_struct.cortical_surface.vertex_normals(3, src_idx);

% Normalize just to be safe
norm_mag = norm([nx_src, ny_src, nz_src]);
nx_src = nx_src / norm_mag;
ny_src = ny_src / norm_mag;
nz_src = nz_src / norm_mag;

%source signal
Fs = 5000;
t = 0 : 1/Fs : 0.1;
A = 1e-9;
freq_signal = 100;
D_real = A * sin(2*pi*freq_signal*t); 

%% Solution
V_my  = G_my_sens(:, src_idx) * D_real;
V_sym_3 = G_sym_sens(:, src_idx) * D_real;
V_tvb = G_tvb(:, src_idx) * D_real;

V_my  = V_my  - mean(V_my, 1);
V_sym_3 = V_sym_3 - mean(V_sym_3, 1);
V_sym_3 = -V_sym_3;
V_tvb = V_tvb - mean(V_tvb, 1);

%% Plot temp series
[~, best_sensor] = max(var(V_tvb, 0, 2));

fig = figure;
fig.WindowState = 'maximized';
plot(t, V_my(best_sensor, :), 'b', 'LineWidth', 1.5); hold on;
plot(t, V_sym_3(best_sensor, :), 'r--', 'LineWidth', 1.5);
plot(t, V_tvb(best_sensor, :), 'k:', 'LineWidth', 1.5);
title(sprintf('Best sensor signal with single activation'));
xlabel('Time (s)'); ylabel('Pot (V)');
legend('Basic BEM', 'Symmetric 3-Layer BEM', 'TVB');
grid on;
exportgraphics(gcf, fullfile(results_dir,'G_formulations_comparison.pdf'), 'ContentType', 'vector');

%% 3D Spatial Maps

[~, peak_idx] = max(V_sym_3(best_sensor, :));

V_my_peak  = V_my(:, peak_idx);
V_sym_3_peak = V_sym_3(:, peak_idx);
V_tvb_peak = V_tvb(:, peak_idx);

c_max = max([max(abs(V_my_peak)), max(abs(V_sym_3_peak)), max(abs(V_tvb_peak))]);
c_min = -c_max;
arrow_scale = 0.03; 

% --- plot 1: First BEM ---
fig = figure;
fig.WindowState = 'maximized';

annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
    sprintf('Time peak Basic BEM'), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

subplot(1,2,1);
trisurf(bnd_f, bnd_v(1,:), bnd_v(2,:), bnd_v(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
hold on; axis equal; rotate3d on;
scatter3(sensors(1,:), sensors(2,:), sensors(3,:), 150, V_my_peak, 'filled', 'MarkerEdgeColor', 'k');
scatter3(x_src, y_src, z_src, 100, 'k', 'filled', 'p', 'MarkerEdgeColor', 'w'); 
quiver3(x_src, y_src, z_src, nx_src, ny_src, nz_src, arrow_scale, 'Color', 'k', 'LineWidth', 3, 'MaxHeadSize', 2);
colormap(jet); clim([c_min, c_max]); colorbar; title('Basic BEM top view'); view(2);

subplot(1,2,2);
trisurf(bnd_f, bnd_v(1,:), bnd_v(2,:), bnd_v(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
hold on; axis equal; rotate3d on;
scatter3(sensors(1,:), sensors(2,:), sensors(3,:), 150, V_my_peak, 'filled', 'MarkerEdgeColor', 'k');
scatter3(x_src, y_src, z_src, 100, 'k', 'filled', 'p', 'MarkerEdgeColor', 'w'); 
quiver3(x_src, y_src, z_src, nx_src, ny_src, nz_src, arrow_scale, 'Color', 'k', 'LineWidth', 3, 'MaxHeadSize', 2);
colormap(jet); clim([c_min, c_max]); colorbar; title('Basic BEM'); view(3);

exportgraphics(gcf, fullfile(results_dir,'G_formulations_comparison_3D_firstBEM.pdf'), 'ContentType', 'vector');

% --- plot 2: Symmetric BEM ---
fig = figure;
fig.WindowState = 'maximized';

annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
    sprintf('Time peak Symmetric 3-Layer BEM'), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

subplot(1,2,1);
trisurf(bnd_f, bnd_v(1,:), bnd_v(2,:), bnd_v(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
hold on; axis equal; rotate3d on;
scatter3(sensors(1,:), sensors(2,:), sensors(3,:), 150, V_sym_3_peak, 'filled', 'MarkerEdgeColor', 'k');
scatter3(x_src, y_src, z_src, 100, 'k', 'filled', 'p', 'MarkerEdgeColor', 'w');
quiver3(x_src, y_src, z_src, nx_src, ny_src, nz_src, arrow_scale, 'Color', 'k', 'LineWidth', 3, 'MaxHeadSize', 2);
colormap(jet); clim([c_min, c_max]); colorbar; title('Symmetric 3-Layer BEM top view'); view(2);

subplot(1,2,2);
trisurf(bnd_f, bnd_v(1,:), bnd_v(2,:), bnd_v(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
hold on; axis equal; rotate3d on;
scatter3(sensors(1,:), sensors(2,:), sensors(3,:), 150, V_sym_3_peak, 'filled', 'MarkerEdgeColor', 'k');
scatter3(x_src, y_src, z_src, 100, 'k', 'filled', 'p', 'MarkerEdgeColor', 'w');
quiver3(x_src, y_src, z_src, nx_src, ny_src, nz_src, arrow_scale, 'Color', 'k', 'LineWidth', 3, 'MaxHeadSize', 2);
colormap(jet); clim([c_min, c_max]); colorbar; title('Symmetric 3-Layer BEM'); view(3);

exportgraphics(gcf, fullfile(results_dir,'G_formulations_comparison_3D_three_layer.pdf'), 'ContentType', 'vector');

% --- plot 3: TVB ---
fig = figure;
fig.WindowState = 'maximized';

annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
    sprintf('Time peak TVB'), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

subplot(1,2,1);
trisurf(bnd_f, bnd_v(1,:), bnd_v(2,:), bnd_v(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
hold on; axis equal; rotate3d on;
scatter3(sensors(1,:), sensors(2,:), sensors(3,:), 150, V_tvb_peak, 'filled', 'MarkerEdgeColor', 'k');
scatter3(x_src, y_src, z_src, 100, 'k', 'filled', 'p', 'MarkerEdgeColor', 'w');
quiver3(x_src, y_src, z_src, nx_src, ny_src, nz_src, arrow_scale, 'Color', 'k', 'LineWidth', 3, 'MaxHeadSize', 2);
colormap(jet); clim([c_min, c_max]); colorbar; title('TVB top view'); view(2);

subplot(1,2,2);
trisurf(bnd_f, bnd_v(1,:), bnd_v(2,:), bnd_v(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
hold on; axis equal; rotate3d on;
scatter3(sensors(1,:), sensors(2,:), sensors(3,:), 150, V_tvb_peak, 'filled', 'MarkerEdgeColor', 'k');
scatter3(x_src, y_src, z_src, 100, 'k', 'filled', 'p', 'MarkerEdgeColor', 'w');
quiver3(x_src, y_src, z_src, nx_src, ny_src, nz_src, arrow_scale, 'Color', 'k', 'LineWidth', 3, 'MaxHeadSize', 2);
colormap(jet); clim([c_min, c_max]); colorbar; title('TVB'); view(3);

exportgraphics(gcf, fullfile(results_dir,'G_formulations_comparison_3D_tvb.pdf'), 'ContentType', 'vector');

%% Errors
RDM_sym_vs_tvb = norm(V_sym_3_peak/norm(V_sym_3_peak) - V_tvb_peak/norm(V_tvb_peak));
MAG_sym_vs_tvb = norm(V_sym_3_peak) / norm(V_tvb_peak);

fprintf('\n--- Comparison with TVB Baseline (Source %d) ---\n', src_idx);
fprintf('Topogr. Error (RDM, 0 = same): %.4f\n', RDM_sym_vs_tvb);
fprintf('Amplitude relation (MAG, 1 = same): %.4f\n\n', MAG_sym_vs_tvb);