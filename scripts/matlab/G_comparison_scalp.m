%% Load Data
dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
h5_dir  = fullfile(dataset_dir, 'h5_files');
mat_dir = fullfile(dataset_dir, 'mat_files');
results_dir = fullfile(pwd,'simulations','matlab');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

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
bnd_v = surf_struct.skin_air.vertices; 
bnd_f = surf_struct.skin_air.triangles' + 1;

n_sensors = size(sensors, 2);
n_vertices = size(bnd_v, 2);
W_map = zeros(n_sensors, n_vertices);

for i = 1:n_sensors
    % Calculate distances from sensor 'i' to all scalp vertices
    dists = vecnorm(bnd_v - sensors(:, i), 2, 1);
    
    [sorted_dists, idx] = sort(dists);
    k_nearest = 3; % Inverse Distance Weighting using 3 nearest neighbors
    nearest_idx = idx(1:k_nearest);
    
    weights = 1 ./ sorted_dists(1:k_nearest);
    weights = weights / sum(weights); % Normalize weights to sum to 1
    
    W_map(i, nearest_idx) = weights;
end

G_my_sens  = W_map * G_first_BEM;
G_sym_sens = W_map * G_three_layers;

%% Source Setup
src_idx = 500; %index of activated source

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
n_vec = [nx_src, ny_src, nz_src] / norm_mag;

% Generate Source Signal
Fs = 5000;
t = 0 : 1/Fs : 0.1;
amplitude = 1e-9; % 1 nA*m
freq_signal = 100;
dipole_moment = amplitude * sin(2*pi*freq_signal*t);

%% Solution
V_my = G_my_sens(:, src_idx) * dipole_moment;
V_sym_3 = G_sym_sens(:, src_idx) * dipole_moment;
V_tvb = G_tvb(:, src_idx) * dipole_moment;

V_my = V_my - mean(V_my, 1);
V_sym_3 = V_sym_3 - mean(V_sym_3, 1);
V_tvb = V_tvb - mean(V_tvb, 1);

%% errors

[~, peak_idx] = max(V_sym_3(1, :));
V_my_peak = V_my(:, peak_idx);
V_sym_3_peak = V_sym_3(:, peak_idx);
V_tvb_peak = V_tvb(:, peak_idx);

% Relative Difference Measure (Topography Error)
RDM_sym_vs_tvb = norm(V_sym_3_peak/norm(V_sym_3_peak) - V_tvb_peak/norm(V_tvb_peak));
% Magnitude Error
MAG_sym_vs_tvb = norm(V_sym_3_peak) / norm(V_tvb_peak);

%% Plot temp series
[~, best_sensor] = max(var(V_tvb, 0, 2));

fig1 = figure;
plot(t, V_my(best_sensor, :), 'b', 'LineWidth', 1.5); hold on;
plot(t, V_sym_3(best_sensor, :), 'r--', 'LineWidth', 1.5);
plot(t, V_tvb(best_sensor, :), 'k:', 'LineWidth', 1.5);
title(sprintf('Sensor %d Signal Comparison (Dipole %d)', best_sensor, src_idx));
xlabel('Time (s)'); ylabel('Potential (V)');
legend('Constant BEM', 'Symmetric BEM', 'TVB (Reference)', 'Location', 'best');
grid on;
exportgraphics(fig1, fullfile(results_dir,'Timeseries_Comparison.pdf'), 'ContentType', 'vector');

%% 3D Topographic Maps
c_max = max([max(abs(V_my_peak)), max(abs(V_sym_3_peak)), max(abs(V_tvb_peak))]);
c_min = -c_max;

fig2 = figure('Name', 'Topographic Comparison (Smooth)', 'Color', 'w', 'Position', [100, 100, 1400, 450]);

peak_data = {V_my_peak, V_sym_3_peak, V_tvb_peak};
titles = {'Constant BEM', 'Symmetric BEM', 'TVB (Reference)'};

for i = 1:3
    subplot(1, 3, i);
    
    F = scatteredInterpolant(sensors(1,:)', sensors(2,:)', sensors(3,:)', ...
                             peak_data{i}, 'natural', 'nearest');
    
    V_scalp_interp = F(bnd_v(1,:)', bnd_v(2,:)', bnd_v(3,:)');
    
    patch('Vertices', bnd_v', 'Faces', bnd_f, ...
          'FaceVertexCData', V_scalp_interp, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 0.95);
    
    hold on; axis equal off;
    
    % Overlay the actual sensor positions as tiny black dots for reference
    scatter3(sensors(1,:), sensors(2,:), sensors(3,:), 15, 'k', 'filled');
    

    scatter3(x_src, y_src, z_src, 120, 'w', 'filled', 'p', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    quiver3(x_src, y_src, z_src, n_vec(1), n_vec(2), n_vec(3), 0.02, ...
            'Color', 'w', 'LineWidth', 3, 'MaxHeadSize', 2);
    quiver3(x_src, y_src, z_src, n_vec(1), n_vec(2), n_vec(3), 0.02, ...
            'Color', 'k', 'LineWidth', 1.5, 'MaxHeadSize', 2); % Black border for contrast
    
    colormap('jet'); 
    clim([c_min, c_max]); 
    title(titles{i}, 'FontSize', 14, 'FontWeight', 'bold'); 
    
    view(-20, 30); 

    material dull; 
    camlight('headlight'); 
    camlight('left');
    lighting gouraud;
end

cb = colorbar('southoutside');
cb.Position = [0.3 0.08 0.4 0.03];
cb.Label.String = 'Electrical Potential (V)';
cb.Label.FontSize = 12;

exportgraphics(fig2, fullfile(results_dir,'Topography_Comparison_Unified.pdf'), 'ContentType', 'vector');