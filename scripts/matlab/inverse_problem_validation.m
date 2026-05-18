%% Load data
dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
h5_dir  = fullfile(dataset_dir, 'h5_files');
mat_dir = fullfile(dataset_dir, 'mat_files');
results_dir = fullfile(pwd,'simulations','matlab');

if ~exist(results_dir, 'dir'), mkdir(results_dir); end

struct_path = fullfile(mat_dir, "Surfaces_structure.mat");
Three_layer_G_path = fullfile(mat_dir, "G_three_layers.mat");
G_tvb_path = fullfile(h5_dir, "ProjectionMatrix_EEG.h5");

load(struct_path);
load(Three_layer_G_path, 'G_three_layers'); 
G_tvb = h5read(G_tvb_path, "/projection_data")';
sensors = surf_struct.sensors;

%% Map to sensors
bnd_v = surf_struct.skin_air.vertices; 
n_sensors = size(sensors, 2);
n_vertices = size(bnd_v, 2);

W_map = zeros(n_sensors, n_vertices);
for i = 1:n_sensors
    dists = vecnorm(bnd_v - sensors(:, i), 2, 1);
    [sorted_dists, idx] = sort(dists);
    k_nearest = 3; 
    nearest_idx = idx(1:k_nearest);
    weights = 1 ./ sorted_dists(1:k_nearest);
    W_map(i, nearest_idx) = weights / sum(weights);
end
G_sym_sens = W_map * G_three_layers; 

G_sym_sens = G_sym_sens - mean(G_sym_sens, 1);
G_tvb_centered = G_tvb - mean(G_tvb, 1);

%% Pre-calc

SNR_levels = [0, 5, 10, 15, 25]; 
n_snr = length(SNR_levels);

n_sources = size(surf_struct.cortical_surface.vertices, 2);
pos_sources = surf_struct.cortical_surface.vertices;
I = eye(size(G_sym_sens, 1));
method_names = {'PINV', 'MNE', 'WMNE', 'sLORETA'};

% Arrays to store median metrics for the final Trend Plots
median_ed1_trend = zeros(n_snr, 4);
median_ed2_trend = zeros(n_snr, 4);

% Prepare mesh data for the 3D plots
c_v = double(surf_struct.cortical_surface.vertices);
c_f = double(surf_struct.cortical_surface.triangles' + 1);
if size(c_v,2) ~= 3, c_v = c_v.'; end
if size(c_f,2) ~= 3, c_f = c_f.'; end
TR = triangulation(c_f, c_v);
VN = -vertexNormal(TR); 
view_angles = {[-90, 0], [90, 0], [0, 90], [0, -90]};
view_names = {'Rear', 'Front', 'Superior (Top)', 'Inferior (Bottom)'};

%% Calc loop

for snr_idx = 1:n_snr
    current_SNR = SNR_levels(snr_idx);
    
    T_pinv = pinv(full(G_sym_sens));
    
    A = G_sym_sens * G_sym_sens';
    alpha_opt_MNE = (trace(A) / n_sensors) * (10^(-current_SNR / 10));
    T_MNE = G_sym_sens' / (A + (alpha_opt_MNE * I));
    
    col_norms = vecnorm(G_sym_sens, 2, 1)'; 
    P = spdiags(1 ./ (col_norms.^2), 0, n_sources, n_sources);
    A_WMNE = G_sym_sens * P * G_sym_sens';
    alpha_opt_WMNE = (trace(A_WMNE) / n_sensors) * (10^(-current_SNR / 10));
    T_WMNE = (P * G_sym_sens') / (A_WMNE + (alpha_opt_WMNE * I));
    
    diag_R = sum(T_MNE .* G_sym_sens', 2); 
    
    % Initialize error arrays
    ed1_curr = zeros(n_sources, 4); 
    ed2_curr = zeros(n_sources, 4);

    % --- Source Loop ---
    for n = 1:n_sources
        M_clean = G_tvb_centered(:, n); 
        noise_power = var(M_clean) / (10^(current_SNR/10));
        noise = sqrt(noise_power) * randn(size(M_clean));
        M = M_clean + noise;
        
        real_source = pos_sources(:, n);
        dist_all = vecnorm(pos_sources - real_source, 2, 1)';
        
        % Inversions
        D_pinv_signed = T_pinv * M;
        D_MNE_signed = T_MNE * M;
        D_WMNE_signed = T_WMNE * M;
        D_sLOR_signed = D_MNE_signed ./ sqrt(diag_R);

        D_pinv = D_pinv_signed.^2;
        D_MNE = D_MNE_signed.^2;
        D_WMNE = D_WMNE_signed.^2;
        D_sLOR = D_sLOR_signed.^2;
        
        estimates = {D_pinv, D_MNE, D_WMNE, D_sLOR};
        
        for m = 1:4
            D_vec = estimates{m};
            [~, idx_max] = max(D_vec); 
            ed1_curr(n, m) = dist_all(idx_max);

            energy_sum = sum(D_vec);
            if energy_sum == 0, energy_sum = eps; end
            ed2_curr(n, m) = sum(dist_all .* (D_vec / energy_sum));
        end
        
        if mod(n, 2000) == 0
            fprintf(' Processed Source %d/%d\n', n, n_sources); 
        end
    end
    
    ed1_curr = ed1_curr * 1000;
    ed2_curr = ed2_curr * 1000;
    
    median_ed1_trend(snr_idx, :) = median(ed1_curr, 1);
    median_ed2_trend(snr_idx, :) = median(ed2_curr, 1);
    
    for m = 1:4
        curr_ed1 = ed1_curr(:, m);
        method_name = method_names{m};
        
        max_ed1_method = max(curr_ed1);
        if max_ed1_method < 1, max_ed1_method = 1; end
        
        fig = figure;
        fig.WindowState = 'maximized';
        
        annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
            sprintf('Method: %s (ED1, Mean: %.1f mm | Median: %.1f mm | Max: %.1f mm), SNR = %d dB', ...
            method_name, mean(curr_ed1), median(curr_ed1), max_ed1_method, current_SNR), ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
            
        for v = 1:4
            subplot(2, 2, v);
            
            plot_brain_surface(c_v, c_f, VN, log10(curr_ed1+1), view_names{v}, view_angles{v}, log10(max_ed1_method+1));
            
            if v == 4
                cb = colorbar;
                tick_candidates = [0, 1, 2, 5, 10, 20, 50, 100, 150];
                real_ticks = tick_candidates(tick_candidates < max_ed1_method);
                real_ticks = [real_ticks, round(max_ed1_method)];
                real_ticks = unique(real_ticks);
                
                cb.Ticks = log10(real_ticks + 1);
                cb.TickLabels = string(real_ticks);
                ylabel(cb, 'ED1 (mm)');
            end
        end
        tmp_filename = sprintf('Inverse_ED1_%s_SNR_%02d.pdf', method_name, current_SNR);
        exportgraphics(fig, fullfile(results_dir, tmp_filename), 'ContentType', 'image', 'Resolution', 600);
    end
end

%% Global Performance

fig_tr = figure;
fig_tr.WindowState = 'maximized';
colors = {'k', 'b', 'm', 'r'};
markers = {'o', 's', '^', 'd'};

% Trend for ED1 (Localization Error)
subplot(1,2,1);
hold on; grid on;
for m = 1:4
    plot(SNR_levels, median_ed1_trend(:, m), 'Color', colors{m}, ...
         'Marker', markers{m}, 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', method_names{m});
end
title('Localization Error');
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Median ED1 (mm)');
legend('Location', 'northeast');
set(gca, 'XDir', 'reverse'); % Show decreasing noise left-to-right

% Trend for ED2 (Spatial Dispersion)
subplot(1,2,2);
hold on; grid on;
for m = 1:4
    plot(SNR_levels, median_ed2_trend(:, m), 'Color', colors{m}, ...
         'Marker', markers{m}, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', method_names{m});
end
title('Spatial Dispersion');
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Median ED2 (mm)');
legend('Location', 'northeast');
set(gca, 'XDir', 'reverse');

exportgraphics(fig_tr, fullfile(results_dir, 'Performance_Trends_vs_SNR.pdf'), 'ContentType', 'vector');