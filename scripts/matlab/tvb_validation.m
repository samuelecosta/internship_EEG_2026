clear all; close all; clc;

%% Load Data
dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
h5_dir  = fullfile(dataset_dir, 'h5_files');
mat_dir = fullfile(dataset_dir, 'mat_files');
results_dir = fullfile(pwd, 'simulations', 'matlab');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

load(fullfile(mat_dir, "Surfaces_structure.mat"), 'surf_struct');
load(fullfile(mat_dir, "G_three_layers.mat"), 'G_three_layers'); 
G_tvb = h5read(fullfile(h5_dir, "ProjectionMatrix_EEG.h5"), "/projection_data")';

tvb_data = load(fullfile(mat_dir, "tvb_validation_data.mat"));
J_true_raw = double(tvb_data.J_true); 
time = double(tvb_data.time);
n_time = length(time);

J_true_raw = J_true_raw - mean(J_true_raw, 2);
J_true_norm_factor = max(abs(J_true_raw(:)));
if J_true_norm_factor > 0
    J_true_raw = J_true_raw / J_true_norm_factor;
end

sensors = double(surf_struct.sensors);
bnd_v = double(surf_struct.skin_air.vertices); 
pos_sources = double(surf_struct.cortical_surface.vertices);

n_sensors = size(sensors, 2);
n_vertices = size(bnd_v, 2);
n_sources = size(pos_sources, 2);

%% Spatial Alignment
tvb_vertices = tvb_data.tvb_vertices; 
if size(tvb_vertices, 1) ~= 3, tvb_vertices = tvb_vertices'; end
if max(abs(tvb_vertices(:))) > 10 && max(abs(pos_sources(:))) < 1
    tvb_vertices = tvb_vertices / 1000; 
end

[idx_map_to_tvb, dists] = knnsearch(tvb_vertices', pos_sources');
J_true_mapped = J_true_raw(idx_map_to_tvb, :);

if isfield(tvb_data, 'activation_center')
    target_idx_py = double(tvb_data.activation_center) + 1; 
    true_pos = tvb_vertices(:, target_idx_py);
    [~, target_idx] = min(vecnorm(pos_sources - true_pos, 2, 1));
else
    target_idx = 1; 
end

activation_center_pos = pos_sources(:, target_idx);
dist_from_target = vecnorm(pos_sources - activation_center_pos, 2, 1)'; 

%% Forward Problem
G_tvb_centered = G_tvb - mean(G_tvb, 1);
M_clean = G_tvb_centered * J_true_raw; 
[~, peak_time_idx] = max(abs(J_true_mapped(target_idx, :)));

%% Sensors Mapping
W_map = zeros(n_sensors, n_vertices);
for i = 1:n_sensors
    dists_sens = vecnorm(bnd_v - sensors(:, i), 2, 1);
    [sorted_dists, idx] = sort(dists_sens);
    k_nearest = 3; 
    nearest_idx = idx(1:k_nearest);
    weights = 1 ./ sorted_dists(1:k_nearest);
    W_map(i, nearest_idx) = weights / sum(weights);
end
G_sym_sens = W_map * G_three_layers; 
G_sym_sens = G_sym_sens - mean(G_sym_sens, 1);

%% Validation Loop (sLORETA)
SNR_levels = [25, 15, 10, 5, 0]; 
n_snr = length(SNR_levels);

n_trials = 50; %monte carlo

ED1_trend = zeros(n_snr, 1); ED2_trend = zeros(n_snr, 1); Temp_Corr_trend = zeros(n_snr, 1);
% Arrays per limiti di errore (IQR)
err_neg_ED1 = zeros(n_snr, 1); err_pos_ED1 = zeros(n_snr, 1);
err_neg_ED2 = zeros(n_snr, 1); err_pos_ED2 = zeros(n_snr, 1);
err_neg_Corr = zeros(n_snr, 1); err_pos_Corr = zeros(n_snr, 1);

I = eye(size(G_sym_sens, 1));
c_v = double(pos_sources);
c_f = double(surf_struct.cortical_surface.triangles' + 1);
if size(c_v,2) ~= 3, c_v = c_v.'; end
if size(c_f,2) ~= 3, c_f = c_f.'; end
TR = triangulation(c_f, c_v);
VN = -vertexNormal(TR); 

rng('default');

for snr_idx = 1:n_snr
    current_SNR = SNR_levels(snr_idx);
    
    A = G_sym_sens * G_sym_sens';
    alpha_opt = (trace(A) / n_sensors) * (10^(-current_SNR / 10));
    T_MNE = G_sym_sens' / (A + (alpha_opt * I));
    diag_R = sum(T_MNE .* G_sym_sens', 2);
    diag_R(diag_R < eps) = eps; 
    T_sLOR = T_MNE ./ sqrt(diag_R);

    t_ed1 = zeros(n_trials, 1); t_ed2 = zeros(n_trials, 1); t_corr = zeros(n_trials, 1);
    
    for trial = 1:n_trials
        noise_power = var(M_clean(:)) / (10^(current_SNR/10));
        noise = sqrt(noise_power) * randn(size(M_clean));
        M_noisy = M_clean + noise;
        
        % Inversione
        J_est = -(T_sLOR * M_noisy);
        J_est = J_est - mean(J_est, 2); 
        J_est_norm_factor = max(abs(J_est(:)));
        if J_est_norm_factor > 0
            J_est = J_est / J_est_norm_factor;
        end
        
        J_est_peak = abs(J_est(:, peak_time_idx));
        [~, est_max_idx] = max(J_est_peak);
        
        t_ed1(trial) = dist_from_target(est_max_idx) * 1000; 
        
        energy_sum = sum(J_est_peak.^2);
        if energy_sum == 0, energy_sum = eps; end
        t_ed2(trial) = sum(dist_from_target .* (J_est_peak.^2 / energy_sum)) * 1000;
        
        ts_true = J_true_mapped(target_idx, :);
        ts_est  = J_est(est_max_idx, :);
        r = corrcoef(ts_true, ts_est);
        t_corr(trial) = r(1,2);
    end
    
    % Metrics
    ED1_trend(snr_idx) = median(t_ed1);
    ED2_trend(snr_idx) = median(t_ed2);
    Temp_Corr_trend(snr_idx) = median(t_corr);
    
    err_neg_ED1(snr_idx) = ED1_trend(snr_idx) - prctile(t_ed1, 25);
    err_pos_ED1(snr_idx) = prctile(t_ed1, 75) - ED1_trend(snr_idx);
    err_neg_ED2(snr_idx) = ED2_trend(snr_idx) - prctile(t_ed2, 25);
    err_pos_ED2(snr_idx) = prctile(t_ed2, 75) - ED2_trend(snr_idx);
    err_neg_Corr(snr_idx) = Temp_Corr_trend(snr_idx) - prctile(t_corr, 25);
    err_pos_Corr(snr_idx) = prctile(t_corr, 75) - Temp_Corr_trend(snr_idx);
    
    if current_SNR == 10
        fig_3d = figure;
        fig_3d.WindowState = 'maximized';
        
        main_title_str = sprintf('EEG Source Localization Comparison | Method: sLORETA | SNR = %d dB | Median ED1 = %.1f mm', ...
            current_SNR, ED1_trend(snr_idx));
        annotation('textbox', [0, 0.95, 1, 0.05], 'String', main_title_str, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 16, 'FontWeight', 'bold');
            
        annotation('textbox', [0, 0.88, 0.5, 0.05], 'String', 'TVB Ground Truth Activity', ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 14, 'FontWeight', 'bold');
            
        annotation('textbox', [0.5, 0.88, 0.5, 0.05], 'String', 'sLORETA Estimate (Single Trial)', ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        J_true_peak = abs(J_true_mapped(:, peak_time_idx));
        cmax = 1.0; 
        
        ax1 = nexttile;
        plot_brain_surface(c_v, c_f, VN, J_true_peak, '', [-90, 0], cmax, activation_center_pos);
        
        ax2 = nexttile;
        plot_brain_surface(c_v, c_f, VN, J_est_peak, '', [-90, 0], cmax, pos_sources(:, est_max_idx));
        
        cb = colorbar(ax2);
        cb.Layout.Tile = 'east'; 
        cb.Ticks = 0:0.2:1.0;
        ylabel(cb, 'Normalized Amplitude', 'FontSize', 12);
        
        linkaxes([ax1, ax2]); 
        linkprop([ax1, ax2], {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'});

        tmp_filename = sprintf('tvb_validation_SNR_%02d.pdf', current_SNR);
        exportgraphics(fig_3d, fullfile(results_dir, tmp_filename), 'ContentType', 'image', 'Resolution', 600);
    end
end

%% Metrics
fig_tr = figure;
fig_tr.WindowState = 'maximized';

max_ed1_lim = max(ED1_trend + err_pos_ED1) + 5; 
if max_ed1_lim < 15, max_ed1_lim = 15; end

subplot(1, 3, 1); hold on; grid on;
errorbar(SNR_levels, ED1_trend, err_neg_ED1, err_pos_ED1, 'k-o', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r', 'CapSize', 5);
set(gca, 'XDir', 'reverse'); 
ylim([0, max_ed1_lim]); % FORZA LO ZERO SULL'ASSE Y
title('Localization Error (ED1)'); xlabel('SNR (dB)'); ylabel('Median ED1 (mm)'); 

subplot(1, 3, 2); hold on; grid on;
errorbar(SNR_levels, ED2_trend, err_neg_ED2, err_pos_ED2, 'b-s', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'c', 'CapSize', 5);
set(gca, 'XDir', 'reverse'); 
ylim([0, max(ED2_trend + err_pos_ED2) + 5]); % Forza lo zero
title('Spatial Dispersion (ED2)'); xlabel('SNR (dB)'); ylabel('Median ED2 (mm)'); 

subplot(1, 3, 3); hold on; grid on;
errorbar(SNR_levels, Temp_Corr_trend, err_neg_Corr, err_pos_Corr, 'g-d', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'y', 'CapSize', 5);
set(gca, 'XDir', 'reverse'); 
title('Temporal Fidelity'); xlabel('SNR (dB)'); ylabel('Pearson Correlation'); 
ylim([0 1.1]);

exportgraphics(fig_tr, fullfile(results_dir,'tvb_validation_metrics.pdf'), 'ContentType', 'vector');