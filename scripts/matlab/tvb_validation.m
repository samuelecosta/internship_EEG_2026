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

tvb_data = load(fullfile(mat_dir, "tvb_epilepsy_data.mat")); %change this to select the data
J_true_raw = double(tvb_data.J_true); 
time = double(tvb_data.time);
n_time = length(time);

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

W_src_map = zeros(n_sources, size(tvb_vertices, 2));
for i = 1:n_sources
    dists = vecnorm(tvb_vertices - pos_sources(:, i), 2, 1);
    [sorted_dists, idx] = sort(dists);
    k_nearest = 3; 
    
    if sorted_dists(1) < 1e-6
        W_src_map(i, idx(1)) = 1;
    else
        weights = 1 ./ sorted_dists(1:k_nearest);
        W_src_map(i, idx(1:k_nearest)) = weights / sum(weights);
    end
end

J_true_mapped = W_src_map * J_true_raw;

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

true_target_power = J_true_mapped(target_idx, :).^2;
power_weights = true_target_power / sum(true_target_power);
[~, absolute_peak_t] = max(true_target_power);

threshold = 0.25; %for time window for metrics
active_mask = (true_target_power > (threshold * max(true_target_power)));

SNR_levels = [25, 15, 10, 5, 0]; 
n_snr = length(SNR_levels);

n_trials = 50; %monte carlo

ED1_trend = zeros(n_snr, 1); err_neg_ED1 = zeros(n_snr, 1); err_pos_ED1 = zeros(n_snr, 1);
ED2_trend = zeros(n_snr, 1); err_neg_ED2 = zeros(n_snr, 1); err_pos_ED2 = zeros(n_snr, 1);
Env_Corr_trend = zeros(n_snr, 1); err_neg_Corr = zeros(n_snr, 1); err_pos_Corr = zeros(n_snr, 1);
Spat_Corr_trend = zeros(n_snr, 1); err_neg_Spat = zeros(n_snr, 1); err_pos_Spat = zeros(n_snr, 1);

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
    
    t_ed1 = zeros(n_trials, 1); 
    t_ed2 = zeros(n_trials, 1); 
    t_env_corr = zeros(n_trials, 1);
    t_spat_corr = zeros(n_trials, 1);
    
    for trial = 1:n_trials
        noise_power = var(M_clean(:)) / (10^(current_SNR/10));
        noise = sqrt(noise_power) * randn(size(M_clean));
        M_noisy = M_clean + noise;
        
        % Inverse sLORETA
        J_est_MNE = T_MNE * M_noisy;
        J_est_energy = J_est_MNE.^2;
        J_est_stat = J_est_energy ./ diag_R; 
        
        % 1. ED1 Istant
        [~, max_verts_all_t] = max(J_est_stat, [], 1);
        ed1_inst = dist_from_target(max_verts_all_t)' * 1000;
        
        % 2. ED2 Istant
        energy_sum_all_t = sum(J_est_stat, 1);
        energy_sum_all_t(energy_sum_all_t == 0) = eps;
        ed2_inst = (dist_from_target' * J_est_stat) ./ energy_sum_all_t * 1000; 
        
        t_ed1(trial) = sum(ed1_inst(active_mask) .* power_weights(active_mask)) / sum(power_weights(active_mask));
        t_ed2(trial) = sum(ed2_inst(active_mask) .* power_weights(active_mask)) / sum(power_weights(active_mask));
        
        % 3. Envelope Temporal Correlation
        [~, best_overall_vert] = max(J_est_stat * power_weights');
        env_true = abs(hilbert(J_true_mapped(target_idx, :)'))'; 
        env_est  = abs(hilbert(J_est_MNE(best_overall_vert, :)'))';
        r_env = corrcoef(env_true, env_est);
        t_env_corr(trial) = r_env(1,2);
        
        % 4. Spatial Correlation at peak
        r_spat = corrcoef(J_true_mapped(:, absolute_peak_t).^2, J_est_stat(:, absolute_peak_t));
        t_spat_corr(trial) = r_spat(1,2);
    end
    
    % Metrics
    ED1_trend(snr_idx) = median(t_ed1);
    err_neg_ED1(snr_idx) = ED1_trend(snr_idx) - prctile(t_ed1, 25);
    err_pos_ED1(snr_idx) = prctile(t_ed1, 75) - ED1_trend(snr_idx);
    
    ED2_trend(snr_idx) = median(t_ed2);
    err_neg_ED2(snr_idx) = ED2_trend(snr_idx) - prctile(t_ed2, 25);
    err_pos_ED2(snr_idx) = prctile(t_ed2, 75) - ED2_trend(snr_idx);
    
    Env_Corr_trend(snr_idx) = median(t_env_corr);
    err_neg_Corr(snr_idx) = Env_Corr_trend(snr_idx) - prctile(t_env_corr, 25);
    err_pos_Corr(snr_idx) = prctile(t_env_corr, 75) - Env_Corr_trend(snr_idx);
    
    Spat_Corr_trend(snr_idx) = median(t_spat_corr);
    err_neg_Spat(snr_idx) = Spat_Corr_trend(snr_idx) - prctile(t_spat_corr, 25);
    err_pos_Spat(snr_idx) = prctile(t_spat_corr, 75) - Spat_Corr_trend(snr_idx);
    
    if current_SNR == 10

        fig_3d = figure; fig_3d.WindowState = 'maximized';
        
        main_title_str = sprintf('Local peak (T = %.1f ms) | sLORETA | SNR = %d dB', ...
            time(absolute_peak_t), current_SNR);
        annotation('textbox', [0, 0.95, 1, 0.05], 'String', main_title_str, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
            
        J_true_snap = J_true_mapped(:, absolute_peak_t).^2;
        J_true_snap = J_true_snap / max(J_true_snap); 
        ax1 = nexttile;
        plot_brain_surface(c_v, c_f, VN, J_true_snap, 'Ground Truth (Peak Power)', [-75, 0], 1.0, activation_center_pos);
        
        [~, snap_est_max] = max(J_est_stat(:, absolute_peak_t));
        J_est_snap = J_est_stat(:, absolute_peak_t);
        J_est_snap = J_est_snap / max(J_est_snap);
        ax2 = nexttile;
        plot_brain_surface(c_v, c_f, VN, J_est_snap, 'sLORETA Estimate', [-75, 0], 1.0, pos_sources(:, snap_est_max));
        
        cb = colorbar(ax2); cb.Layout.Tile = 'east'; ylabel(cb, 'Normalized Power', 'FontSize', 12);
        linkaxes([ax1, ax2]); linkprop([ax1, ax2], {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'});
        exportgraphics(fig_3d, fullfile(results_dir, sprintf('tvb_validation_epilepsy_SNR_%02d_v2.pdf', current_SNR)), 'ContentType', 'image', 'Resolution', 600);
    end
end

%% Metrics
fig_tr = figure; fig_tr.WindowState = 'maximized';
max_ed1_lim = max(ED1_trend + err_pos_ED1) + 5; 
if max_ed1_lim < 15, max_ed1_lim = 15; end

% ED1 (Weighted)
subplot(2, 2, 1); hold on; grid on;
errorbar(SNR_levels, ED1_trend, err_neg_ED1, err_pos_ED1, 'k-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'CapSize', 5);
set(gca, 'XDir', 'reverse'); ylim([0, max_ed1_lim]); 
title('Energy-Weighted Localization Error (ED1)'); xlabel('SNR (dB)'); ylabel('Median Weighted ED1 (mm)'); 

% ED2 (Weighted)
subplot(2, 2, 2); hold on; grid on;
errorbar(SNR_levels, ED2_trend, err_neg_ED2, err_pos_ED2, 'b-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'c', 'CapSize', 5);
set(gca, 'XDir', 'reverse'); ylim([0, max(ED2_trend + err_pos_ED2) + 5]); 
title('Energy-Weighted Spatial Dispersion (ED2)'); xlabel('SNR (dB)'); ylabel('Median Weighted ED2 (mm)'); 

% Envelope Temporal Correlation
subplot(2, 2, 3); hold on; grid on;
errorbar(SNR_levels, Env_Corr_trend, err_neg_Corr, err_pos_Corr, 'g-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'CapSize', 5);
set(gca, 'XDir', 'reverse'); ylim([0 1.1]);
title('Envelope Temporal Correlation'); xlabel('SNR (dB)'); ylabel('Pearson r (Hilbert Envelope)'); 

% Peak Spatial Correlation
subplot(2, 2, 4); hold on; grid on;
errorbar(SNR_levels, Spat_Corr_trend, err_neg_Spat, err_pos_Spat, 'm-^', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'CapSize', 5);
set(gca, 'XDir', 'reverse'); ylim([0 1.1]);
title('Peak-Frame Spatial Correlation'); xlabel('SNR (dB)'); ylabel('Pearson r'); 

exportgraphics(fig_tr, fullfile(results_dir,'tvb_validation_epilepsy_metrics_v2.pdf'), 'ContentType', 'vector');