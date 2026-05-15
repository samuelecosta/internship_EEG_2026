%% Load data
dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
h5_dir  = fullfile(dataset_dir, 'h5_files');
mat_dir = fullfile(dataset_dir, 'mat_files');
results_dir = fullfile(pwd,'simulations','matlab');
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

%% Pre-calc

SNR_dB = 5;

n_sources = size(surf_struct.cortical_surface.vertices, 2);
pos_sources = surf_struct.cortical_surface.vertices;

G_tvb_centered = G_tvb - mean(G_tvb, 1);

I = eye(size(G_sym_sens, 1));

% --- Inverse per pinv ---
T_pinv = pinv(full(G_sym_sens));

% --- SVD per MNE ---
A = G_sym_sens * G_sym_sens';
[U, S_mat, ~] = svd(A); 
s2 = diag(S_mat);

%alpha_opt_MNE = get_global_alpha(U, s2, G_tvb_centered);
%alpha_opt_MNE = 0;
tr_A = trace(A);
alpha_opt_MNE = (tr_A / n_sensors) * (10^(-SNR_dB / 10));

T_MNE = G_sym_sens' / (A + (alpha_opt_MNE * I));

% --- SVD per WMNE ---
col_norms = vecnorm(G_sym_sens, 2, 1)'; 
P = spdiags(1 ./ (col_norms.^2), 0, n_sources, n_sources);
A_WMNE = G_sym_sens * P * G_sym_sens';
[U_WMNE, S_mat_WMNE, ~] = svd(A_WMNE);
s2_WMNE = diag(S_mat_WMNE);

%alpha_opt_WMNE = get_global_alpha(U_WMNE, s2_WMNE, G_tvb_centered);
%alpha_opt_WMNE = 0;
tr_A_WMNE = trace(A_WMNE);
alpha_opt_WMNE = (tr_A_WMNE / n_sensors) * (10^(-SNR_dB / 10));

T_WMNE = (P * G_sym_sens') / (A_WMNE + (alpha_opt_WMNE * I));

% --- sLORETA ---
diag_R = sum(T_MNE .* G_sym_sens', 2); 

ed1_pinv = zeros(n_sources, 1);
ed1_MNE  = zeros(n_sources, 1);
ed1_sLOR = zeros(n_sources, 1);
ed1_WMNE = zeros(n_sources, 1);

ed2_pinv = zeros(n_sources, 1);
ed2_MNE  = zeros(n_sources, 1);
ed2_sLOR = zeros(n_sources, 1);
ed2_WMNE = zeros(n_sources, 1);

%% Source Loop

for n = 1:n_sources
    M_clean = G_tvb_centered(:, n); 

    signal_power = var(M_clean);
    noise_power = signal_power / (10^(SNR_dB/10));

    noise = sqrt(noise_power) * randn(size(M_clean));

    M = M_clean + noise;

    real_source = pos_sources(:, n);

    dist_all = vecnorm(pos_sources - real_source, 2, 1)';
    
    % PINV
    D_pinv_vec = T_pinv * M;
    [~, idx_pinv] = max(abs(D_pinv_vec));
    ed1_pinv(n) = dist_all(idx_pinv);
    ed2_pinv(n) = sum(dist_all .* (abs(D_pinv_vec) / max(abs(D_pinv_vec))));
    
    % MNE
    D_MNE_vec = T_MNE * M;
    [~, idx_MNE] = max(abs(D_MNE_vec));
    ed1_MNE(n) = dist_all(idx_MNE);
    ed2_MNE(n) = sum(dist_all .* (abs(D_MNE_vec) / max(abs(D_MNE_vec))));
    
    % WMNE
    D_WMNE_vec = T_WMNE * M;
    [~, idx_WMNE] = max(abs(D_WMNE_vec));
    ed1_WMNE(n) = dist_all(idx_WMNE);
    ed2_WMNE(n) = sum(dist_all .* (abs(D_WMNE_vec) / max(abs(D_WMNE_vec))));
    
    % sLORETA
    D_sLOR_vec = D_MNE_vec ./ sqrt(diag_R);
    [~, idx_sLOR] = max(abs(D_sLOR_vec));
    ed1_sLOR(n) = dist_all(idx_sLOR);
    ed2_sLOR(n) = sum(dist_all .* (abs(D_sLOR_vec) / max(abs(D_sLOR_vec))));
    
    if mod(n, 1000) == 0
        fprintf('Source %d/%d\n', n, n_sources); 
    end
end

%% Error 3D visualization

method_names = {'PINV', 'MNE', 'WMNE', 'sLORETA'};
method_ed1 = {ed1_pinv, ed1_MNE, ed1_WMNE, ed1_sLOR};
method_ed2 = {ed2_pinv, ed2_MNE, ed2_WMNE, ed2_sLOR};
view_angles = {[-90, 0], [90, 0], [0, 90], [0, -90]};
view_names = {'Rear', 'Front', 'Superior (Top)', 'Inferior (Bottom)'};

c_v = double(surf_struct.cortical_surface.vertices);
c_f = double(surf_struct.cortical_surface.triangles' + 1);
if size(c_v,2) ~= 3, c_v = c_v.'; end
if size(c_f,2) ~= 3, c_f = c_f.'; end

% lighting normals
TR = triangulation(c_f, c_v);
VN = -vertexNormal(TR); 

for m = 1:4
    curr_ed1 = method_ed1{m} * 1000;
    method_name = method_names{m};
    
    max_ed1_method = max(curr_ed1);
    if max_ed1_method < 1, max_ed1_method = 1; end
    
    fig = figure;
    fig.WindowState = 'maximized';
    
    annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
        sprintf('Method: %s (ED1, Mean: %.1f mm | Median: %.1f mm | Max: %.1f mm), SNR = %.1f dB', method_name, mean(curr_ed1),median(curr_ed1), max_ed1_method, SNR_dB), ...
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

    tmp_filename = sprintf('inverse_ed1_%s_%d_dB.pdf', method_name, SNR_dB);
    exportgraphics(gcf, fullfile(results_dir, tmp_filename), 'ContentType', 'image', 'Resolution', 600);
end

% for m = 1:4
%     curr_ed2 = method_ed2{m};
%     method_name = method_names{m};
% 
%     max_ed2_method = max(curr_ed2);
%     if max_ed2_method < 1, max_ed2_method = 1; end
% 
%     fig = figure;
%     fig.WindowState = 'maximized';
% 
%     annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
%         sprintf('Method: %s (ED2, Mean: %.1f mm | Median: %.1f mm | Max: %.1f mm)', method_name, mean(curr_ed2),median(curr_ed2), max_ed2_method), ...
%         'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
% 
%     for v = 1:4
%         subplot(2, 2, v);
% 
%         plot_brain_surface(c_v, c_f, VN, log10(curr_ed2+1), view_names{v}, view_angles{v}, log10(max_ed2_method+1));
% 
%         if v == 4
%             cb = colorbar;
% 
%             tick_candidates = [0, 1, 2, 5, 10, 20, 50, 100, 150];
%             real_ticks = tick_candidates(tick_candidates < max_ed2_method);
%             real_ticks = [real_ticks, round(max_ed2_method)];
% 
%             cb.Ticks = log10(real_ticks + 1);
%             cb.TickLabels = string(real_ticks);
%             ylabel(cb, 'ED2 (mm)');
%         end
%     end
% 
%     tmp_filename = sprintf('inverse_ed2_%s.pdf', method_name);
%     %exportgraphics(gcf, fullfile(results_dir, tmp_filename), 'ContentType', 'image', 'Resolution', 600);
% end

%% Max performance calc

n_sources = size(surf_struct.cortical_surface.vertices, 2);
pos_sources = surf_struct.cortical_surface.vertices;

G_tvb_centered = G_tvb - mean(G_tvb, 1);

I = eye(size(G_tvb, 1));

% --- Inverse per pinv ---
T_pinv = pinv(full(G_tvb));

% --- SVD per MNE ---
A = G_tvb * G_tvb';
[U, S_mat, ~] = svd(A); 
s2 = diag(S_mat);

%alpha_opt_MNE = get_global_alpha(U, s2, G_tvb_centered);
%alpha_opt_MNE = 0;
tr_A = trace(A);
alpha_opt_MNE = (tr_A / n_sensors) * (10^(-SNR_dB / 10));

T_MNE = G_tvb' / (A + (alpha_opt_MNE * I));

% --- SVD per WMNE ---
col_norms = vecnorm(G_tvb, 2, 1)'; 
P = spdiags(1 ./ (col_norms.^2), 0, n_sources, n_sources);
A_WMNE = G_tvb * P * G_tvb';
[U_WMNE, S_mat_WMNE, ~] = svd(A_WMNE);
s2_WMNE = diag(S_mat_WMNE);

%alpha_opt_WMNE = get_global_alpha(U_WMNE, s2_WMNE, G_tvb_centered);
%alpha_opt_WMNE = 0;
tr_A_WMNE = trace(A_WMNE);
alpha_opt_WMNE = (tr_A_WMNE / n_sensors) * (10^(-SNR_dB / 10));

T_WMNE = (P * G_tvb') / (A_WMNE + (alpha_opt_WMNE * I));

% --- sLORETA ---
diag_R = sum(T_MNE .* G_tvb', 2); 

ed1_pinv = zeros(n_sources, 1);
ed1_MNE  = zeros(n_sources, 1);
ed1_sLOR = zeros(n_sources, 1);
ed1_WMNE = zeros(n_sources, 1);

ed2_pinv = zeros(n_sources, 1);
ed2_MNE  = zeros(n_sources, 1);
ed2_sLOR = zeros(n_sources, 1);
ed2_WMNE = zeros(n_sources, 1);

%% Source Loop

for n = 1:n_sources
    M_clean = G_tvb_centered(:, n); 

    signal_power = var(M_clean);
    noise_power = signal_power / (10^(SNR_dB/10));

    noise = sqrt(noise_power) * randn(size(M_clean));

    M = M_clean + noise;

    real_source = pos_sources(:, n);

    dist_all = vecnorm(pos_sources - real_source, 2, 1)';
    
    % PINV
    D_pinv_vec = T_pinv * M;
    [~, idx_pinv] = max(abs(D_pinv_vec));
    ed1_pinv(n) = dist_all(idx_pinv);
    ed2_pinv(n) = sum(dist_all .* (abs(D_pinv_vec) / max(abs(D_pinv_vec))));
    
    % MNE
    D_MNE_vec = T_MNE * M;
    [~, idx_MNE] = max(abs(D_MNE_vec));
    ed1_MNE(n) = dist_all(idx_MNE);
    ed2_MNE(n) = sum(dist_all .* (abs(D_MNE_vec) / max(abs(D_MNE_vec))));
    
    % WMNE
    D_WMNE_vec = T_WMNE * M;
    [~, idx_WMNE] = max(abs(D_WMNE_vec));
    ed1_WMNE(n) = dist_all(idx_WMNE);
    ed2_WMNE(n) = sum(dist_all .* (abs(D_WMNE_vec) / max(abs(D_WMNE_vec))));
    
    % sLORETA
    D_sLOR_vec = D_MNE_vec ./ sqrt(diag_R);
    [~, idx_sLOR] = max(abs(D_sLOR_vec));
    ed1_sLOR(n) = dist_all(idx_sLOR);
    ed2_sLOR(n) = sum(dist_all .* (abs(D_sLOR_vec) / max(abs(D_sLOR_vec))));
    
    if mod(n, 1000) == 0
        fprintf('Source %d/%d\n', n, n_sources); 
    end
end

%% Error 3D visualization

method_names = {'PINV', 'MNE', 'WMNE', 'sLORETA'};
method_ed1 = {ed1_pinv, ed1_MNE, ed1_WMNE, ed1_sLOR};
method_ed2 = {ed2_pinv, ed2_MNE, ed2_WMNE, ed2_sLOR};
view_angles = {[-90, 0], [90, 0], [0, 90], [0, -90]};
view_names = {'Rear', 'Front', 'Superior (Top)', 'Inferior (Bottom)'};

c_v = double(surf_struct.cortical_surface.vertices);
c_f = double(surf_struct.cortical_surface.triangles' + 1);
if size(c_v,2) ~= 3, c_v = c_v.'; end
if size(c_f,2) ~= 3, c_f = c_f.'; end

% lighting normals
TR = triangulation(c_f, c_v);
VN = -vertexNormal(TR); 

for m = 1:4
    curr_ed1 = method_ed1{m} * 1000;
    method_name = method_names{m};
    
    max_ed1_method = max(curr_ed1);
    if max_ed1_method < 1, max_ed1_method = 1; end
    
    fig = figure;
    fig.WindowState = 'maximized';
    
    annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
        sprintf('Max performance, method: %s (ED1, Mean: %.1f mm | Median: %.1f mm | Max: %.1f mm), SNR = %.1f dB', method_name, mean(curr_ed1),median(curr_ed1), max_ed1_method, SNR_dB), ...
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

    tmp_filename = sprintf('inverse_ed1_maxperf_%s_%d_dB.pdf', method_name, SNR_dB);
    exportgraphics(gcf, fullfile(results_dir, tmp_filename), 'ContentType', 'image', 'Resolution', 600);
end

% for m = 1:4
%     curr_ed2 = method_ed2{m};
%     method_name = method_names{m};
% 
%     max_ed2_method = max(curr_ed2);
%     if max_ed2_method < 1, max_ed2_method = 1; end
% 
%     fig = figure;
%     fig.WindowState = 'maximized';
% 
%     annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
%         sprintf('Method: %s (ED2, Mean: %.1f mm | Median: %.1f mm | Max: %.1f mm)', method_name, mean(curr_ed2),median(curr_ed2), max_ed2_method), ...
%         'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
% 
%     for v = 1:4
%         subplot(2, 2, v);
% 
%         plot_brain_surface(c_v, c_f, VN, log10(curr_ed2+1), view_names{v}, view_angles{v}, log10(max_ed2_method+1));
% 
%         if v == 4
%             cb = colorbar;
% 
%             tick_candidates = [0, 1, 2, 5, 10, 20, 50, 100, 150];
%             real_ticks = tick_candidates(tick_candidates < max_ed2_method);
%             real_ticks = [real_ticks, round(max_ed2_method)];
% 
%             cb.Ticks = log10(real_ticks + 1);
%             cb.TickLabels = string(real_ticks);
%             ylabel(cb, 'ED2 (mm)');
%         end
%     end
% 
%     tmp_filename = sprintf('inverse_ed2_%s.pdf', method_name);
%     %exportgraphics(gcf, fullfile(results_dir, tmp_filename), 'ContentType', 'image', 'Resolution', 600);
% end