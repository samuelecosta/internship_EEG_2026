%% Load data
dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
h5_dir  = fullfile(dataset_dir, 'h5_files');
mat_dir = fullfile(dataset_dir, 'mat_files');
results_dir = fullfile(pwd,'simulations','matlab');
struct_path = fullfile(mat_dir, "Surfaces_structure.mat");
G_tvb_path = fullfile(h5_dir, "ProjectionMatrix_EEG.h5");

load(struct_path);
G_tvb = h5read(G_tvb_path, "/projection_data")';
sensors = surf_struct.sensors;

%% Pre-calc
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

alpha_opt_MNE = get_global_alpha(U, s2, G_tvb_centered);

T_MNE = G_tvb' / (A + (alpha_opt_MNE * I));

% --- SVD per WMNE ---
col_norms = vecnorm(G_tvb, 2, 1)'; 
P = spdiags(1 ./ (col_norms.^2), 0, n_sources, n_sources);
A_WMNE = G_tvb * P * G_tvb';
[U_WMNE, S_mat_WMNE, ~] = svd(A_WMNE);
s2_WMNE = diag(S_mat_WMNE);

alpha_opt_WMNE = get_global_alpha(U_WMNE, s2_WMNE, G_tvb_centered);

T_WMNE = (P * G_tvb') / (A_WMNE + (alpha_opt_WMNE * I));

% --- sLORETA ---
diag_R = sum(T_MNE .* G_tvb', 2); 

ed1_pinv = zeros(n_sources, 1);
ed1_MNE  = zeros(n_sources, 1);
ed1_sLOR = zeros(n_sources, 1);
ed1_WMNE = zeros(n_sources, 1);

%% Source Loop

for n = 1:n_sources
    M = G_tvb_centered(:, n); 
    real_source = pos_sources(:, n);
    
    % PINV
    D_pinv_vec = T_pinv * M;
    [~, idx_pinv] = max(abs(D_pinv_vec));
    ed1_pinv(n) = norm(real_source - pos_sources(:, idx_pinv)) * 1000;
    
    % MNE
    D_MNE_vec = T_MNE * M;
    [~, idx_MNE] = max(abs(D_MNE_vec));
    ed1_MNE(n) = norm(real_source - pos_sources(:, idx_MNE)) * 1000;
    
    % WMNE
    D_WMNE_target = T_WMNE * M;
    [~, idx_WMNE] = max(abs(D_WMNE_target));
    ed1_WMNE(n) = norm(real_source - pos_sources(:, idx_WMNE)) * 1000;
    
    % sLORETA
    D_sLOR_vec = D_MNE_vec ./ sqrt(diag_R);
    [~, idx_sLOR] = max(abs(D_sLOR_vec));
    ed1_sLOR(n) = norm(real_source - pos_sources(:, idx_sLOR)) * 1000;
    
    if mod(n, 1000) == 0
        fprintf('Source %d/%d\n', n, n_sources); 
    end
end

%% Error 3D visualization

method_names = {'PINV', 'MNE', 'WMNE', 'sLORETA'};
method_ed1 = {ed1_pinv, ed1_MNE, ed1_WMNE, ed1_sLOR};
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
    curr_ed1 = method_ed1{m};
    method_name = method_names{m};
    
    max_ed1_method = max(curr_ed1);
    if max_ed1_method < 1, max_ed1_method = 1; end
    
    fig = figure;
    fig.WindowState = 'maximized';
    
    annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
        sprintf('Method: %s (ED1, Mean: %.1f mm | Median: %.1f mm | Max: %.1f mm)', method_name, mean(curr_ed1),median(curr_ed1), max_ed1_method), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        
    for v = 1:4
        subplot(2, 2, v);
        
        plot_brain_surface(c_v, c_f, VN, log10(curr_ed1+1), view_names{v}, view_angles{v}, log10(max_ed1_method+1));
        
        if v == 4
            cb = colorbar;
            
            tick_candidates = [0, 1, 2, 5, 10, 20, 50, 100, 150];
            real_ticks = tick_candidates(tick_candidates < max_ed1_method);
            real_ticks = [real_ticks, round(max_ed1_method)];
            
            cb.Ticks = log10(real_ticks + 1);
            cb.TickLabels = string(real_ticks);
            ylabel(cb, 'ED1 (mm)');
        end
    end

    tmp_filename = sprintf('inverse_maxperf_%s.pdf', method_name);
    exportgraphics(gcf, fullfile(results_dir, tmp_filename), 'ContentType', 'image', 'Resolution', 600);
end
%% Single source error calc

source_target = 100; %source selection
pos_true = pos_sources(:, source_target);
M = G_tvb_centered(:, source_target);

% Calcolo i vettori D usando direttamente gli operatori globali:
D_pinv_target = T_pinv * M;
D_MNE_target  = T_MNE * M;
D_WMNE_target = T_WMNE * M;
D_sLOR_target = D_MNE_target ./ sqrt(diag_R);

[~, idx_pinv] = max(abs(D_pinv_target));
ed1_pinv = norm(pos_true - pos_sources(:, idx_pinv)) * 1000;
[~, idx_MNE] = max(abs(D_MNE_target));
ed1_MNE = norm(pos_true - pos_sources(:, idx_MNE)) * 1000;
[~, idx_WMNE] = max(abs(D_WMNE_target));
ed1_WMNE = norm(pos_true - pos_sources(:, idx_WMNE)) * 1000;
[~, idx_sLOR] = max(abs(D_sLOR_target));
ed1_sLOR = norm(pos_true - pos_sources(:, idx_sLOR)) * 1000;

%normalization
D_pinv_plot = abs(D_pinv_target) / max(abs(D_pinv_target));
D_MNE_plot  = abs(D_MNE_target) / max(abs(D_MNE_target));
D_sLOR_plot = abs(D_sLOR_target) / max(abs(D_sLOR_target));
D_WMNE_plot = abs(D_WMNE_target) / max(abs(D_WMNE_target));

%% Single source error visual

method_plots = {D_pinv_plot, D_MNE_plot, D_WMNE_plot, D_sLOR_plot};
method_names = {'PINV', 'MNE', 'WMNE', 'sLORETA'};
method_ed1 = {ed1_pinv, ed1_MNE, ed1_WMNE, ed1_sLOR};

for m = 1:4
    curr_plot = method_plots{m};
    method_name = method_names{m};
    curr_ed1 = method_ed1{m};
    
    fig = figure;
    fig.WindowState = 'maximized';
    
    annotation('textbox', [0, 0.9, 1, 0.1], 'String', ...
        sprintf('Reconstructed Activity: %s, ED1: %.1f mm', method_name, curr_ed1), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        
    plot_brain_surface(c_v, c_f, VN, curr_plot, [], [-90 0], 1, pos_true);

    cb = colorbar;
    ylabel(cb, 'Normalized Amplitude');
 
    
    filename = sprintf('inverse_single_maxperf_%s.pdf', method_name);
    exportgraphics(gcf, fullfile(results_dir, filename), 'ContentType', 'image', 'Resolution', 600);
end