clear all; close all; clc;

%% Load data
dataset_dir = fullfile(pwd, 'datasets', 'tvb_default');
h5_dir  = fullfile(dataset_dir, 'h5_files');
mat_dir = fullfile(dataset_dir, 'mat_files');
py_dir = fullfile(dataset_dir, 'py_generated');
results_dir = fullfile(pwd,'simulations','matlab');

struct_path = fullfile(mat_dir, "Surfaces_structure.mat");
Three_layer_G_path = fullfile(mat_dir, "G_three_layers.mat");
G_tvb_path = fullfile(h5_dir, "ProjectionMatrix_EEG.h5");
val_data_path = fullfile(py_dir, "tvb_validation_data.mat");

load(struct_path);
load(Three_layer_G_path, 'G_three_layers'); 
load(val_data_path);
G_tvb = h5read(G_tvb_path, "/projection_data")';
sensors = surf_struct.sensors;

%% visual

figure
plot(time, M_tvb)

figure
plot(time, J_true)