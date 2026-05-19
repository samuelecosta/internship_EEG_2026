import os
import numpy as np
import scipy.io
from tvb.simulator.lab import *
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
mat_dir = os.path.join(project_root, "datasets", "tvb_default", "mat_files")
results_dir = os.path.join(project_root, "simulations", "python")

if not os.path.exists(results_dir):
    os.makedirs(results_dir)

white_matter = connectivity.Connectivity.from_file()
white_matter.speed = np.array([4.0])
white_matter_coupling = coupling.Difference(a=np.array([1e-3])) 

noise_vector = np.array([1e-6, 1e-6, 0.0, 1e-6, 1e-6, 0.0])
hiss = noise.Additive(nsig=noise_vector)

heunint = integrators.HeunStochastic(dt=0.1, noise=hiss)

default_cortex = cortex.Cortex.from_file()
default_cortex.region_mapping_data.connectivity = white_matter
default_cortex.coupling_strength = np.array([1e-2]) 
default_cortex.local_connectivity = local_connectivity.LocalConnectivity.from_file()

default_cortex.configure()
surface_geometry = default_cortex.region_mapping_data.surface
surface_geometry.configure()

target_vertex = 499
n_vertices = surface_geometry.vertices.shape[0]
target_coords = surface_geometry.vertices[target_vertex]

actual_simulation_length = 1000.0  
burn_in_length = 500.0            
total_length = actual_simulation_length + burn_in_length

x0_spatial = np.full(n_vertices, -2.4)
patch_radius_mm = 5.0
distances = np.linalg.norm(surface_geometry.vertices - target_coords, axis=1)

ez_indices = np.where(distances <= patch_radius_mm)[0]
x0_spatial[ez_indices] = -1.6  

epileptor = models.Epileptor(x0=x0_spatial)

mon_source = monitors.TemporalAverage(period=2.0)
mon_eeg = monitors.EEG.from_file(
    sensors_fname='eeg_unitvector_62.txt.bz2',
    projection_fname='projection_eeg_62_surface_16k.mat',
    period=2.0
)
what_to_watch = (mon_source, mon_eeg)

sim = simulator.Simulator(
    model=epileptor, 
    connectivity=white_matter, 
    coupling=white_matter_coupling, 
    integrator=heunint, 
    monitors=what_to_watch, 
    surface=default_cortex, 
    simulation_length=total_length
)

sim.configure()
print("TVB Epileptor simulation...")

all_sources = []
all_eeg = []
time_vector = []

for source_data, eeg_data in sim():
    if source_data is not None:
        all_sources.append(source_data[1][0, :, 0])
        time_vector.append(source_data[0])
    if eeg_data is not None:
        all_eeg.append(eeg_data[1][0, :, 0])

SOURCES_RAW = np.array(all_sources) 
EEG_RAW = np.array(all_eeg)         
TIME_RAW = np.array(time_vector)

burn_in_idx = np.argmax(TIME_RAW >= burn_in_length)

SOURCES = SOURCES_RAW[burn_in_idx:, :]
EEG = EEG_RAW[burn_in_idx:, :]
TIME = TIME_RAW[burn_in_idx:] - TIME_RAW[burn_in_idx]

SOURCES = SOURCES - np.mean(SOURCES, axis=0)
EEG = EEG - np.mean(EEG, axis=0)

mat_filepath = os.path.join(mat_dir, "tvb_epilepsy_data.mat")

scipy.io.savemat(mat_filepath, {
    'J_true': SOURCES.T,   
    'M_tvb': EEG.T,        
    'time': TIME,
    'target_vertex': target_vertex,
    'target_coords': target_coords,
    'tvb_vertices': surface_geometry.vertices,
    'activation_center': target_vertex
})

print(f"Data saved successfully to {mat_filepath}")

plt.figure(1, figsize=(10, 4))
plt.plot(TIME, SOURCES.mean(axis=1), 'k', label='Global Mean Background', alpha=0.5)
plt.plot(TIME, SOURCES[:, target_vertex], 'r', label='Epileptogenic Zone (Fast Discharges)')
plt.title("Cortical Seizure Dynamics")
plt.xlabel("Time (ms)")
plt.ylabel("Amplitude")
plt.legend()
plt.savefig(os.path.join(results_dir, "epilepsy_activation_signal.png"), dpi=600, bbox_inches='tight')

plt.figure(2, figsize=(10, 4))
plt.plot(TIME, EEG)
plt.title("EEG signal during Seizure (62 sensors)")
plt.xlabel("Time (ms)")
plt.ylabel("Amplitude")
plt.savefig(os.path.join(results_dir, "epilepsy_eeg_signal.png"), dpi=600, bbox_inches='tight')

plt.show()