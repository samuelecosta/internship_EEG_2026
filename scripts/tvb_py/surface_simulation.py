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

# Initialise a Model, Coupling, and Connectivity.
oscillator = models.Generic2dOscillator(a=np.array([1.0]), b=np.array([-1.0]))

white_matter = connectivity.Connectivity.from_file()
white_matter.speed = np.array([4.0])
white_matter_coupling = coupling.Linear(a=np.array([0.005]))

# Initialise an Integrator
heunint = integrators.HeunStochastic(dt=0.5, noise=noise.Additive(nsig=np.array([1e-3])))

# Initialise a surface
default_cortex = cortex.Cortex.from_file()
default_cortex.region_mapping_data.connectivity = white_matter
default_cortex.coupling_strength = np.array([2**-10])
default_cortex.local_connectivity = local_connectivity.LocalConnectivity.from_file()

default_cortex.configure()
surface_geometry = default_cortex.region_mapping_data.surface
surface_geometry.configure()

# Target vertex configuration
target_vertex = 100
n_vertices = surface_geometry.vertices.shape[0]
target_coords = surface_geometry.vertices[target_vertex]

# Find the first triangle that contains the target_vertex
target_triangle = np.where(surface_geometry.triangles == target_vertex)[0][0]

# Spatial Gaussian equation (spread = 15.0 mm).
eqn_s = equations.Gaussian(parameters={"amp": 1.0, "sigma": 15.0, "midpoint": 0.0, "offset": 0.0})

# Temporal profile: A 50ms spike in the middle of the simulation.
eqn_t = equations.Gaussian(parameters={"amp": 5.0, "sigma": 10.0, "midpoint": 250.0, "offset": 0.0})

# Initialize the stimulus using the correct kwargs and TVB Equation objects
# Force dtype=int for the triangle index to be safe
stimulus = patterns.StimuliSurface(
    surface=surface_geometry, 
    focal_points_triangles=np.array([target_triangle], dtype=int), 
    spatial=eqn_s, 
    temporal=eqn_t
)

# Initialise some Monitors with period in physical time
mon_source = monitors.TemporalAverage(period=1.0)
mon_eeg = monitors.EEG.from_file(
    sensors_fname='eeg_unitvector_62.txt.bz2',
    projection_fname='projection_eeg_62_surface_16k.mat',
    period=1.0
)
what_to_watch = (mon_source, mon_eeg)

# Initialise Simulator
sim = simulator.Simulator(
    model=oscillator, 
    connectivity=white_matter, 
    coupling=white_matter_coupling, 
    integrator=heunint, 
    monitors=what_to_watch, 
    surface=default_cortex, 
    stimulus=stimulus,
    simulation_length=500.0
)

sim.configure()
print("Starting TVB simulation...")

# Perform the simulation
all_sources = []
all_eeg = []
time_vector = []

for source_data, eeg_data in sim():
    if source_data is not None:
        # Extract the state variable 0 (V)
        all_sources.append(source_data[1][0, :, 0])
        time_vector.append(source_data[0])
    if eeg_data is not None:
        all_eeg.append(eeg_data[1][0, :, 0])


SOURCES = np.array(all_sources) # Shape: (time, vertices)
EEG = np.array(all_eeg)         # Shape: (time, sensors)
TIME = np.array(time_vector)

mat_filepath = os.path.join(mat_dir, "tvb_validation_data.mat")

scipy.io.savemat(mat_filepath, {
    'J_true': SOURCES.T,   # Transpose to [Vertices x Time] for MATLAB
    'M_tvb': EEG.T,        # Transpose to [Sensors x Time] for MATLAB
    'time': TIME,
    'target_vertex': target_vertex,
    'target_coords': target_coords,
    'tvb_vertices': surface_geometry.vertices
})

print(f"Data saved successfully to {mat_filepath}")

# Plot region averaged time series
plt.figure(1, figsize=(10, 4))
plt.plot(TIME, SOURCES.mean(axis=1), 'k', label='Global Mean Background', alpha=0.5)
plt.plot(TIME, SOURCES[:, target_vertex], 'r', label='Target Focal Region')
plt.title("Cortical Source Dynamics (Focal Stimulus)")
plt.xlabel("Time (ms)")
plt.ylabel("Amplitude")
plt.legend()
plt.savefig(os.path.join(results_dir, "activation_signal.png"), dpi=600, bbox_inches='tight')

# Plot EEG time series
plt.figure(2, figsize=(10, 4))
plt.plot(TIME, EEG)
plt.title("EEG signal (62 sensors)")
plt.xlabel("Time (ms)")
plt.ylabel("Amplitude")
plt.savefig(os.path.join(results_dir, "eeg_signal.png"), dpi=600, bbox_inches='tight')

plt.show()