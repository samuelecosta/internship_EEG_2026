import os
import numpy as np
import scipy.io
from tvb.simulator.lab import *
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
py_dir = os.path.join(project_root, "datasets", "tvb_default", "py_generated")
results_dir = os.path.join(project_root, "simulations", "python")

#Initialise a Model, Coupling, and Connectivity.
oscillator = models.Generic2dOscillator()

white_matter = connectivity.Connectivity.from_file()
white_matter.speed = np.array([4.0])
white_matter_coupling = coupling.Linear(a = np.array([0.014]))

#Initialise an Integrator
heunint = integrators.HeunDeterministic(dt = 0.1)

#Initialise a surface
default_cortex = cortex.Cortex.from_file()
default_cortex.region_mapping_data.connectivity = white_matter
default_cortex.coupling_strength = np.array([2**-10])
default_cortex.local_connectivity = local_connectivity.LocalConnectivity.from_file()

#Initialise some Monitors with period in physical time
mon_source = monitors.TemporalAverage(period=1.0)
# load the default region mapping
mon_eeg = monitors.EEG.from_file(
    sensors_fname='eeg_unitvector_62.txt.bz2',
    projection_fname='projection_eeg_62_surface_16k.mat',
    period=1.0
)
#Bundle them
what_to_watch = (mon_source, mon_eeg)

#Initialise Simulator -- Model, Connectivity, Integrator, Monitors, and surface.
sim = simulator.Simulator(model = oscillator, connectivity = white_matter, coupling = white_matter_coupling, integrator = heunint, monitors = what_to_watch, surface = default_cortex, simulation_length=500.0)

sim.configure()

#Perform the simulation
all_sources = []
all_eeg = []
time_vector = []

for source_data, eeg_data in sim():
    if source_data is not None:
        # Estraiamo la variabile 0 (V) e salviamo
        all_sources.append(source_data[1][0, :, 0])
        time_vector.append(source_data[0])
    if eeg_data is not None:
        all_eeg.append(eeg_data[1][0, :, 0])

#Make the lists numpy.arrays for easier use.
SOURCES = np.array(all_sources)
EEG = np.array(all_eeg)
TIME = np.array(time_vector)

mat_filepath = os.path.join(py_dir, "tvb_validation_data.mat")

scipy.io.savemat(mat_filepath, {
    'J_true': SOURCES.T,
    'M_tvb': EEG.T,
    'time': TIME
})

#Plot region averaged time series
plt.figure(1)
plt.plot(TIME, SOURCES[:, :55])
plt.title("Sources activation (first 55)")
plt.xlabel("Time (ms)")
plt.ylabel("Amplitude")

plt.savefig(os.path.join(results_dir, "activation_signal.png"), dpi=600, bbox_inches='tight')

#Plot EEG time series
plt.figure(2)
plt.plot(TIME, EEG)
plt.title("EEG signal (63 sensors)")
plt.xlabel("Time (ms)")
plt.ylabel("Amplitude")

plt.savefig(os.path.join(results_dir, "eeg_signal.png"), dpi=600, bbox_inches='tight')

#Show them
plt.show()