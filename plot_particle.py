import yt
import numpy as np
import argparse

def read_simulation_data(filename):
    data = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), len(bodies) + 3):  # Skip 'Step' and header lines
            step_data = []
            for line in lines[i + 2:i + 2 + len(bodies)]:
                x, y, mass = map(float, line.split())
                step_data.append([x, y, mass])
            data.append(step_data)
    return np.array(data)

def plot_simulation(data, idx_start, idx_end, didx, output_prefix):
    for idx in range(idx_start, idx_end + 1, didx):
        step_data = data[idx]
        ds = yt.load_particles(step_data, length_unit="Mpc", mass_unit="Msun", n_ref=64)
        p = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y', 'particle_mass', center='c')
        p.set_background_color('particle_mass')
        p.annotate
