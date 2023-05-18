!pip install h5py numpy

import h5py
import numpy as np

class ParticleDataUtils:
    class ParticleInfo:
        def __init__(self, momentum, mass, energy, refractiveIndex, ckov, map_data, filledBins):
            self.momentum = momentum
            self.mass = mass
            self.energy = energy
            self.refractiveIndex = refractiveIndex
            self.ckov = ckov
            self.map = map_data
            self.filledBins = filledBins

def save_particle_info_to_hdf5(particle_vector, filename):
    with h5py.File(filename, 'w') as file:
        for i, particle in enumerate(particle_vector):
            # Create a group for each particle
            group = file.create_group(f'Particle{i}')

            # Store scalar values
            group.attrs['Momentum'] = particle.momentum
            group.attrs['Mass'] = particle.mass
            group.attrs['Energy'] = particle.energy
            group.attrs['RefractiveIndex'] = particle.refractiveIndex
            group.attrs['Ckov'] = particle.ckov

            # Store 2D array "map"
            group.create_dataset("Map", data=np.array(particle.map))

            # Write filledBins to HDF5 file
            group.create_dataset("FilledBins", data=np.array(particle.filledBins, dtype=[('x', 'f'), ('y', 'f')]))

def load_particle_info_from_hdf5(filename):
    particle_vector = []
    
    with h5py.File(filename, 'r') as file:
        for i, group_name in enumerate(file):
            group = file[group_name]

            # Read scalar values
            momentum = group.attrs['Momentum']
            mass = group.attrs['Mass']
            energy = group.attrs['Energy']
            refractiveIndex = group.attrs['RefractiveIndex']
            ckov = group.attrs['Ckov']

            # Read the map
            map_dataset = group['Map']
            map_data = np.array(map_dataset)

            # Read filledBins
            filled_bins_dataset = group['FilledBins']
            filled_bins_data = np.array(filled_bins_dataset)

            particle_info = ParticleDataUtils.ParticleInfo(
                momentum, mass, energy, refractiveIndex, ckov, map_data, filled_bins_data)
            
            particle_vector.append(particle_info)

    return particle_vector

