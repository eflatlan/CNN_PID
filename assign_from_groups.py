def assign_attributes_from_group(data_dict, group, key, index_particle):
    """Assign attributes from the HDF5 group to the corresponding data_dict."""
    #print(f" reading assign_attributes_from_group:  index_particle {index_particle } key {key} | FILL {group.attrs[key]}")
    data_dict[key][index_particle] = group.attrs[key]

    #assign_attributes_from_group(event_data_dict, group, key, num_tracks_in_event)

#assign_values_from_group(values_data_dict, group, key, index_particle, actual_length)

def assign_values_from_group(data_dict, group, key, index_particle, actual_length):

    if index_particle >= data_dict[key].shape[0]:
        raise ValueError(f"Index {index_particle} is out of bounds for data_dict[{key}] with shape {data_dict[key].shape}")


    # Assign the values
    #print(f"assign_values_from_group group{key} {group[key][...]}")
    try:
        data_dict[key][index_particle, :actual_length] = group[key][...]
    except Exception as e:
        print(f"An error occurred while assigning values to data_dict[{key}]: {e}")
        print(f"group[key][...].shape {group[key][...].shape}")
        print(f"index_particle {index_particle} | actual_length {actual_length}")

        raise
    #print(f"assign_values_from_group  data_dict{key} { data_dict[key][index_particle, :actual_length]}")

# assign_particle_dict(particle_dict, group, key, index_particle)


def assign_particle_dict(particle_dict, group, key, index_particle):
    #print(f" reading assign_particle_dict:  index_particle {index_particle } key {key}")
    #print(f" reading assign_particle_dict:   {group.attrs[key]}")

    particle_dict[key][index_particle] = group.attrs[key]
    #print(f"  particle_dict[key][index_particle]:   {particle_dict[key][index_particle]}")

