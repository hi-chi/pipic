use std::os::raw::{c_double, c_int};

#[no_mangle]
pub unsafe extern "C" fn remove_particles_near_border(
    i: *mut c_int,
    d: *mut c_double,
    f_data: *mut c_double,
    p_data: *mut c_double,
    np_data: *mut c_double,
    data_double: *mut c_double,
    data_int: *mut c_int,
) {
    // The unused arguments are part of the framework's fixed C ABI.
    let _ = (d, f_data, data_double);

    if i.is_null() || p_data.is_null() || np_data.is_null() || data_int.is_null() {
        return;
    }

    let i_ptr = i;
    // I[] carries cell metadata; the framework stores its particle layout in doubles.
    let i = std::slice::from_raw_parts(i_ptr as *const c_int, 15);
    let particle_subset_size = i[9] as isize;
    let number_of_attributes = i[7] as isize;
    let particle_type_index = i[8] as isize;
    let buffer_capacity = i[10] as isize;
    let iz = i[2] as isize;
    let nz = i[5] as isize;
    let n_border = *data_int as isize;

    if particle_subset_size <= 0 || buffer_capacity <= 0 {
        return;
    }

    if particle_type_index != 0 {
        return;
    }

    let near_border = iz < n_border || iz >= nz - n_border;
    if !near_border {
        return;
    }

    // Each particle occupies 8 doubles plus any user-defined attributes.
    let particle_stride = 8 + number_of_attributes;
    let particles = std::slice::from_raw_parts_mut(
        p_data,
        (particle_subset_size * particle_stride) as usize,
    );
    let new_particles = std::slice::from_raw_parts_mut(
        np_data,
        (buffer_capacity * particle_stride) as usize,
    );

    let mut added = 0usize;
    for ip in 0..particle_subset_size as usize {
        let base = ip * particle_stride as usize;
        // Copy the full record before zeroing the source weight.
        let removed_particle = particles[base..base + particle_stride as usize].to_vec();

        // Slot 6 is the particle weight.
        particles[base + 6] = 0.0;

        if added < buffer_capacity as usize {
            let new_base = added * particle_stride as usize;
            new_particles[new_base..new_base + particle_stride as usize]
                .copy_from_slice(&removed_particle);
            // Slot 7 is the type tag stored in the particle id field.
            *(new_particles.as_mut_ptr().add(new_base + 7) as *mut u64) = 1;
            added += 1;
        }
    }

    // Tell the framework how many particles were placed into NP_data.
    *i_ptr.add(11) = added as c_int;
}
