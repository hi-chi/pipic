# Extensions

Here we list extensions with short descriptions, references and contacts of developers, which is followed by the list of parameters (if applicable the default values are specified in parenthesis). For an extension named `<name>` an extended description and the usage example can be found in `/src/extensions/<name>/<name>.cpp` and `/examples/<name>_test.py`.

- **qed_volokitin2023** is an optimized QED event generator that performs minimal possible number of rate computations per QED event [[V.&nbsp;Volokitin et al. JCS **74**, 102170 (2023)](https://doi.org/10.1016/j.jocs.2023.102170)]<br/>
*implemented by Joel Magnusson* (joel.magnusson@physics.gu.se) *, based on the implementation of pyHiChi*
    - `electron_type`, `positron_type` and `photon_type` are the type indices of electrons, positrons and photons; these types must be declared for the simulation container using `add_particles()` and then indices can be retrieved by `get_type_index(type_name)`

</br>

- **qed_gonoskov2015** is a simple implementation of QED event generator based on rejection sampling and subcycling [[A.&nbsp;Gonoskov et al. PRE **92**, 023305 (2015)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.023305); [arXiv:1412.6426](https://arxiv.org/abs/1412.6426)]; (note that `qed_volokitin2023` is an equivalent but faster version)<br/>
*implemented by Arkady Gonoskov* (arkady.gonoskov@physics.gu.se)
    - `electron_type`, `positron_type` and `photon_type` are the type indices of electrons, positrons and photons; these types must be declared for the simulation container using `add_particles()` and then indices can be retrieved by `get_type_index(type_name)`
    - `probability_threshold`($10^{-3}$) is an optimization threshold: if the estimated probability of an event is below the threshold, we randomly either halt the computation or accordingly boost the probability
    - `probability_subcycle`(0.1) is the target probability of an event occurrence within a substep (used to compute the substep)

</br>

- **downsampler_gonoskov2022** is an implementation of agnostic (no flattening) conservative downsampling based on the method described in [[A.&nbsp;Gonoskov, CPC **271**,108200 (2022)](https://doi.org/10.1016/j.cpc.2021.108200); [arXiv:1607.03755](https://arxiv.org/abs/1607.03755)]. The method provides the possibility to dynamically downsample the ensemble of particles while introducing no flattening or any other systematic changes to any distributions at any scale (agnostic downsampling), and exactly preserving any desired number of conserved quantities. The extension retrieves each subset of particles that all give CIC contributions to each set of 2/4/8 nearby nodes, depending on dimensionality. The implementation is configured to always preserve the total weight of particles in each subset. In addition, the energy, momentum and CIC contributions are preserved (can be switched off). The extension can act on several types of particles; this is to be configured with `add_assignment()` using the same parameters as for the initialization (except `ensemble_data` is not needed). </br>
*implemented by Arkady Gonoskov* (arkady.gonoskov@physics.gu.se)
    - `ensemble_data` is the address of the ensemble data to be retrieved from the simulation container by calling `ensemble_data()`
    - `type_index` is the type index of particles to be affected by the downsampling
    - `preserve_energy`(`True`), `preserve_momentum`(`True`), `preserve_cic_weight`(`True`) are the flags that define whether the respective quantities are to be preserved
    - `cap`(15) is the threshold for the number of particles in a subset (effectively per cell), after exceeding which downsampling is applied
    - `target_ratio`(1.0) defines how many particles are removed in each occurrence: once called for a given subset, downsampling is configured to reduce the number of particles to `target_ratio*cap`

</br>

- **landau_lifshitz** provides a way to account for radiation reaction according to leading terms of Landau-Lifshitz model.</br>
*implemented by Joel Magnusson* (joel.magnusson@physics.gu.se)

- **focused_pulse** is an extension for precomputing and setting the initial state of a tightly focused (up to $4\pi$) electromagnetic pulse. The extension takes specific structure of the pulse at a distance from the focus and computes the resultant electromagnetic field at the instance when the pulse gets into the computation box. The method is based on "folding" initial space into simulation region under assumption of space periodicity and then advancing the field state by a single step using Fourier field solver (the idea follows [[E.&nbsp;Panova et al. Appl. Sci., 11(3), 956 (2021)]](https://www.mdpi.com/2076-3417/11/3/956) with a further development of permitting pulse overlaying in the initial instance of time).<br/>
*implemented by Christoffer Olofsson* (christoffer.olofsson@physics.gu.se) *and Arkady Gonoskov* </br>
    - `set_box(nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax)` should be used to set the computational region of the extension
    - `set_path(x, y, z)` is to set the path to go (defines the axis)
    - `set_e_axis(e_axis_x, e_axis_y, e_axis_z)` is to set the polarization axis (electric field)
    - `set_l_size(l_size)` is to set the longitudinal size of the pulse
    - `set_shape(pulse_shape.address)` is to set the shape of the pulse using a callback (can be used to set, e.g. amplitude, polarization and phased variations across the pulse):
      ```
      @cfunc(numba.types.double(numba.types.CPointer(numba.types.double)))
      def pulse_shape(par): # initial field amplitude (at point r) as a function of the following parameters
        # par[0] = dist - abs(r), longitudinal coordinate relative to dist (\in [-l_size/2, l_size/2])
        # par[1] is the angle between r and -path, theta angle in polar coordinates (\in [0, pi))
        # par[2] is the angle between e_axis and r in the plane perpendicular to path (\in [0, 2*pi))
        # par[3] is the angle between r and e_axis (\in [0, pi))
      return "any function of par"
      ```
      The use of the extension is examplified with the case of setting electric dipole wave in 'examples/focused_pulse_test.py'.
