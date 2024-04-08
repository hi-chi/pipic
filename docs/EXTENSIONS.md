# Extensions

Here we list extensions with short descriptions, references and contacts of developers, which is followed by the list of parameters (if applicable the default values are specified in parenthesis). For an extension named `<name>` an extended description and the usage example can be found in `/src/extensions/<name>/<name>.cpp` and `/examples/<name>_test.py`.

- **qed_volokitin2023** is an optimized QED event generator that performs minimal possible number of rate computations per QED event [[V.&nbsp;Volokitin et al. JCS **74**, 102170 (2023)](https://doi.org/10.1016/j.jocs.2023.102170)]<br/>
*implemented by Joel Magnusson* (joel.magnusson@physics.gu.se) *, based on the implementation of pyHiChi*
    - `electron_type`, `positron_type` and `photon_type` are the type indices of electrons, positrons and photons; these types must be declared for the simulation container using `add_particles()` and then indices can be retrieved by `get_type_index(type_name)`

</br>

- **qed_gonoskov2015** is a simple implementation of QED event generator based on rejection sampling and subcycling [[A.&nbsp;Gonoskov et al. PRE **92**, 023305 (2015)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.023305); [arXiv:1412.6426](https://arxiv.org/abs/1412.6426)]; (note that `qed_volokitin2023` is an equivivalent but faster version)<br/>
*implemented by Arkady Gonoskov* (arkady.gonoskov@physics.gu.se)
    - `electron_type`, `positron_type` and `photon_type` are the type indices of electrons, positrons and photons; these types must be declared for the simulation container using `add_particles()` and then indices can be retrieved by `get_type_index(type_name)`
    - `probability_threshold`($10^{-3}$) is an optimization threshold: if the estimated probability of an event is below the threshold, we randomly either halt the computation or accordingly boost the probability
    - `probability_subcycle`(0.1) is the target probability of an event occurrence within a substep (used to compute the substep)

</br>

- **downsampler_gonoskov2022** is an implementation of agnostic (no flattening) conservative downsampling based on the method described in [[A.&nbsp;Gonoskov, CPC **271**,108200(2022)](https://doi.org/10.1016/j.cpc.2021.108200); [arXiv:1607.03755](https://arxiv.org/abs/1607.03755)]. The method provides the possibility to dynamically downsample the ensemble of particles while introducing no flattening or any other systematic changes to any distributions at any scale (agnostic downsampling), and exactly preserving any desired number of conserved quantities. The extension retrieves each subset of particles that all give CIC contributions to each set of 2/4/8 nearby nodes, depending on dimensionality. The implementation is configured to always preserve the total weight of particles in each subset. In addition, the energy, momentum and CIC contributions are preserved (can be switched off). The extension can act on several types of particles; this is to be configured with `add_assignment()` using the same parameters as for the initialization (except `ensemble_data` is not needed). </br>
*implemented by Arkady Gonoskov* (arkady.gonoskov@physics.gu.se)
    - `ensemble_data` is the address of the ensemble data to be retrieved from the simulation container by calling `ensemble_data()`
    - `type_index` is the type index of particles to be affected by the downsampling
    - `preserve_energy`(`True`), `preserve_momentum`(`True`), `preserve_cic_weight`(`True`) are the flags that define whether the respective quantities are to be preserved
    - `cap`(15) is the threshold for the number of particles in a subset (effectively per cell), after exceeding which downsampling is applied
    - `target_ratio`(1.0) defines how many particles are removed in each occurance: once called for a given subset, downsampling is configured to reduce the number of particles to `target_ratio*cap`

</br>

- **landau_lifshitz** provides a way to account for radiation reaction according to leading terms of Landau-Lifshitz model.</br>
*implemented by Joel Magnusson* (joel.magnusson@physics.gu.se)
