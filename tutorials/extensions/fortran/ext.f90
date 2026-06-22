module particle_removal_extension
	use iso_c_binding, only: c_double, c_int
	implicit none
contains
	!
	! C-compatible handler for pi-PIC cell callbacks.
	!
	! The handler removes all particles in cells that are within n cells of any
	! simulation boundary, where n is provided through dataInt(1).
	!
	! The project stores each particle as a flat block of 8 + numberOfAttributes
	! doubles in P_data. The particle weight is the seventh double in each block,
	! so setting that entry to zero marks the particle for removal.
	!
	! Expected C signature:
	!   void handler(int *I, double *D, double *F_data, double *P_data,
	!                double *NP_data, double *dataDouble, int *dataInt)
	!
	subroutine remove_particles_near_border(I, D, F_data, P_data, NP_data, dataDouble, dataInt) &
		bind(C, name="remove_particles_near_border")
		integer(c_int), intent(in)    :: I(*)
		real(c_double), intent(inout) :: D(*)
		real(c_double), intent(inout) :: F_data(*)
		real(c_double), intent(inout) :: P_data(*)
		real(c_double), intent(inout) :: NP_data(*)
		real(c_double), intent(inout) :: dataDouble(*)
		integer(c_int), intent(in)    :: dataInt(*)

		integer(c_int) :: particle_subset_size
		integer(c_int) :: number_of_attributes
		integer(c_int) :: particle_stride
		integer(c_int) :: ix, iy, iz
		integer(c_int) :: nx, ny, nz
		integer(c_int) :: n_border
		integer(c_int) :: ip
		integer(c_int) :: base
		logical :: near_border

		particle_subset_size = I(10)
		number_of_attributes = I(8)
		particle_stride = 8 + number_of_attributes
		ix = I(1)
		iy = I(2)
		iz = I(3)
		nx = I(4)
		ny = I(5)
		nz = I(6)
		n_border = dataInt(1)

		near_border = (iz < n_border) .or. (iz >= nz - n_border)

		if (near_border) then
			do ip = 0, particle_subset_size - 1
				base = ip * particle_stride
				P_data(base + 7) = 0.0_c_double
			end do
		end if
	end subroutine remove_particles_near_border
end module particle_removal_extension
