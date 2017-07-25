Module arrays

  use healpix_types

  Integer :: status1,status2,status3,status4,status5,status6
  Integer*8, allocatable, dimension(:) :: jackknife_data_indices

  Real*8, allocatable, dimension(:) :: redshift,luminosity_distance,galactic_longitude,galactic_latitude
  Real*8, allocatable, dimension(:) :: jackknife_dipole_amplitude,jackknife_galactic_longitude,jackknife_galactic_latitude
  Real(kind=DP), allocatable, dimension(:,:) :: map, map_nc, shot_noise, mask
  Real(kind=DP), allocatable, dimension(:,:) :: cl, nl

End module arrays
