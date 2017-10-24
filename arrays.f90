Module arrays

  use healpix_types

  Integer :: status1,status2,status3,status4,status5,status6
  Integer*8, allocatable, dimension(:) :: jackknife_data_indices

  Real*8, allocatable, dimension(:) :: redshift,luminosity_distance,galactic_longitude,galactic_latitude,redshift_corrected
  Real*8, allocatable, dimension(:) :: redshift_jla,galactic_longitude_jla,galactic_latitude_jla
  Real*8, allocatable, dimension(:) :: jackknife_dipole_amplitude,jackknife_galactic_longitude,jackknife_galactic_latitude
  Real*8, allocatable, dimension(:) :: jack_noise_dipole_amplitude,jack_noise_galactic_longitude,jack_noise_galactic_latitude
  Real(kind=DP), allocatable, dimension(:,:) :: map, map_nc, shot_noise, mask
  Real(kind=DP), allocatable, dimension(:,:) :: cl, nl

End module arrays
