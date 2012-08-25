!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine read_parameter_file( NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        OCEANS,TOPOGRAPHY,ANISOTROPY,ABSORBING_CONDITIONS, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,NTSTEP_BETWEEN_OUTPUT_INFO, &
                        SIMULATION_TYPE,SAVE_FORWARD, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY,IMODEL )

  implicit none

  include "constants.h"

  integer NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,SIMULATION_TYPE, NTSTEP_BETWEEN_READ_ADJSRC
  integer NSOURCES,NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,UTM_PROJECTION_ZONE
  integer NOISE_TOMOGRAPHY
  integer IMODEL

  double precision DT,HDUR_MOVIE

  logical ATTENUATION,USE_OLSEN_ATTENUATION,OCEANS,TOPOGRAPHY,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT,USE_HIGHRES_FOR_MOVIES
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,SUPPRESS_UTM_PROJECTION

  character(len=256) LOCAL_PATH,CMTSOLUTION

! local variables
  integer ::ios,icounter,isource,idummy,nproc_eta_old,nproc_xi_old
  double precision :: hdur,minval_hdur
  character(len=256) :: dummystring
  integer, external :: err_occurred

  character(len=150) MODEL
  integer :: i,irange

  ! opens file Par_file
  call open_parameter_file()

  ! reads in parameters
  call read_value_integer(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if(err_occurred() /= 0) return
  call read_value_integer(NOISE_TOMOGRAPHY, 'solver.NOISE_TOMOGRAPHY')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if(err_occurred() /= 0) return
  call read_value_integer(UTM_PROJECTION_ZONE, 'mesher.UTM_PROJECTION_ZONE')
  if(err_occurred() /= 0) return
  call read_value_logical(SUPPRESS_UTM_PROJECTION, 'mesher.SUPPRESS_UTM_PROJECTION')
  if(err_occurred() /= 0) return
  ! total number of processors
  call read_value_integer(NPROC, 'mesher.NPROC')
  if(err_occurred() /= 0) then
    ! checks if it's using an old Par_file format
    call read_value_integer(nproc_eta_old, 'mesher.NPROC_ETA')
    if( err_occurred() /= 0 ) then
      print*,'please specify the number of processes in Par_file as:'
      print*,'NPROC           =    <my_number_of_desired_processes> '
      return
    endif
    ! checks if it's using an old Par_file format
    call read_value_integer(nproc_xi_old, 'mesher.NPROC_XI')
    if( err_occurred() /= 0 ) then
      print*,'please specify the number of processes in Par_file as:'
      print*,'NPROC           =    <my_number_of_desired_processes> '
      return
    endif
    NPROC = nproc_eta_old * nproc_xi_old
  endif
  call read_value_integer(NSTEP, 'solver.NSTEP')
  if(err_occurred() /= 0) return
  call read_value_double_precision(DT, 'solver.DT')
  if(err_occurred() /= 0) return

  ! define the velocity model
  call read_value_string(MODEL, 'model.MODEL')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MODEL'

  call read_value_logical(OCEANS, 'model.OCEANS')
  if(err_occurred() /= 0) return
  call read_value_logical(TOPOGRAPHY, 'model.TOPOGRAPHY')
  if(err_occurred() /= 0) return
  call read_value_logical(ATTENUATION, 'model.ATTENUATION')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_OLSEN_ATTENUATION, 'model.USE_OLSEN_ATTENUATION')
  if(err_occurred() /= 0) return
  call read_value_logical(ANISOTROPY, 'model.ANISOTROPY')
  if(err_occurred() /= 0) return
  call read_value_logical(ABSORBING_CONDITIONS, 'solver.ABSORBING_CONDITIONS')
  if(err_occurred() /= 0) return
  call read_value_logical(MOVIE_SURFACE, 'solver.MOVIE_SURFACE')
  if(err_occurred() /= 0) return
  call read_value_logical(MOVIE_VOLUME, 'solver.MOVIE_VOLUME')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'solver.NTSTEP_BETWEEN_FRAMES')
  if(err_occurred() /= 0) return
  call read_value_logical(CREATE_SHAKEMAP, 'solver.CREATE_SHAKEMAP')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_DISPLACEMENT, 'solver.SAVE_DISPLACEMENT')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_HIGHRES_FOR_MOVIES, 'solver.USE_HIGHRES_FOR_MOVIES')
  if(err_occurred() /= 0) return
  call read_value_double_precision(HDUR_MOVIE, 'solver.HDUR_MOVIE')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_MESH_FILES, 'mesher.SAVE_MESH_FILES')
  if(err_occurred() /= 0) return
  call read_value_string(LOCAL_PATH, 'LOCAL_PATH')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'solver.NTSTEP_BETWEEN_OUTPUT_INFO')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'solver.NTSTEP_BETWEEN_OUTPUT_SEISMOS')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_READ_ADJSRC, 'solver.NTSTEP_BETWEEN_READ_ADJSRC')
  if(err_occurred() /= 0) return
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'solver.PRINT_SOURCE_TIME_FUNCTION')
  if(err_occurred() /= 0) return

  ! close parameter file
  call close_parameter_file()

  ! noise simulations:
  ! double the number of time steps, if running noise simulations (+/- branches)
  if ( NOISE_TOMOGRAPHY /= 0 )   NSTEP = 2*NSTEP-1

  ! the default value of NTSTEP_BETWEEN_READ_ADJSRC (0) is to read the whole trace at the same time
  if ( NTSTEP_BETWEEN_READ_ADJSRC == 0 )  NTSTEP_BETWEEN_READ_ADJSRC = NSTEP

  ! total times steps must be dividable by adjoint source chunks/blocks
  if ( mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) /= 0 ) then
    print*,'error: mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) must be zero!'
    print*,'      change your Par_file (when NOISE_TOMOGRAPHY is not equal to zero, ACTUAL_NSTEP=2*NSTEP-1)'
    stop 'mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) must be zero!'
  endif

  ! for noise simulations, we need to save movies at the surface (where the noise is generated)
  ! and thus we force MOVIE_SURFACE to be .true., in order to use variables defined for surface movies later
  if ( NOISE_TOMOGRAPHY /= 0 ) then
    MOVIE_SURFACE = .true.
    CREATE_SHAKEMAP = .false.           ! CREATE_SHAKEMAP and MOVIE_SURFACE cannot be both .true.
    USE_HIGHRES_FOR_MOVIES = .true.     ! we need to save surface movie everywhere, i.e. at all GLL points on the surface
    ! since there are several flags involving surface movies, check compatability
    if ( EXTERNAL_MESH_MOVIE_SURFACE .or. EXTERNAL_MESH_CREATE_SHAKEMAP ) then
        print*, 'error: when running noise simulations ( NOISE_TOMOGRAPHY /= 0 ),'
        print*, '       we can NOT use EXTERNAL_MESH_MOVIE_SURFACE or EXTERNAL_MESH_CREATE_SHAKEMAP'
        print*, '       change EXTERNAL_MESH_MOVIE_SURFACE & EXTERNAL_MESH_CREATE_SHAKEMAP in constant.h'
        stop 'incompatible NOISE_TOMOGRAPHY, EXTERNAL_MESH_MOVIE_SURFACE, EXTERNAL_MESH_CREATE_SHAKEMAP'
    endif
  endif

  ! compute the total number of sources in the CMTSOLUTION file
  ! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file
  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION',&
       IN_DATA_FILES_PATH(1:len_trim(IN_DATA_FILES_PATH))//'CMTSOLUTION')

  open(unit=21,file=trim(CMTSOLUTION),iostat=ios,status='old',action='read')
  if(ios /= 0) stop 'error opening CMTSOLUTION file'

  icounter = 0
  do while(ios == 0)
    read(21,"(a)",iostat=ios) dummystring
    if(ios == 0) icounter = icounter + 1
  enddo
  close(21)

  if(mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
    stop 'total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'

  NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE
  if(NSOURCES < 1) stop 'need at least one source in CMTSOLUTION file'

  ! compute the minimum value of hdur in CMTSOLUTION file
  open(unit=21,file=trim(CMTSOLUTION),status='old',action='read')
  minval_hdur = HUGEVAL
  do isource = 1,NSOURCES

    ! skip other information
    do idummy = 1,3
      read(21,"(a)") dummystring
    enddo

    ! read half duration and compute minimum
    read(21,"(a)") dummystring
    read(dummystring(15:len_trim(dummystring)),*) hdur
    minval_hdur = min(minval_hdur,hdur)

    ! skip other information
    do idummy = 1,9
      read(21,"(a)") dummystring
    enddo

  enddo
  close(21)

  ! one cannot use a Heaviside source for the movies
  if((MOVIE_SURFACE .or. MOVIE_VOLUME) .and. sqrt(minval_hdur**2 + HDUR_MOVIE**2) < TINYVAL) &
    stop 'hdur too small for movie creation, movies do not make sense for Heaviside source'

  ! converts all string characters to lowercase
  irange = iachar('a') - iachar('A')
  do i = 1,len_trim(MODEL)
    if( lge(MODEL(i:i),'A') .and. lle(MODEL(i:i),'Z') ) then
      MODEL(i:i) = achar( iachar(MODEL(i:i)) + irange )
    endif
  enddo

  ! determines velocity model
  select case( trim(MODEL) )

  ! default mesh model
  case( 'default' )
    IMODEL = IMODEL_DEFAULT

  ! 1-D models
  case( '1d_prem' )
    IMODEL = IMODEL_1D_PREM
  case( '1d_socal' )
    IMODEL = IMODEL_1D_SOCAL
  case( '1d_cascadia')
    IMODEL = IMODEL_1D_CASCADIA

  ! user models
  case( '1d_prem_pb' )
    IMODEL = IMODEL_1D_PREM_PB
  case( 'aniso' )
    IMODEL = IMODEL_DEFAULT
    ANISOTROPY = .true.
  case( 'external' )
    IMODEL = IMODEL_USER_EXTERNAL
  case( 'ipati' )
    IMODEL = IMODEL_IPATI
  case( 'gll' )
    IMODEL = IMODEL_GLL
  case( 'salton_trough')
    IMODEL = IMODEL_SALTON_TROUGH
  case( 'tomo' )
    IMODEL = IMODEL_TOMO

  case default
    print*
    print*,'********** model not recognized: ',trim(MODEL),' **************'
    print*,'********** using model: default',' **************'
    print*
    IMODEL = IMODEL_DEFAULT
  end select

  ! check
  if( IMODEL == IMODEL_IPATI ) then
    if( USE_RICKER_IPATI .eqv. .false. ) stop 'please set USE_RICKER_IPATI to true in shared/constants.h and recompile'
  endif

  end subroutine read_parameter_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_gpu_mode(GPU_MODE,GRAVITY)

  implicit none
  include "constants.h"

  logical :: GPU_MODE
  logical :: GRAVITY

  ! initializes flags
  GPU_MODE = .false.
  GRAVITY = .false.

  ! opens file Par_file
  call open_parameter_file()

  call read_value_logical(GPU_MODE, 'solver.GPU_MODE')
  call read_value_logical(GRAVITY, 'solver.GRAVITY')

  ! close parameter file
  call close_parameter_file()

  end subroutine read_gpu_mode
