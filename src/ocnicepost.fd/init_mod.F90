module init_mod

  implicit none

  public

  real, parameter :: maxvars = 70              !< The maximum number of fields expected in a source file

  type :: valid_ranges
     real, allocatable :: rng2d(:)
     real, allocatable :: rng3d(:,:)
  end type valid_ranges

  type :: vardefs
     type(valid_ranges)   :: ranges            !< Nested valid_ranges for each vardefs
     character(len= 20)   :: var_name          !< A variable's variable name
     character(len= 20)   :: var_remapmethod   !< A variable's mapping method
     character(len=120)   :: long_name         !< A variable's long name
     character(len= 20)   :: units             !< A variable's unit
     integer              :: var_dimen         !< A variable's dimensionality
     character(len=  4)   :: var_grid          !< A variable's input grid location; all output locations are on cell centers
     character(len= 20)   :: var_pair          !< A variable's pair
     character(len=  4)   :: var_pair_grid     !< A pair variable grid
     logical              :: isvector          !< A logical indicating a vector field
     real                 :: var_fillvalue     !< A variable's fillvalue
     character(len= 20)   :: name_gb2          !< A variable's grib2 variable name
     character(len=120)   :: discription_gb2   !< A variable's discription
     character(len= 20)   :: unit_gb2          !< A variable's unit
!!!may need to add a fillvalue for grib2 file
     integer              ::var_g1             !< Variables' grib2 coefficients g1-Dissipline
     integer              ::var_g2             !< Variables' grib2 coefficients g2-Master Tables Version Number
     integer              ::var_g3             !< Variables' grib2 coefficients g3-Section 1 originating center, used for local tables
     integer              ::var_g4             !< Variables' grib2 coefficients g4-Section 1 Local Tables Version Number
     integer              ::var_g5             !< Variables' grib2 coefficients g5-Section 4 Template 4.0 Parameter category
     integer              ::var_g6             !< Variables' grib2 coefficients g6-Section 4 Template 4.0 Parameter number
     integer              ::var_g7             !< Variables' grib2 coefficients g7-Level ID
     integer              ::var_g8             !< Variables' grib2 coefficients g8-
  end type vardefs

  type(vardefs) :: outvars(maxvars)            !< An empty structure filled by reading a csv file describing the fields

  character(len=10)  :: ftype      !< The type of tripole grid (ocean or ice)
  character(len=10)  :: fsrc       !< A character string for tripole grid
  character(len=10)  :: fdst       !< A character string for the destination grid
  character(len=120) :: wgtsdir    !< The directory containing the regridding weights
  character(len=20)  :: input_file !< The input file name
  character(len=10)  :: maskvar    !< The variable in the source file used to create the interpolation mask

  ! rotation angles
  character(len=10)  :: sinvar    !< The variable in the source file containing the cosine of the rotation angle (ocean only)
  character(len=10)  :: cosvar    !< The variable in the source file containing the sine of the rotation angle (ocean only)
  character(len=10)  :: angvar    !< The variable in the source file containing the rotation angle (ice only)

  integer :: nxt        !< The x-dimension of the source tripole grid
  integer :: nyt        !< The y-dimension of the source tripole grid
  integer :: nlevs      !< The vertical dimension of the source tripole grid

  integer :: nxr        !< The x-dimension of the destination rectilinear grid
  integer :: nyr        !< The y-dimension of the destination rectilinear grid

  integer :: logunit          !< The log unit
  logical :: write_grib2      !< If true, write grib2 message
  logical :: write_netcdf     !< If true, write netCDF flies
  logical :: debug            !< If true, print debug messages and intermediate files
  logical :: do_ocnpost       !< If true, the source file is ocean, otherwise ice

  integer :: nvpairb2d         !< The number of 2d vector pairs, bilinear
  integer :: nvpairc2d         !< The number of 2d vector pairs, conservative
  integer :: nvpairb3d         !< The number of 2d vector pairs, bilinear

contains

  subroutine readnml

    ! local variable
    character(len=40) :: fname
    integer :: ierr, iounit
    integer :: srcdims(2), dstdims(2)

    namelist /ocnicepost_nml/ ftype, srcdims, wgtsdir, dstdims, maskvar, sinvar, cosvar, &
         angvar, write_grib2, write_netcdf, debug

    ! --------------------------------------------------------
    ! read the name list
    ! --------------------------------------------------------

    srcdims = 0; dstdims = 0
    sinvar=''; cosvar=''; angvar=''

    fname = 'ocnicepost.nml'
    inquire (file=trim(fname), iostat=ierr)
    if (ierr /= 0) then
       write (0, '(3a)') 'FATAL ERROR: input file "', trim(fname), '" does not exist.'
       stop 1
    end if

    ! Open and read namelist file.
    open (action='read', file=trim(fname), iostat=ierr, newunit=iounit)
    read (nml=ocnicepost_nml, iostat=ierr, unit=iounit)
    if (ierr /= 0) then
       write (6, '(a)') 'Error: invalid namelist format.'
    end if
    close (iounit)
    nxt = srcdims(1); nyt = srcdims(2)
    nxr = dstdims(1); nyr = dstdims(2)

    ! initialize the source file type and variables
    if (trim(ftype) == 'ocean') then
       do_ocnpost = .true.
    else
       do_ocnpost = .false.
    end if
    input_file = trim(ftype)//'.nc'

    open(newunit=logunit, file=trim(ftype)//'.post.log',form='formatted')
    if (debug) write(logunit, '(a)')'input file: '//trim(input_file)

    ! set grid names
    fsrc = ''
    if (nxt == 1440 .and. nyt == 1080) fsrc = 'mx025'    ! 1/4deg tripole
    if (nxt == 720  .and. nyt == 576) fsrc = 'mx050'     ! 1/2deg tripole
    if (nxt == 360  .and. nyt == 320) fsrc = 'mx100'     ! 1deg tripole
    if (nxt == 72   .and. nyt == 35) fsrc = 'mx500'      ! 5deg tripole
    if (len_trim(fsrc) == 0) then
       write(0,'(a)')'FATAL ERROR: source grid dimensions unknown'
       stop 2
    end if

    fdst = ''
    if (nxr == 1440 .and. nyr == 721) fdst = '0p25'      ! 1/4deg rectilinear
    if (nxr == 720  .and. nyr == 361) fdst = '0p50'      ! 1/2 deg rectilinear
    if (nxr == 360  .and. nyr == 181) fdst = '1p00'      ! 1 deg rectilinear
    if (nxr == 72   .and. nyr == 36) fdst = '5p00'       ! 5 deg rectilinear
    if (len_trim(fdst) == 0) then
       write(0,'(a)')'FATAL ERROR: destination grid dimensions unknown'
       stop 3
    end if

    !TODO: test for consistency of source/destination resolution
  end subroutine readnml

  subroutine readcsv(nvalid)

    integer, intent(out) :: nvalid

    character(len= 40) :: fname
    character(len=100) :: chead
    character(len= 20) :: c1,c3,c4,c5,c6,c7,c8,c9
    integer :: i2,i10,i11,i12,i13,i14,i15,i16,i17
    integer :: nn,n,ierr,iounit

    ! --------------------------------------------------------
    ! Open and read list of variables
    ! --------------------------------------------------------

    fname=trim(ftype)//'.csv'
    open(newunit=iounit, file=trim(fname), status='old', iostat=ierr)
    if (ierr /= 0) then
       write (0, '(3a)') 'FATAL ERROR: input file "', trim(fname), '" does not exist.'
       stop 4
    end if

    ! initialize
    outvars%var_name = ''
    outvars%var_dimen = 0
    outvars%var_grid = ''
    outvars%var_remapmethod = ''
    outvars%var_pair = ''
    outvars%var_pair_grid = ''

    read(iounit,*)chead
    nn=0
    do n = 1,maxvars
       read(iounit,*,iostat=ierr)c1,i2,c3,c4,c5,c6,c7,c8,c9,i10,i11,i12,i13,i14,i15,i16,i17
       if (ierr .ne. 0) exit
       if (len_trim(c1) > 0) then
          nn = nn+1
          outvars(nn)%var_name = trim(c1)
          outvars(nn)%var_dimen = i2
          outvars(nn)%var_grid = trim(c3)
          outvars(nn)%var_remapmethod = trim(c4)
          outvars(nn)%var_pair = trim(c5)
          outvars(nn)%var_pair_grid = trim(c6)
          outvars(nn)%name_gb2 = trim(c7)
          outvars(nn)%discription_gb2 = trim(c8)
          outvars(nn)%unit_gb2 = trim(c9)
          outvars(nn)%var_g1 = i10
          outvars(nn)%var_g2 = i11
          outvars(nn)%var_g3 = i12
          outvars(nn)%var_g4 = i13
          outvars(nn)%var_g5 = i14
          outvars(nn)%var_g6 = i15
          outvars(nn)%var_g7 = i16
          outvars(nn)%var_g8 = i17
       end if
    end do
    close(iounit)
    nvalid = nn
    outvars%isvector = .false.
    where (len_trim(outvars%var_pair) > 0) outvars%isvector = .true.
    do n = 1,maxvars-1
       if (outvars(n)%isvector)then
          print *,'XX '//trim(outvars(n)%var_name)//'  '//trim(outvars(n+1)%var_name),outvars(n+1)%isvector

          !if (trim(outvars(n)%var_pair) .ne. trim(outvars(n+1)%var_name)) then
          !   write(0,'(a,i0)') 'FATAL ERROR: vector pairs must occur sequentially in list'
          !   stop 6
          !end if
       end if
    end do

    do n = 1,maxvars,2
       if (outvars(n)%isvector) then
       !   print *,'XX '//trim(outvars(n)%var_name)//'  '//trim(outvars(n+1)%var_name),outvars(n+1)%isvector
       end if
       !if (outvars(n)%isvector) .and. .not. outvars(n+1)%isvector) then
       !   write(0,'(a,i0)') 'FATAL ERROR: vector pairs must occur sequentially in list'
       !   stop 6
       !end if
    end do

    nvpairb2d = countvpairs(outvars, 'bilinear', 2)
    nvpairc2d = countvpairs(outvars, 'conserve', 2)
    nvpairb3d = countvpairs(outvars, 'bilinear', 3)

    !seqpairs = checkvpairs(outvars)

    !print *,'XX ',nvpairb2d, nvpairc2d, nvpairb3d
  end subroutine readcsv

  integer function countvpairs(vars, maptype, size) result(count)

    type(vardefs),    intent(in) :: vars(:)
    character(len=*), intent(in) :: maptype
    integer,          intent(in) :: size

    integer :: n

    count = 0
    do n = 1,maxvars
       if (outvars(n)%isvector) then
          if (vars(n)%var_dimen == size .and. trim(vars(n)%var_remapmethod)  == trim(maptype)) count = count+1
       end if
    end do

    if (mod(count,2) /= 0) then
       write(0,'(a,i0)') 'FATAL ERROR: non-even number of vector pairs map = '//trim(maptype)//' dimension = ',size
       stop 5
    else
       count = count / 2
    end if
  end function countvpairs

  ! logical function checkvpairs(vars) result(seqpairs)

  !   type(vardefs),    intent(in) :: vars(:)

  !   integer :: i,n

  !   do n = 1,maxvars,2
  !      if (vars(n)%isvector .and. .not.vars(n+1)%isvector) then
  !         write(0,'(a,i0)') 'FATAL ERROR: vector pairs must occur sequentially in list'
  !         stop 6
  !end function checkvpairs
end module init_mod
