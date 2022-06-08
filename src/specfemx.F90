! This is a main program SPECFEMX
! REVISION:
!  HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
! NOTE:
!  - Shear strain components in the strain tensor are engineering strain.
!  - Ordinary shear strain components can be obtained by multiplying the
!    engineering strain with 0.5 if necessary.
program specfemx
! Import necessary modules.
use dimensionless
use global
use package_version
use string_library, only : parse_file
use input
use mesh_spec
use element
use infinite_element,only:classify_finite_infinite_elements
use dof
use integration
use model
use gll_library,only:precompute_gll1d,cleanup_gll1d
use free_surface
#if (USE_MPI)
use mpi_library
use math_library_mpi
use output_to_user
#else
use serial_library
use math_library_serial
#endif
use nondimensionalisation
use save_mesh
use postprocess,only:write_scalar_to_file_freesurf
use write_ensight

implicit none
integer :: iounit,iounit_inf,iounit_fs,i,ios,j,k
integer :: i_elmt

integer :: gnum_hex8(8),node_hex8(8)
integer :: gnum_quad4(4),node_quad4(4)

character(len=250) :: arg1,arg2,inp_fname,prog
character(len=150) :: path
character(len=80) :: buffer ! this must be 80 characters long
character(len=20) :: ext,format_str
character(len=250) :: case_file,geo_file
character(len=250) :: infcase_file,infgeo_file,trinfcase_file,trinfgeo_file
character(len=250) :: fscase_file,fsgeo_file
character(len=250) :: fspcase_file,fspgeo_file
! switch to check if the geometry file changes with time steps, for example,
! multistage excavation
logical :: isgeo_change
integer :: npart,twidth
! Ensight Gold "Time Section" variables
! ts: time set
! ns: number of steps
! fs: filename start number
! fi: filename increment
integer :: ns,fi,fs,ts ! ts: time set for ensight gold
integer,allocatable :: ipart(:)
character(len=80) :: spart_fs(1) ! this must be 80 characters long
character(len=80),allocatable :: spart(:) ! this must be 80 characters long
real(kind=kreal) :: absmaxx,absmaxy,absmaxz
real(kind=kreal) :: cpu_tstart,cpu_tend,telap,max_telap,mean_telap

integer :: tot_nelmt,max_nelmt,min_nelmt,tot_nnode,max_nnode,min_nnode

character(len=250) :: cmd ! command line
character(len=8) :: tdate ! date
character(len=10) :: ttime ! time
character(len=5) :: tzone ! time zone

character(len=250) :: errtag ! error message
integer :: errcode

character(len=60) :: add_tag
! flag to check whether the file is opened
logical :: isopen
! if the following flag is true program will stop after saving the mesh files
logical :: ismesh_only
myrank=0; nproc=1;
errtag=""; errcode=-1



call start_process()


call get_command_argument(0, prog)
if (command_argument_count() <= 0) then
  errcode=-1
  errtag='ERROR: no input file!'
  call control_error(errcode,errtag,stdout,myrank)
endif

call get_command_argument(1, arg1)
if(trim(arg1)==('--help'))then
  if(myrank==0)then
    write(stdout,'(a)')'Usage: '//trim(prog)//' [Options] [input_file]'
    write(stdout,'(a)')'Options:'
    write(stdout,'(a)')'    --help        : Display this information.'
    write(stdout,'(a)')'    --version     : Display version information.'
  endif
  !call sync_process
  call close_process()
elseif(trim(arg1)==('--version'))then
  if(myrank==0)then
    write(stdout,'(a)')trim(packname)//' '//trim(packver)//' '//trim(packtype)
    write(stdout,'(a)')'This is free software; see the source for copying '
    write(stdout,'(a)')'conditions.  There is NO warranty; not even for '
    write(stdout,'(a)')'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
  endif
  !call sync_process
  call close_process()
endif
ismesh_only=.false.
call get_command_argument(2, arg2)
if(trim(arg2)==('--mesh-only'))then
  ismesh_only=.true.
endif

! get input file name
call get_command_argument(1, inp_fname)

! parse input file name for future use
call parse_file(inp_fname,path,file_head,ext)

! open log file
! log file must be opened after parsing the input file name since it uses
! file_head!
! log file must be stored in the current directory since out_path isn't known
! at this stage!
if(myrank==0)then
  log_file = trim(file_head)//'.log'
  open(unit=logunit,file=trim(log_file),status='replace',action='write',iostat=ios)
  if(ios.ne.0)then
    print*,ios,trim(log_file)
    write(errtag,'(a)')'ERROR: cannot open log file: '//trim(log_file)
    call control_error(errcode,errtag,stdout,myrank)
  endif
  write(logunit,'(a)')'--------------------------------------------'
  write(logunit,'(a)')'Result summary produced by '//&
  trim(packname)//' '//trim(packver)//' '//trim(packtype)
  write(logunit,'(a)')'--------------------------------------------'
  write(logunit,'(a)')'DATE(CCYYMMDD) TIME(HHMMSS.SSS) ZONE(+-HHMM)'
  call date_and_time(tdate,ttime,tzone)
  write(logunit,'(3x,a,7x,a,7x,a)')tdate,ttime,tzone
  call get_command(cmd)
  write(logunit,'(a)')trim(cmd)
  flush(logunit)
endif

! starting timer
call cpu_time(cpu_tstart)

! get processor tag
ptail=proc_tag()
if(ismpi.and.nproc.gt.1)then
  ptail_inp=trim(ptail)
else
  ptail_inp=''
endif


proc_str=''
if(ismpi.and.nproc.gt.1)then
  write(format_str,*)ceiling(log10(real(nproc)+1.))
  format_str='(i'//trim(adjustl(format_str))//')'
  write(proc_str,fmt=format_str)myrank
endif


! read input data
call read_input(inp_fname,errcode,errtag)
call sync_process()
call control_error(errcode,errtag,stdout,myrank)


! check method
if (trim(method)/='sem')then
  write(errtag,'(a)')'ERROR: wrong input for SPECFEM3D!'
  call control_error(errcode,errtag,stdout,myrank)
else
  write(logunit, '(a)')'Correct method: sem'
endif


! Calculate model extents for individual processors/whole model
call calc_model_coord_extents(tot_nelmt,max_nelmt,min_nelmt, &
                              tot_nnode,max_nnode,min_nnode, &
                              absmaxx,absmaxy,absmaxz)



! Initialize model - allocates shearmod/bulkmod/massdensity arrays
call initialize_model(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)



isgeo_change=.false.
ts=1 ! time set
fs=0; fi=1
!if(nexcav==0)then
!  nt=nsrf
!else
!  nt=nexcav+1 ! include 0 excavation stage (i.e., initial)
!  tstart=0
!  isgeo_change=.true.
!endif
! nt for ensight file format must be > 0
ns=max(1,nstep)
twidth=ceiling(log10(real(ns)+1.))

! NEED TO EDIT LOCAL VARIABLES STILL 
call write_original_mesh


! If only saving mesh (not running simulation) then quit program at this point. 
if(ismesh_only)then
  call close_process
endif

! ___________________________________________________________________________________


! store orginal connectivity which helps to identify ghost interfaces
allocate(g_num0(ngnode,nelmt))
g_num0=g_num


! precompute gll 1D
call precompute_gll1d()


! create spectral elements
log_msg = trim('creating spectral elements...') ;   call write_ifproc0()

call hex2spec(ndim,ngnode,nelmt,nnode,ngllx,nglly,ngllz,errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)

log_msg = trim('completed creating spectral elements') ;   call write_ifproc0()


tot_nelmt=sumscal(nelmt); tot_nnode=sumscal(nnode)
max_nelmt=maxscal(nelmt); max_nnode=maxscal(nnode)
min_nelmt=minscal(nelmt); min_nnode=minscal(nnode)
if(myrank==0)then
  write(logunit,'(a,i0,1x,a,i0,1x,a,i0)')' spectral elements => total:',tot_nelmt, &
  ' max:',max_nelmt,' min:',min_nelmt
  write(logunit,'(a,i0,1x,a,i0,1x,a,i0)')' spectral nodes    => total:',tot_nnode, &
  ' max:',max_nnode,' min:',min_nnode
  flush(logunit)
endif


! number of elemental nodes (nodes per element)
nenode=ngll !(ngllx*nglly*ngllz)

! Reclassify (in)finite elements if there is an infinite boundary condition
if(infbc)then
  ! free memory for modified arrays
  deallocate(g_num_finite,g_num_trinfinite,g_num_infinite, &
             node_finite,node_trinfinite,node_infinite)
  call classify_finite_infinite_elements(1)
endif

! initialize and set element DOFs
call initialize_dof()
call set_element_dof()


! Write model details to log file for user 
call print_model_details()

! prepare hexes
call prepare_hex(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)

! prepare hex faces
call prepare_hexface(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)

! prepare integration
call prepare_integration(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)

! prepare surface (2D) integration
call prepare_integration2d(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)

! Set model properties
call set_model_properties(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)

! Nondimensionalisation 
call set_nondimensional_params
call calc_nondimensionalisation_vals



! Read and prepare free surface file.
! Information is later used to determine the elevation of the source point and
! to plot the free surface files.
call prepare_free_surface(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)
case_file=trim(out_path)//trim(file_head)//trim(ptail)//'.case'
if(nexcav==0)then
  geo_file=trim(file_head)//trim(ptail)//'.geo'
else
  isgeo_change=.true.
  geo_file=trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.geo'
endif

add_tag=''
call write_ensight_casefile_long(case_file,geo_file,add_tag,isgeo_change, &
ts,ns,fs,fi,twidth,errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)




call save_mesh_ensight(infcase_file,infgeo_file,trinfcase_file, &
                       trinfgeo_file,isgeo_change,add_tag, twidth, &
                       fscase_file,fsgeo_file, fspcase_file,fspgeo_file, &
                       ns,fi,fs,ts, errcode, errtag, ext,format_str, &
                       case_file,geo_file, ipart, spart,spart_fs, buffer, node_hex8, gnum_hex8, &
                       gnum_quad4,node_quad4)
call sync_process


! Actually apply non-dimensionalisation
call apply_nondimensionalisation

! compute this after nondimensionlization
call compute_max_elementsize()



if(myrank==0)then
  write(logunit,'(a)')'--------------------------------------------'
  flush(logunit)
endif

! solver type
if(solver_type.eq.smart_solver)then
  if(ismpi)then
    solver_type=petsc_solver
  else
    solver_type=builtin_solver
  endif
endif
if(ismpi)then
  if(solver_type.eq.builtin_solver)then
    if(myrank==0)write(logunit,'(a)')'solver type: builtin parallel solver'
  elseif(solver_type.eq.petsc_solver)then
    if(myrank==0)write(logunit,'(a)')'solver type: PETSc parallel solver'
  else
    write(errtag,'(a,i4)')'ERROR: invalid solver type:',solver_type
    call control_error(errcode,errtag,logunit,myrank)
  endif
else
  if(solver_type.eq.builtin_solver)then
    write(logunit,'(a)')'solver type: builtin serial solver'
  elseif(solver_type.eq.petsc_solver)then
    write(errtag,'(a,i4)')'ERROR: PETSc solver must be run with MPI:',solver_type
    call control_error(errcode,errtag,logunit,myrank)
  else
    write(errtag,'(a,i4)')'ERROR: invalid solver type:',solver_type
    call control_error(errcode,errtag,logunit,myrank)
  endif
endif


write(logunit,*)'CALLING SPECFEM3D'
! call main routines
call specfem3d()
!if(nexcav==0)then
!  ! slope stability
!  !call semslope3d()
!else
!  ! excavation
!  !call semexcav3d()
!endif
!-----------------------------------


! Clean up and deallocate 
if(ISDISP_DOF)then
  deallocate(edofu)
endif
if(ISPOT_DOF)then
  deallocate(edofphi)
endif
! clean up                                                                       
call cleanup_model(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)

call cleanup_gll1d()

call cleanup_hexface(errcode,errtag)                                             
call control_error(errcode,errtag,stdout,myrank)

call cleanup_integration(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)

call cleanup_integration2d(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)

call cleanup_free_surface()

! compute elapsed time
call cpu_time(cpu_tend)
telap=cpu_tend-cpu_tstart
max_telap=maxscal(telap)
mean_telap=sumscal(telap)/real(nproc,kreal)

if(myrank==0)then
  write(format_str,*)ceiling(log10(real(max_telap)+1.))+5 ! 1 . and 4 decimals
  format_str='(3(f'//trim(adjustl(format_str))//'.4,1X))'
  write(logunit,'(a)')'ELAPSED TIME, MAX ELAPSED TIME, MEAN ELAPSED TIME'
  write(logunit,fmt=format_str)telap,max_telap,mean_telap
  write(logunit,'(a)')'--------------------------------------------'
  flush(logunit)
  close(logunit)
endif
!-----------------------------------

if(myrank==0)then
  inquire(stdout,opened=isopen)
  if(isopen)close(stdout)
endif

errcode=0

call sync_process
call close_process()
contains
!-------------------------------------------------------------------------------







! Routine must be called before converting hex8 to spectral elements
subroutine compute_max_elementsize
use math_library,only:distance
implicit none
integer :: i_elmt,mdomain,num(8)
real(kind=kreal),dimension(ndim) :: x1,x2,x3,x4,x5,x6,x7,x8
real(kind=kreal) :: d1,d2,d3,d4
real(kind=kreal) :: maxsize,maxdiag
! find largest size of the element
maxsize=zero
do i_elmt=1,nelmt
  ! discard transition and infinite elements
  mdomain=mat_domain(mat_id(i_elmt))
  if(mdomain==ELASTIC_TRINFDOMAIN .or.      &
     mdomain==ELASTIC_INFDOMAIN   .or.      &
     mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
     mdomain==VISCOELASTIC_INFDOMAIN)cycle

  num=g_num(hex8_gnode,i_elmt)
  x1=g_coord(:,num(1))
  x2=g_coord(:,num(2))
  x3=g_coord(:,num(3))
  x4=g_coord(:,num(4))
  x5=g_coord(:,num(5))
  x6=g_coord(:,num(6))
  x7=g_coord(:,num(7))
  x8=g_coord(:,num(8))
  d1=distance(x1,x7,3)
  d2=distance(x2,x8,3)
  d3=distance(x3,x5,3)
  d4=distance(x4,x6,3)
  maxdiag=max(d1,d2,d3,d4)
  if(maxdiag.gt.maxsize)maxsize=maxdiag
enddo
maxsize_elmt=maxscal(maxsize)
sqmaxsize_elmt=maxsize_elmt*maxsize_elmt
if(myrank==0)then
  write(logunit,'(a,g0.6)')'maximum element size across the diagonal: ', &
  maxsize_elmt*DIM_L
  flush(logunit)
endif
end subroutine compute_max_elementsize
!===============================================================================

end program specfemx
!===============================================================================
