! This is a main program SPECFEMX
! REVISION:
!  HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
!  WE Jun 8 2022 - cleaned up driver code into series of functions
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
use sea_level
use ice
#if (USE_MPI)
use mpi_library
use math_library_mpi
use output_to_user
#else
use serial_library
use math_library_serial
#endif
use nondimensionalisation
use user_input
use save_mesh
use postprocess,only:write_scalar_to_file_freesurf
use write_ensight
use print_w_MPI
use cleanup
use det_solver


!******************************************************************************
!******************************************************************************

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

!******************************************************************************
!******************************************************************************



! Start up MPI 
call start_process()

! Read cmd line/process input file etc ... 
call process_user_input(cmd, tdate, ttime, tzone, ios, path, &
                        ext, format_str, errcode, errtag, &
                        ismesh_only,arg1,arg2,inp_fname,prog, &
                        cpu_tstart)


ismesh_only = .false.





! Print info from read input for SL/Ice 
if(is_SL)then 
  call print_SL_read()
  call print_ice_read()
endif 



! Calculate model extents for individual processors/whole model
call calc_model_coord_extents(tot_nelmt,max_nelmt,min_nelmt, &
                              tot_nnode,max_nnode,min_nnode, &
                              absmaxx,absmaxy,absmaxz)


! Initialize model - allocates shearmod/bulkmod/massdensity arrays
call initialize_model(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)


! STILL NOT SURE WHAT THIS IS FOR! 
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


! Write original mesh to Ensight file 
call write_original_mesh(ipart, spart, npart, &
                         geo_file, infcase_file,& 
                         infgeo_file, &
                         add_tag, errcode, errtag, case_file,&
                         trinfcase_file,trinfgeo_file)


! If only saving mesh then quit program at this point. 
if(ismesh_only)then
  call close_process
endif

! ___________________________________________________________________
! store orginal connectivity which helps to identify ghost interfaces
allocate(g_num0(ngnode,nelmt))
g_num0=g_num

! precompute gll 1D
call precompute_gll1d()


! create spectral elements
call create_spec_elem(tot_nelmt,max_nelmt,min_nelmt, &
                      tot_nnode,max_nnode,min_nnode, &
                      errcode,errtag)




! number of elemental nodes (nodes per element = ngllx*nglly*ngllz
nenode=ngll 

! Reclassify (in)finite elements if there is an infinite BC
if(infbc)then
  ! free memory for modified arrays
  deallocate(g_num_finite,g_num_trinfinite,g_num_infinite, &
             node_finite,node_trinfinite,node_infinite)
  call classify_finite_infinite_elements(1)
endif

! initialize and set element DOFs
call initialize_dof()
call set_element_dof_uphi()


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
! Information is later used to determine the elevation 
! of the source point and to plot the free surface files.
call prepare_free_surface(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)
call sync_process()
! Calculate all of the nodes 
allnodesfs = sumscal(nnode_fs)
if(myrank.eq.0)then 
  write(*,*)'Number of FS nodes: ', allnodesfs
endif 

! Prepare SEA LEVEL IF NECESSARY
! Must be called after initialisation of Free Surface
! Reason for doing this separately is to keep all of the SL vars at the end
! by themselves because its only for surface elements. 
if(ISSL_DOF)then 
  call sea_level_dof()
endif 


case_file=trim(out_path)//trim(file_head)//trim(ptail)//'.case'
if(nexcav==0)then
  geo_file=trim(file_head)//trim(ptail)//'.geo'
else
  isgeo_change=.true.
  geo_file=trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.geo'
endif


! Add 1 time step to plot elastic and plastic results together
if(isplastic.and.nstep.le.1)then
  ns=ns+1
  dstep=one
endif
add_tag=''
call write_ensight_casefile_long(case_file,geo_file,add_tag,isgeo_change, &
ts,ns,fs,fi,twidth,errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)



! Save the mesh to the EnSight files. 
call save_mesh_ensight(infcase_file,infgeo_file,trinfcase_file, &
                       trinfgeo_file,isgeo_change,add_tag, twidth, &
                       fscase_file,fsgeo_file, fspcase_file,fspgeo_file, &
                       ns,fi,fs,ts, errcode, errtag, format_str, &
                       case_file,geo_file, ipart, spart,spart_fs, buffer, &
                       node_hex8, gnum_hex8, gnum_quad4,node_quad4)
call sync_process()


! Actually apply non-dimensionalisation
call apply_nondimensionalisation()

! compute (nondimensionalised) max element size. 
call compute_max_elementsize()

! Work out which solver to use
call determine_solver(errcode, errtag)

 
! Now, call main routine...
if(myrank.eq.0)then 
  write(logunit,*) 
  write(logunit,*) '**************** Starting SPECFEM3D  *****************'
  write(logunit,*) 
  write(logunit,*) 
endif 

call specfem3d()


! Print confirmation of completion/cpu time etc and cleanup
call print_completion_details(cpu_tstart,cpu_tend,telap,&
max_telap,mean_telap, format_str)
call run_cleanup_specfemx(errtag, errcode)
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
  write(logunit,'(a,g0.6)')'* Maximum element size across the diagonal: ', &
  maxsize_elmt*DIM_L
  write(logunit,*)
  flush(logunit)
endif
end subroutine compute_max_elementsize
!===============================================================================

end program specfemx
!===============================================================================
