! DESCRIPTION
!  This module conains the routines to set model properties.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!  HNG, Feb 19, 2016; HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
!  WE Jun 8 2022 - added calc_model_coord_extents
!TODO
!  -
module model
use set_precision
implicit none
! private
character(len=250),private :: myfname=' => model.f90'
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

subroutine initialize_model(errcode,errtag)
use global,only:myrank,NDIM,ngll,nelmt,ISDISP_DOF,ISPOT_DOF, &
                POT_TYPE,PGRAVITY,PMAGNETIC,PELECTRIC,PCHARGE, &
                isbulkmod,isshearmod,ismassdens, &
                ismagnetization,iselectric,ischarge, &
                bulkmod_elmt,shearmod_elmt,massdens_elmt,magnetization_elmt,&
                econductivity1_elmt,econductivity2_elmt,econductivity3_elmt, &
                econalpha_elmt,econbeta_elmt,econgamma_elmt, &
                econductivity_aniso, &
                charge_density_elmt,logunit
use math_constants,only:ZERO
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

errtag="ERROR: unknown!"
errcode=-1


isbulkmod=.false.
isshearmod=.false.
ismassdens=.false.
ismagnetization=.false.
iselectric=.false.
ischarge=.false.
if(ISDISP_DOF)then
  isbulkmod=.true.
  isshearmod=.true.
  allocate(bulkmod_elmt(ngll,nelmt),shearmod_elmt(ngll,nelmt))
  bulkmod_elmt=ZERO
  shearmod_elmt=ZERO
endif

! Mass density
if(ISDISP_DOF .or. (ISPOT_DOF.and.POT_TYPE==PGRAVITY))then
  ismassdens=.true.
  allocate(massdens_elmt(ngll,nelmt))
  massdens_elmt=ZERO
endif

! Magnetization
if(ISPOT_DOF.and.POT_TYPE==PMAGNETIC)then
  ismagnetization=.true.
  allocate(magnetization_elmt(NDIM,ngll,nelmt))
  magnetization_elmt=ZERO
endif
! Electrical conductivity
if(ISPOT_DOF.and.POT_TYPE==PELECTRIC)then
  iselectric=.true.
  allocate(econductivity1_elmt(ngll,nelmt))
  econductivity1_elmt=ZERO
  if(econductivity_aniso)then
    ! anisotropic
    allocate(econductivity2_elmt(ngll,nelmt),econductivity3_elmt(ngll,nelmt))
    econductivity2_elmt=ZERO
    econductivity3_elmt=ZERO
    allocate(econalpha_elmt(ngll,nelmt),econbeta_elmt(ngll,nelmt), &
    econgamma_elmt(ngll,nelmt))
    econalpha_elmt=ZERO
    econbeta_elmt=ZERO
    econgamma_elmt=ZERO
  endif
endif
! Charge density
if(ISPOT_DOF.and.POT_TYPE==PCHARGE)then
  ischarge=.true.
  allocate(charge_density_elmt(ngll,nelmt))
  charge_density_elmt=ZERO
endif
errcode=0

if(myrank==0)then
  write(logunit,*)'✓ Completed model initialisation'
  write(logunit,*) 
endif 
end subroutine initialize_model
!===============================================================================

subroutine cleanup_model(errcode,errtag)
use global,only:NDIM,ngll,nelmt,ISDISP_DOF,ISPOT_DOF, &
                POT_TYPE,PGRAVITY,PMAGNETIC, &
                isbulkmod,isshearmod,ismassdens, &
                ismagnetization,iselectric,ischarge, &
                bulkmod_elmt,shearmod_elmt,massdens_elmt,magnetization_elmt, &
                econductivity1_elmt,econductivity2_elmt,econductivity3_elmt, &
                econalpha_elmt,econbeta_elmt,econgamma_elmt, &
                econductivity_aniso, &
                charge_density_elmt
use math_constants,only:ZERO

implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

errtag="ERROR: unknown!"
errcode=-1

if(isbulkmod)deallocate(bulkmod_elmt)
if(isshearmod)deallocate(shearmod_elmt)
if(ismassdens)deallocate(massdens_elmt)
if(ismagnetization)deallocate(magnetization_elmt)
if(iselectric)then
  deallocate(econductivity1_elmt)
  if(econductivity_aniso)then
    deallocate(econductivity2_elmt,econductivity3_elmt)
    deallocate(econalpha_elmt,econbeta_elmt,econgamma_elmt)
  endif
endif
if(ischarge)deallocate(charge_density_elmt)
errcode=0
end subroutine cleanup_model
!===============================================================================

subroutine write_model_cell(errcode,errtag)
use math_constants,only:INFTOL,FOUR_THIRD,HALF,ONE,TWO,ZERO
use global
use write_ensight,only:ensight_hex8, &
                write_ensight_perelementSCALAS, &
                write_ensight_perelementSCALAS_part1, &
                write_ensight_perelementVECAS, &
                write_ensight_perelementVECAS_part1
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt
integer :: iblk,imatmag

integer :: npart
integer,allocatable :: ipart(:)
character(len=80),allocatable :: spart(:) ! this must be 80 characters long
character(len=250) :: out_fname

integer :: i
real(kind=kreal),allocatable :: rho_elmt(:),M_elmt(:,:),Mmag_elmt(:), &
econd_elmt(:),erho_elmt(:)

errtag="ERROR: unknown!"
errcode=-1

! Mass density
if(ismassdens)then
  allocate(rho_elmt(nelmt))
  rho_elmt=-9999_kreal

  do i_elmt=1,nelmt
    iblk=mat_id(i_elmt)
    rho_elmt(i_elmt)=rho_blk(iblk)
  enddo
endif

! Magnetization
if(ismagnetization)then
  !allocate(M_elmt(NDIM,nelmt))
  !M_elmt=ZERO
  allocate(Mmag_elmt(nelmt))
  Mmag_elmt=ZERO

  do i_elmt=1,nelmt
    iblk=mat_id(i_elmt)
    if(ismagnet_blk(iblk))then
      imatmag=imat_to_imatmag(iblk)
      Mmag_elmt(i_elmt)=Mmag_blk(imatmag)
    endif
  enddo
endif

! Electrical Conductivity
if(iselectric)then
  allocate(econd_elmt(nelmt))
  econd_elmt=ZERO

  do i_elmt=1,nelmt
    iblk=mat_id(i_elmt)
    if(iselectric_blk(iblk))then
      econd_elmt(i_elmt)=econductivity1_blk(iblk)
    endif
  enddo
endif

! Charge Density
if(ischarge)then
  allocate(erho_elmt(nelmt))
  erho_elmt=ZERO

  do i_elmt=1,nelmt
    iblk=mat_id(i_elmt)
    if(ischarge_blk(iblk))then
      erho_elmt(i_elmt)=charge_density_blk(iblk)
    endif
  enddo
endif
! save model
if(infbc)then
  npart=2
  allocate(ipart(npart),spart(npart))
  ipart=(/ (i,i=1,npart) /)
  spart(1)='finite region'
  spart(2)='infinite region'
else
  npart=1
  allocate(ipart(npart),spart(npart))
  ipart=(/ (i,i=1,npart) /)
  spart(1)='finite region'
endif

! Mass density
if(ismassdens)then
  if(infbc)then
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.rho'
    call write_ensight_perelementSCALAS_part1(out_fname,ensight_hex8,ipart,spart,1, &
    nelmt_finite,elmt_finite,nelmt,real(rho_elmt))
  else
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.rho'
    call write_ensight_perelementSCALAS(out_fname,ensight_hex8,ipart,spart, &
    nelmt,real(rho_elmt))
  endif
  deallocate(rho_elmt)
endif
! Magnetization
if(ismagnetization)then
  if(infbc)then
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.mag'
    call write_ensight_perelementSCALAS_part1(out_fname,ensight_hex8,ipart,spart,1, &
    nelmt_finite,elmt_finite,nelmt,real(Mmag_elmt))
  else
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.mag'
    call write_ensight_perelementSCALAS(out_fname,ensight_hex8,ipart,spart, &
    nelmt,real(Mmag_elmt))
  endif
  deallocate(Mmag_elmt)
  !if(infbc)then
  !  write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.mag'
  !  call write_ensight_perelementVECAS_part1(out_fname,ensight_hex8,ipart,spart,1, &
  !  nelmt_finite,elmt_finite,NDIM,nelmt,real(M_elmt))
  !else
  !  write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.mag'
  !  call write_ensight_perelementVECAS(out_fname,ensight_hex8,ipart,spart, &
  !  NDIM,nelmt,real(M_elmt))
  !endif
  !deallocate(M_elmt)
endif
! Electrical Conductivity
if(iselectric)then
  if(infbc)then
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.econ'
    call write_ensight_perelementSCALAS_part1(out_fname,ensight_hex8,ipart,spart,1, &
    nelmt_finite,elmt_finite,nelmt,real(econd_elmt))
  else
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.econ'
    call write_ensight_perelementSCALAS(out_fname,ensight_hex8,ipart,spart, &
    nelmt,real(econd_elmt))
  endif
  deallocate(econd_elmt)
endif
! Charge Density
if(ischarge)then
  if(infbc)then
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.erho'
    call write_ensight_perelementSCALAS_part1(out_fname,ensight_hex8,ipart,spart,1, &
    nelmt_finite,elmt_finite,nelmt,real(erho_elmt))
  else
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.erho'
    call write_ensight_perelementSCALAS(out_fname,ensight_hex8,ipart,spart, &
    nelmt,real(erho_elmt))
  endif
  deallocate(econd_elmt)
endif
errcode=0
end subroutine write_model_cell
!===============================================================================

subroutine set_model_properties(errcode,errtag)
use math_constants,only:INFTOL,FOUR_THIRD,HALF,ONE,TWO,ZERO
use global
use shape_library,only:shape_function_hex8p
use write_ensight,only:write_ensight_pernodeSCALAS,write_ensight_pernodeSCALAS_part1
#if (USE_MPI)
use ghost_library_mpi
use output_to_user
#else
use serial_library
#endif
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i,i_blk,i_dim,i_elmt,i_gll, ios
integer :: iblk,ielmt
integer :: imat,imatmag
integer,allocatable :: ielmts(:),block_nelmt(:)

integer :: npart
integer,allocatable :: ipart(:)
character(len=80),allocatable :: spart(:) ! this must be 80 characters long
character(len=250) :: out_fname

! custome model
real(kind=kreal) :: drho,drho03,depthp,fac
real(kind=kreal) :: zp 

integer,allocatable :: num(:),nvalency(:)
real(kind=kreal),allocatable :: bulkmod_node(:),shearmod_node(:),rho_node(:)

character(len=250) :: fname

type material_block
  integer,dimension(:),allocatable :: elmt
  ! whether the node on the interface is intact (.true.) or void (.false.)
end type material_block
type(material_block),dimension(:),allocatable :: block

errtag="ERROR: unknown!"
errcode=-1

allocate(num(ngll))
!! Elastic moduli
!if(ISDISP_DOF)then
!  allocate(bulkmod_elmt(ngll,nelmt),shearmod_elmt(ngll,nelmt))
!  bulkmod_elmt=ZERO
!  shearmod_elmt=ZERO
!endif
!! Mass density
!if(ISDISP_DOF .or. (ISPOT_DOF.and.POT_TYPE==PMAGNETIC))then
!  allocate(massdens_elmt(ngll,nelmt))
!  massdens_elmt=ZERO
!endif
!! Magnetization
!if(ISPOT_DOF.and.POT_TYPE==PMAGNETIC)then
!  allocate(magnetization_elmt(NDIM,ngll,nelmt))
!  magnetization_elmt=ZERO
!endif

! classify material blocks
! count material blocks


allocate(block_nelmt(nmatblk))

block_nelmt=0
do i_blk=1,nmatblk
  block_nelmt(i_blk)=count(mat_id==i_blk)
enddo

allocate(block(nmatblk))
do i_blk=1,nmatblk
    allocate(block(i_blk)%elmt(block_nelmt(i_blk)))
enddo

allocate(ielmts(nmatblk))

ielmts=0
do i_elmt=1,nelmt
  iblk=mat_id(i_elmt)
  ielmts(iblk)=ielmts(iblk)+1
  block(iblk)%elmt(ielmts(iblk))=i_elmt

enddo

! Convert block model to point model
matblock: do i_blk=1,nmatblk
  
! block model
  if(type_blk(i_blk).eq.0)then
    ! Shear/bulk modulus
    if(ISDISP_DOF)then
      bulkmod_elmt(:,block(i_blk)%elmt)=bulkmod_blk(i_blk)
      shearmod_elmt(:,block(i_blk)%elmt)=shearmod_blk(i_blk)
    endif

    ! Mass density if using gravity 
    if(ISDISP_DOF .or. (ISPOT_DOF.and.POT_TYPE==PGRAVITY))then
      massdens_elmt(:,block(i_blk)%elmt)=rho_blk(i_blk)

      if(trim(cmodel)=='chakravarthi')then
        ! make sure that the block for gradient anomaly must be assigen block ID 1 
        if(i_blk==1)then
          if(myrank==0)then
            write(logunit,'(a,g0.6,1x,g0.6,1x,g0.6)')'gradient density model: \rho0, \alpha, z0:',&
            cmodel_drho0,cmodel_alpha,cmodel_zref
            flush(logunit)
          endif
          ! assign gradient model
          drho03=cmodel_drho0*cmodel_drho0*cmodel_drho0
          do i_elmt=1,block_nelmt(i_blk)
            ielmt=block(i_blk)%elmt(i_elmt)
            num=g_num(:,ielmt)
            do i_gll=1,ngll
              zp=g_coord(3,num(i_gll))

              depthp=cmodel_zref-zp ! positive down
              fac=cmodel_drho0-cmodel_alpha*depthp
              if(fac.eq.ZERO)then
                write(errtag,'(a,g0.6,a)')'ERROR: gradient density model is &
                &singular at depth ',zp,' !'
                return
              endif
              drho=drho03/(fac*fac)

              massdens_elmt(i_gll,i_elmt)=drho
            enddo
          enddo

        endif ! i_blk
      endif !(trim(cmodel)=='chakravarthi')
    endif

    ! magnetization
    if(ISPOT_DOF.and.POT_TYPE==PMAGNETIC)then
      if(ismagnet_blk(i_blk))then
        imatmag=imat_to_imatmag(i_blk)
        do i_gll=1,NGLL
          do i_dim=1,NDIM
            magnetization_elmt(i_dim,i_gll,block(i_blk)%elmt)= &
              magnetization_blk(i_dim,imatmag)
          enddo ! i_dim
        enddo ! i_gll
      endif
    endif
    
    ! electrial conductivity
    if(ISPOT_DOF.and.POT_TYPE==PELECTRIC)then
      if(iselectric_blk(i_blk))then
        econductivity1_elmt(:,block(i_blk)%elmt) = econductivity1_blk(i_blk)
        if(econductivity_aniso)then
          econductivity2_elmt(:,block(i_blk)%elmt) = econductivity2_blk(i_blk)
          econductivity3_elmt(:,block(i_blk)%elmt) = econductivity3_blk(i_blk)
          econalpha_elmt(:,block(i_blk)%elmt) = econalpha_blk(i_blk)
          econbeta_elmt(:,block(i_blk)%elmt) = econbeta_blk(i_blk)
          econgamma_elmt(:,block(i_blk)%elmt) = econgamma_blk(i_blk)
        endif
      endif
    endif
    
    ! electrial charge density
    if(ISPOT_DOF.and.POT_TYPE==PCHARGE)then
      if(ischarge_blk(i_blk))then
        charge_density_elmt(:,block(i_blk)%elmt)=charge_density_blk(i_blk)
      endif
    endif

  ! tomographic structured grid model
  elseif(type_blk(i_blk).eq.-1)then
    call convert_tomo_to_point_model(i_blk, num , nvalency, ios, &
    bulkmod_node, shearmod_node, rho_node, errcode, errtag)
  else
    print*,'ERROR: unsupported type_blk:',type_blk(i_blk)
  endif

enddo matblock

! model_type=='gll'
if(trim(cmodel).eq.'gll')then
  ! open GLL model files and read
  fname=trim(out_path)//trim(file_head)//'_pdensity_gll_'//trim(ptail)
  open(unit=11,file=trim(fname),access='stream',form='unformatted', &
  status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif
  !read(11,*)bulkmod_elmt
  !read(11,*)shearmod_elmt
  if(.not.ISDISP_DOF.and.ISPOT_DOF)then
    read(11)massdens_elmt
  else
    write(*,*)'ERROR: GLL model is not yet unviersally implemented!'
    stop
  endif
  close(11)
endif

deallocate(block_nelmt,ielmts)
do i_blk=1,nmatblk
  deallocate(block(i_blk)%elmt)
enddo
deallocate(block)


! Write Savedata as a seperate function
! Save model only for finite region.
! NOTE:
!   Plotting model with node is less accurate than the cell. Cell representation
!   can preserve the model discontinuities if any.
if(savedata%model)then
  allocate(bulkmod_node(nnode),shearmod_node(nnode),rho_node(nnode))

  ! Compute nodal values.
  ! This plot is not exact since the values have to be interpolated on the 
  ! nodes for the exact plot we may have to define multiblk data.
  allocate(nvalency(nnode))
  nvalency=0
  bulkmod_node=ZERO
  shearmod_node=ZERO
  rho_node=ZERO
  do i_elmt=1,nelmt
    ! Skip transition and infinite elements
    imat=mat_id(i_elmt)
    if(mat_domain(imat).ge.ELASTIC_TRINFDOMAIN)cycle

    num=g_num(:,i_elmt)
    bulkmod_node(num) = bulkmod_node(num) + bulkmod_elmt(:,i_elmt)
    shearmod_node(num) = shearmod_node(num) + shearmod_elmt(:,i_elmt)
    rho_node(num) = rho_node(num) + massdens_elmt(:,i_elmt)
    nvalency(num)=nvalency(num) + 1
  enddo
  ! Assemble across the processors.
  call assemble_ghosts_nodal_iscalar(nvalency,nvalency)
  call assemble_ghosts_nodal_fscalar(bulkmod_node,bulkmod_node)
  call assemble_ghosts_nodal_fscalar(shearmod_node,shearmod_node)
  call assemble_ghosts_nodal_fscalar(rho_node,rho_node)

  bulkmod_node = bulkmod_node/nvalency
  shearmod_node = shearmod_node/nvalency
  rho_node = rho_node/nvalency
  !print*,myrank,minval(nvalency),maxval(massdens_elmt),maxval(rho_node)
  !print*,myrank,minval(nvalency),maxval(bulkmod_elmt),maxval(bulkmod_node)
  !print*,myrank,minval(nvalency),maxval(shearmod_elmt),maxval(shearmod_node)
  deallocate(nvalency)

  if(infbc)then
    npart=2
    allocate(ipart(npart),spart(npart))
    ipart=(/ (i,i=1,npart) /)
    spart(1)='finite region'
    spart(2)='infinite region'
  else
    npart=1
    allocate(ipart(npart),spart(npart))
    ipart=(/ (i,i=1,npart) /)
    spart(1)='finite region'
  endif
  if(infbc)then
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//trim(ptail)//'.kappa'
    call write_ensight_pernodeSCALAS_part1(out_fname,ipart,spart,1, &
    nnode_finite,node_finite,nnode,real(bulkmod_node))
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//trim(ptail)//'.mu'
    call write_ensight_pernodeSCALAS_part1(out_fname,ipart,spart,1, &
    nnode_finite,node_finite,nnode,real(shearmod_node))

    write(out_fname,'(a)')trim(out_path)//trim(file_head)//trim(ptail)//'.rho'
    call write_ensight_pernodeSCALAS_part1(out_fname,ipart,spart,1, &
    nnode_finite,node_finite,nnode,real(rho_node))
  else
    write(out_fname,'(a)')trim(out_path)//trim(file_head)//trim(ptail)//'.kappa'
    call write_ensight_pernodeSCALAS(out_fname,ipart,spart, &
    nnode,real(bulkmod_node))

    write(out_fname,'(a)')trim(out_path)//trim(file_head)//trim(ptail)//'.mu'
    call write_ensight_pernodeSCALAS(out_fname,ipart,spart, &
    nnode,real(shearmod_node))

    write(out_fname,'(a)')trim(out_path)//trim(file_head)//trim(ptail)//'.rho'
    call write_ensight_pernodeSCALAS(out_fname,ipart,spart, &
    nnode,real(rho_node))
  endif
  deallocate(bulkmod_node,shearmod_node,rho_node)
endif


deallocate(num)
if(myrank.eq.0)write(logunit, *)'Completed set_model_properties...'

errcode=0
end subroutine set_model_properties
!===============================================================================

subroutine convert_tomo_to_point_model(i_blk, num , nvalency, ios, &
  bulkmod_node, shearmod_node, rho_node, errcode, errtag)
  ! WE created this subroutine that was originally part of set_model_properties
  ! May not be functional - may need some edits for input/output variables

  ! USES: 
  use global
  use math_constants
  use shape_library,only:shape_function_hex8p

  ! IO variables
  integer :: i_blk, ios
  integer,allocatable :: num(:),nvalency(:)
  real(kind=kreal), allocatable :: bulkmod_node(:), shearmod_node(:),rho_node(:)
  integer :: errcode
  character(len=250) :: errtag

  ! Local variables
  integer :: i_elmt, i_grid, i_gll
  integer :: ic(8)
  integer :: ix1(3),ix2(3)
  integer :: grid_l1,grid_l2,grid_m1,grid_m2,grid_n1,grid_n2
  integer :: grid_n,grid_nx,grid_ny,grid_nz,grid_nxy
  real(kind=kreal) :: grid_x0(3),grid_x(3),grid_x1(3),grid_dx(3)
  real(kind=kreal) :: xp(3),midx(3)
  real(kind=kreal) :: grid_vpmin,grid_vsmin,grid_rhomin,grid_qpmin,grid_qsmin
  real(kind=kreal) :: grid_vpmax,grid_vsmax,grid_rhomax,grid_qpmax,grid_qsmax
  real(kind=kreal),allocatable :: grid_vp(:),grid_vs(:),grid_rho(:)
  real(kind=kreal) :: vp,vs,rho
  real(kind=kreal) :: shape_hex8(8)
 
  ! CODE: 
    ! Read file 
  write(logunit, *)'Converting tomo. model to pointwise...'

    open(unit=11,file=trim(inp_path)//trim(mfile_blk(i_blk)),status='old',     &
      action='read',iostat = ios)

      ! Err if cant read
      if( ios /= 0 ) then
        write(errtag,'(a)')'ERROR: file "'//trim(mfile_blk(i_blk))//'" cannot be opened!'
        return
      endif

      ! Read tomographic data
      read(11,*)grid_x0,grid_x1
      read(11,*)grid_dx
      read(11,*)grid_nx,grid_ny,grid_nz
      !if(myrank==0)print*,'grid_x0:',grid_x0
      !if(myrank==0)print*,'grid_x1:',grid_x1
      !if(myrank==0)print*,'grid_dx:',grid_dx
      !if(myrank==0)print*,'grid_nx:',grid_nx,grid_ny,grid_nz
      grid_l1=1
      grid_m1=1
      grid_n1=1

      grid_l2=grid_nx
      grid_m2=grid_ny
      grid_n2=grid_nz

      grid_n=grid_nx*grid_ny*grid_nz
      grid_nxy=grid_nx*grid_ny
      allocate(grid_vp(grid_n),grid_vs(grid_n),grid_rho(grid_n))
      grid_vp=-inftol
      grid_vs=-inftol
      grid_rho=-inftol

      grid_vpmax=-inftol
      grid_vpmin=-inftol
      grid_vsmax=-inftol
      grid_vsmin=-inftol
      grid_rhomax=-inftol
      grid_rhomin=-inftol
      grid_qpmax=-inftol
      grid_qpmin=-inftol
      grid_qsmax=-inftol
      grid_qsmin=-inftol

      read(11,*,iostat=ios)grid_vpmin,grid_vpmax,grid_vsmin,grid_vsmax,          &
      grid_rhomin,grid_rhomax,grid_qpmin,grid_qpmax,grid_qsmin,grid_qsmax
      do i_grid=1,grid_n
        read(11,*,iostat=ios)grid_x,grid_vp(i_grid),grid_vs(i_grid),grid_rho(i_grid)
      enddo
      close(11)

      ! check the properties read
      if(minval(grid_vp).lt.grid_vpmin .or. maxval(grid_vp).gt.grid_vpmax)then
        print*,'ERROR: read grid vp is beyond the given range!'
        stop
      endif
      if(minval(grid_vs).lt.grid_vsmin .or. maxval(grid_vs).gt.grid_vsmax)then
        print*,'ERROR: read grid vs is beyond the given range!'
        stop
      endif
      if(minval(grid_rho).lt.grid_rhomin .or. maxval(grid_rho).gt.grid_rhomax)then
        print*,'ERROR: read grid rho is beyond the given range!'
        stop
      endif

      ! interpolate the model
      do i_elmt=1,nelmt
        num=g_num(:,i_elmt)
        do i_gll=1,ngll
          xp=g_coord(:,num(i_gll))

          ix1=floor((xp-grid_x0)/grid_dx)+1

          if(ix1(1).le.grid_l1)ix1(1)=grid_l1
          if(ix1(2).le.grid_m1)ix1(2)=grid_m1
          if(ix1(3).le.grid_n1)ix1(3)=grid_n1

          if(ix1(1).ge.grid_l2)ix1(1)=grid_l2-1
          if(ix1(2).ge.grid_m2)ix1(2)=grid_m2-1
          if(ix1(3).ge.grid_n2)ix1(3)=grid_n2-1

          ix2=ix1+1

          ic(1)=(ix1(3)-1)*grid_nxy+(ix1(2)-1)*grid_nx+ix1(1)
          ic(2)=(ix1(3)-1)*grid_nxy+(ix1(2)-1)*grid_nx+ix2(1)
          ic(3)=(ix1(3)-1)*grid_nxy+(ix2(2)-1)*grid_nx+ix2(1)
          ic(4)=(ix1(3)-1)*grid_nxy+(ix2(2)-1)*grid_nx+ix1(1)
          ic(5)=(ix2(3)-1)*grid_nxy+(ix1(2)-1)*grid_nx+ix1(1)
          ic(6)=(ix2(3)-1)*grid_nxy+(ix1(2)-1)*grid_nx+ix2(1)
          ic(7)=(ix2(3)-1)*grid_nxy+(ix2(2)-1)*grid_nx+ix2(1)
          ic(8)=(ix2(3)-1)*grid_nxy+(ix2(2)-1)*grid_nx+ix1(1)
          if(maxval(ic).gt.grid_n .or. minval(ic).lt.1)then
            print*,'impossible!'
            print*,'grid nx:',grid_nx,grid_ny,grid_nz
            print*,'ix1:',ix1
            print*,'ix2:',ix2
            print*,'ic:',ic
            print*,'n',grid_n
            stop
          endif

          ! normalize point coordinates to natural coordinates
          midx=grid_x0+(ix1-1)*grid_dx+half*grid_dx

          ! shift origin to the cell center
          xp=xp-midx

          ! normalize to [-1, 1] range
          xp=TWO*xp/grid_dx
          where(xp.lt.-ONE)xp=-ONE
          where(xp.gt.ONE)xp=ONE
          ! compute shape function in natural coordinates
          call shape_function_hex8p(8,xp(1),xp(2),xp(3),shape_hex8)

          vp=dot_product(shape_hex8,grid_vp(ic))
          vs=dot_product(shape_hex8,grid_vs(ic))
          rho=dot_product(shape_hex8,grid_rho(ic))

          bulkmod_elmt(i_gll,i_elmt)=rho*(vp*vp-FOUR_THIRD*vs*vs)
          shearmod_elmt(i_gll,i_elmt)=rho*vs*vs
          massdens_elmt(i_gll,i_elmt)=rho
        enddo

      enddo
      deallocate(grid_vp,grid_vs,grid_rho)

end subroutine convert_tomo_to_point_model
!-------------------------------------------------------------------------------

subroutine calc_model_coord_extents(tot_nelmt,max_nelmt,min_nelmt, &
                                    tot_nnode,max_nnode,min_nnode, &
                                    absmaxx,absmaxy,absmaxz)

! USES
use global
use infinite_element
#if (USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif

implicit none
! IO variables
integer :: tot_nelmt,max_nelmt,min_nelmt,tot_nnode,max_nnode,min_nnode
real(kind=kreal) :: absmaxx,absmaxy,absmaxz

! Local variables
! none 

if(myrank.eq.0)then
  write(logunit,'(a)')'------------- Model Dimensions -------------'
endif

! Set coordinate extents of the finite model. This extent is later used in
! location routine.
if(infbc)then
  ! classify (in)finite elements
  call classify_finite_infinite_elements(0)
  ! Finite region
  ! Element and node count
  tot_nelmt=sumscal(nelmt_finite); tot_nnode=sumscal(nnode_finite)
  max_nelmt=maxscal(nelmt_finite); max_nnode=maxscal(nnode_finite)
  min_nelmt=minscal(nelmt_finite); min_nnode=minscal(nnode_finite)
  ! Coordinate extents of partitioned model
  pmodel_minx=minval(g_coord(1,node_finite))
  pmodel_maxx=maxval(g_coord(1,node_finite))
  pmodel_miny=minval(g_coord(2,node_finite))
  pmodel_maxy=maxval(g_coord(2,node_finite))
  pmodel_minz=minval(g_coord(3,node_finite))
  pmodel_maxz=maxval(g_coord(3,node_finite))
  ! Coordinate extents of finite model 
  model_minx=minscal(pmodel_minx)
  model_maxx=maxscal(pmodel_maxx)
  model_miny=minscal(pmodel_miny)
  model_maxy=maxscal(pmodel_maxy)
  model_minz=minscal(pmodel_minz)
  model_maxz=maxscal(pmodel_maxz)
  mincoord=min(model_minx,model_miny,model_minz)
  maxcoord=max(model_maxx,model_maxy,model_maxz)
  absmaxx=maxscal(maxval(abs(g_coord(1,:))))
  absmaxy=maxscal(maxval(abs(g_coord(2,:))))
  absmaxz=maxscal(maxval(abs(g_coord(3,:))))
  absmaxcoord=max(absmaxx,absmaxy,absmaxz)
  if(myrank==0)then
    write(logunit,'(a)')' * Finite region      : '
    write(logunit,'(a,i0)')'   --> total elements             : ',tot_nelmt
    write(logunit,'(a,i0)')'   --> min. elements per processor: ',min_nelmt
    write(logunit,'(a,i0)')'   --> max. elements per processor: ',max_nelmt
    write(logunit,*)
    write(logunit,'(a,i0)')'   --> total nodes                : ',tot_nnode
    write(logunit,'(a,i0)')'   --> min. nodes per processor   : ',min_nnode
    write(logunit,'(a,i0)')'   --> max. nodes per processor   : ',max_nnode
    write(logunit,*)
    write(logunit,'(a,g0.6,1x,g0.6,a)')'   --> Maximum x extent       : [',model_minx,model_maxx, ']'
    write(logunit,'(a,g0.6,1x,g0.6)')' y extent min max: ',model_miny,model_maxy
    write(logunit,'(a,g0.6,1x,g0.6)')' z extent min max: ',model_minz,model_maxz
    write(logunit,'(a,g0.6,1x,g0.6)')' min/max coord: ',mincoord,maxcoord
    write(logunit,'(a,g0.6)')' abs max coord: ',absmaxcoord
    flush(logunit)
  endif
else
  if(myrank.eq.0)then
    write(logunit,*)'No infinite bc (infbc = F)...'
  endif
  pmodel_minx=minval(g_coord(1,:))
  pmodel_maxx=maxval(g_coord(1,:))
  pmodel_miny=minval(g_coord(2,:))
  pmodel_maxy=maxval(g_coord(2,:))
  pmodel_minz=minval(g_coord(3,:))
  pmodel_maxz=maxval(g_coord(3,:))
   
endif

tot_nelmt=sumscal(nelmt); tot_nnode=sumscal(nnode)
max_nelmt=maxscal(nelmt); max_nnode=maxscal(nnode)
min_nelmt=minscal(nelmt); min_nnode=minscal(nnode)


! Coordinate extents of whole model
model_minx=minscal(minval(g_coord(1,:))); model_maxx=maxscal(maxval(g_coord(1,:)))
model_miny=minscal(minval(g_coord(2,:))); model_maxy=maxscal(maxval(g_coord(2,:)))
model_minz=minscal(minval(g_coord(3,:))); model_maxz=maxscal(maxval(g_coord(3,:)))
mincoord=min(model_minx,model_miny,model_minz)
maxcoord=max(model_maxx,model_maxy,model_maxz)
absmaxx=maxscal(maxval(abs(g_coord(1,:))))
absmaxy=maxscal(maxval(abs(g_coord(2,:))))
absmaxz=maxscal(maxval(abs(g_coord(3,:))))
absmaxcoord=max(absmaxx,absmaxy,absmaxz)

if(myrank==0)then
  write(logunit,'(a)')' * Whole model      : '
  write(logunit,'(a,i0)')'   --> total elements             : ',tot_nelmt
  write(logunit,'(a,i0)')'   --> min. elements per processor: ',min_nelmt
  write(logunit,'(a,i0)')'   --> max. elements per processor: ',max_nelmt
  write(logunit,*)
  write(logunit,'(a,i0)')'   --> total nodes                : ',tot_nnode
  write(logunit,'(a,i0)')'   --> min. nodes per processor   : ',min_nnode
  write(logunit,'(a,i0)')'   --> max. nodes per processor   : ',max_nnode
  write(logunit,*)       
  write(logunit,'(a,g0.6,2x,g0.6)')'   --> X extent                   : ',model_minx,model_maxx
  write(logunit,'(a,g0.6,2x,g0.6)')'   --> Y extent                   : ',model_miny,model_maxy
  write(logunit,'(a,g0.6,2x,g0.6)')'   --> Z extent                   : ',model_minz,model_maxz
  write(logunit,'(a,g0.6)')        '   --> Absolute maximum coord     : ',absmaxcoord
  write(logunit,*)
  flush(logunit)
endif

! Reassign pole coordinates if it is "center" of the model.
! NOTE: check if the center should be taken for the finite region only.
if(trim(pole0)=='center')then
  pole_coord0(1)=HALF*(model_minx+model_maxx)
  pole_coord0(2)=HALF*(model_miny+model_maxy)
  pole_coord0(3)=HALF*(model_minz+model_maxz)
endif

end subroutine calc_model_coord_extents 


end module model
!===============================================================================
