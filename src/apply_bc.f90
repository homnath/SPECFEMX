! This routine consists of the Dirichlet boundary condition routines.
module bc
implicit none
character(len=250),private :: myfname=' => apply_bc.f90'
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

! This subroutine applies displacement boundary conditions and determines
! global degrees of freedom. the modification on the RHS vector due to the
! prescribed displacement field is done during stiffnesss computation.
! REVISION
!   HNG, Jul 12,2011; ; HNG, Apr 09,2010; WE 2022 Jun 18
subroutine apply_bc(errcode,errtag)
use global
use math_constants, only:zero
use element,only:hexface
use free_surface,only:nnode_fs,gnode_fs
use nondimensionpar,only:NONDIM_L
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: bctype,i,imat,ios,j,k
integer :: nelpart,i_elpart,i_node
integer :: ielmt,iface,idir
integer :: mdomain
real(kind=kreal) :: val
real(kind=kreal) :: dcoord(NDIM)
real(kind=kreal) :: diffx(NDIM),sqfix_radius
integer :: nfault,nfnode
integer,allocatable :: ifnode(:)
character(len=250) :: fname
character(len=150) :: data_path
logical :: nozero
logical,allocatable :: isnode(:),isnode0(:)

character(len=80) :: char80
integer :: dumi
! ufs must be real because it is read from the Ensight file
real,allocatable :: ufs(:,:)


errtag="ERROR: unknown!"
errcode=-1
! Set data path
if(ismpi.and.nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

nfnode=ngllx*nglly
allocate(ifnode(nfnode))

infinite_iface=.false.

if(ISDISP_DOF.and.isubc)then
  fname=trim(data_path)//trim(uxfile)//trim(ptail_inp)
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif

  allocate(isnode0(nnode),isnode(nnode))
  isnode0=.false.
  nfault=0
  bcux: do
    read(11,*,iostat=ios)bctype,val
    if (ios/=0)exit bcux
    if(val/=zero)nozero=.true.
    if(bctype==0)then ! point
      write(errtag,*)'ERROR: nodal displacement BC not implemented!'
      return
    elseif(bctype==1)then ! edge
      write(errtag,*)'ERROR: edge displacement BC not implemented!'
      return
    elseif(bctype==2)then ! face

      read(11,*)nelpart

      do i_elpart=1,nelpart
        read(11,*)ielmt,iface ! This will read a line and proceed to next line
        gdof(1,g_num(hexface(iface)%node,ielmt))=0

        bcnodalv(1,g_num(hexface(iface)%node,ielmt))=val
      enddo
    elseif(bctype==21)then ! face fault
      read(11,*)nelpart
      isnode=.false.
      nfault=nfault+1
      do i_elpart=1,nelpart
        read(11,*)ielmt,iface ! This will read a line and proceed to next line
        ! undo unsplit node
        ifnode=g_num(hexface(iface)%node,ielmt)
        do i_node=1,nfnode
          if(nfault.eq.2.and.isnode0(ifnode(i_node)))then
            ! undo for unslplit fault nodes
            gdof(1,ifnode(i_node))=1
            bcnodalv(1,ifnode(i_node))=zero
          elseif(.not.isnode(ifnode(i_node)))then
            gdof(1,ifnode(i_node))=1
            isnode(ifnode(i_node))=.true.
            gdof(1,ifnode(i_node))=0
            bcnodalv(1,ifnode(i_node))=val
          endif
        enddo
      enddo
      ! reset
      isnode0=isnode
      print*,'total fault points:',nfault,count(isnode),val
    else
      write(errtag,*)'ERROR: undefined displacement BC type ux!',bctype
      return
    endif
  enddo bcux
  close(11)

  fname=trim(data_path)//trim(uyfile)//trim(ptail_inp)
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif

  isnode0=.false.
  nfault=0
  bcuy: do
    read(11,*,iostat=ios)bctype,val
    if (ios/=0)exit bcuy
    if(val/=zero)nozero=.true.
    if(bctype==0)then ! point
      write(errtag,*)'ERROR: nodal displacement BC not implemented!'
      return
    elseif(bctype==1)then ! edge
      write(errtag,*)'ERROR: edge displacement BC not implemented!'
      return
    elseif(bctype==2)then ! face
      read(11,*)nelpart
      do i_elpart=1,nelpart
        ! This will read a line and proceed to next line
        read(11,*)ielmt,iface
        gdof(2,g_num(hexface(iface)%node,ielmt))=0
        bcnodalv(2,g_num(hexface(iface)%node,ielmt))=val
      enddo
    elseif(bctype==21)then ! face fault
      read(11,*)nelpart
      isnode=.false.
      nfault=nfault+1
      do i_elpart=1,nelpart
        ! This will read a line and proceed to next line
        read(11,*)ielmt,iface
        ! undo unsplit node
        ifnode=g_num(hexface(iface)%node,ielmt)
        do i_node=1,nfnode
          if(nfault.eq.2.and.isnode0(ifnode(i_node)))then
            ! undo for unsplit fault nodes
            gdof(2,ifnode(i_node))=1
            bcnodalv(2,ifnode(i_node))=zero
          elseif(.not.isnode(ifnode(i_node)))then
            gdof(2,ifnode(i_node))=1
            isnode(ifnode(i_node))=.true.
            gdof(2,ifnode(i_node))=0
            bcnodalv(2,ifnode(i_node))=val
          endif
        enddo
      enddo
      ! reset
      isnode0=isnode
      print*,'total fault points:',nfault,count(isnode),val
    else
      write(errtag,*)'ERROR: undefined displacement BC type uy!',bctype
      return
    endif
  enddo bcuy
  close(11)

  fname=trim(data_path)//trim(uzfile)//trim(ptail_inp)
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif
  isnode0=.false.
  nfault=0
  bcuz: do
    read(11,*,iostat=ios)bctype,val
    if (ios/=0)exit bcuz
    if(val/=zero)nozero=.true.
    if(bctype==0)then ! point
      write(errtag,*)'ERROR: nodal displacement BC not implemented!'
      return
    elseif(bctype==1)then ! edge
      write(errtag,*)'ERROR: edge displacement BC not implemented!'
      return
    elseif(bctype==2)then ! face
      read(11,*)nelpart
      do i_elpart=1,nelpart
        ! This will read a line and proceed to next line
        read(11,*)ielmt,iface
        gdof(3,g_num(hexface(iface)%node,ielmt))=0
        bcnodalv(3,g_num(hexface(iface)%node,ielmt))=val
      enddo
    elseif(bctype==21)then ! face fault
      read(11,*)nelpart
      isnode=.false.
      nfault=nfault+1
      do i_elpart=1,nelpart
        ! This will read a line and proceed to next line
        read(11,*)ielmt,iface
        ! undo unsplit node
        ifnode=g_num(hexface(iface)%node,ielmt)
        do i_node=1,nfnode
          if(nfault.eq.2.and.isnode0(ifnode(i_node)))then
            ! undo for unslplit fault nodes
            gdof(3,ifnode(i_node))=1
            bcnodalv(3,ifnode(i_node))=zero
          elseif(.not.isnode(ifnode(i_node)))then
            gdof(3,ifnode(i_node))=1
            isnode(ifnode(i_node))=.true.
            gdof(3,ifnode(i_node))=0
            bcnodalv(3,ifnode(i_node))=val
          endif
        enddo
      enddo
      ! reset
      isnode0=isnode
      print*,'total fault points:',nfault,count(isnode),val
    else
      write(errtag,*)'ERROR: undefined displacement BC type uz!',bctype
      return
    endif
  enddo bcuz
  deallocate(isnode,isnode0)
  close(11)
endif ! if(ISDISP_DOF)

! Fix center of the Earth.
if(fix_center)then
! We assume the center is at (0,0,0)
  if(fix_radius.gt.ZERO)then
    sqfix_radius = (NONDIM_L*fix_radius)*(NONDIM_L*fix_radius)
    do i_node=1,nnode
      ! compute distance
      diffx=g_coord(:,i_node)-center_coord
      if(dot_product(diffx,diffx).le.sqfix_radius)then
        gdof(:,i_node)=0
        bcnodalv(1,i_node)=ZERO
        write(logunit,*)'Center is fixed!'
      endif
    enddo
  else
    do i_node=1,nnode
      if(all(g_coord(:,i_node).eq.center_coord))then
        gdof(:,i_node)=0
        bcnodalv(1,i_node)=ZERO
        write(logunit,*)'Center is fixed!'
      endif
    enddo
  endif
endif

! Surface displacement defined on the surface SEM points 
if(ISDISP_DOF.and.isfsubc.and.nnode_fs>0)then
  fname=trim(ufspath)//trim(ufsfile)//trim(ptail_inp)//'.dis'
  open(unit=11,file=trim(fname),status='old',action='read',access='stream', &
  form='unformatted',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif
  allocate(ufs(nnode_fs,nndofu))
  read(11)char80
  read(11)char80
  read(11)dumi
  read(11)char80
  read(11)ufs
  close(11)
  gdof(idofu,gnode_fs)=0
  bcnodalv(idofu,gnode_fs)=NONDIM_L*transpose(ufs)
  deallocate(ufs)
endif


! Infinite boundary conditions
if(infbc)then
  fname=trim(data_path)//trim(inffile)//trim(ptail_inp)
  !print*,trim(fname)
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif
  bcinf: do
    read(11,*,iostat=ios)bctype,val
    if (ios/=0)exit
    if(val/=zero)nozero=.true.
    if(bctype==200)then ! infinite BC
      ! This will read a line and proceed to next line
      read(11,*)nelpart
      do i_elpart=1,nelpart
        read(11,*)ielmt,iface,idir
        ifnode=g_num(hexface(iface)%node,ielmt)
        mdomain=mat_domain(mat_id(ielmt))
        imat=mat_id(ielmt)
        if(ISDISP_DOF)then
        ! Although, transition/empty elements have no displacement DOFs, we have 
        ! to impose displacement BC to make the gdof IDs consistent across the 
        ! interfaces, in particular the interfaces between transition/empty
        ! elements and the solid elements. Imposing ZERO DOFs is equivalent to
        ! ignoring those DOFs anyway.

          !if(.not.isempty(imat))then
            gdof(idofu,ifnode)=0
            bcnodalv(idofu,ifnode)=val
          !endif
        endif
        if(ISPOT_DOF)then
          if( (idir.eq.1.and.isdxINF) .or. &
              (idir.eq.2.and.isdyINF) .or. &
              (idir.eq.3.and.isdzINF) )then
            ! set nonzero traction values (NOT YET SUPPORTED)
            ! do nothing for zero traction
          else
            ! Set Dirichlet BC on the infinite face
            !gdof(idofphi,g_num(hexface(iface)%node,ielmt))=0
            !bcnodalv(idofphi,g_num(hexface(iface)%node,ielmt))=val
            gdof(idofphi,ifnode)=0
            bcnodalv(idofphi,ifnode)=val
          endif
        endif
        infinite_iface(iface,ielmt)=.true.
        infinite_face_idir(iface,ielmt)=idir
      enddo
    else
      write(errtag,*)'ERROR: undefined displacement BC type for INF!',bctype
      return
    endif
  enddo bcinf
  close(11)
endif

if(myrank.eq.0)then
  write(logunit,*)' ✓  Applied boundary conditions'
  write(logunit,*)
endif
errcode=0

return
end subroutine apply_bc
!===========================================================================

!WE - applies non-zero boundary conditions originally in specfem3d.f90
subroutine apply_nonzero_bc()
  use global
  use local
  use set_precision
  use math_constants
  implicit none 

  ! Local 
  integer :: i_elmt, i,  ielmt, iedof, j_node, i_dof, j_dof, idof


  ! Modify RHS vector for prescribed displacements
  ! i.e. if boundary dispalcements are not equal to zero
  ! WARNING: need to check for nedofu
  ! MAKES ubcload = ubcload - K * U 
  ! i.e. before we had KU = F where F = the ubcload bit. 
  do i_elmt=1,nelmt
    ielmt=i_elmt ! all elements
    num=g_num(:,ielmt)
    egdof=gdof_elmt(:,ielmt)
   
    kmat=storekmat(:,:,ielmt)
    iedof=0
    do j_node=1,nenode
      do i_dof=1,nndof !nndofu
        iedof=iedof+1
        if(bcnodalv(i_dof,num(j_node))/=ZERO)then
          ubcload(egdof)=ubcload(egdof)-kmat(:,iedof)*bcnodalv(i_dof,num(j_node))
        endif
      enddo
    enddo
  enddo ! i_elmt


  ! Copies the non-zero values of bcnodalv to the nodalu vector
  ! For displacement only. 
  if(ISDISP_DOF)then
    do i_dof=1,nndofu !WE e.g. 3 for a 3D disp vector?
      idof=idofu(i_dof)
      do j_dof=1,nnode
        if(bcnodalv(idof,j_dof)/=ZERO)then 
            nodalu(i_dof,j_dof)=bcnodalv(idof,j_dof)
        endif 
      enddo
    enddo
  endif


  ! set BC nodal potential to nodalphi array
  ! DOESNT DO ANYTHING TO THE A 'LOAD' vector (e.g. ubcload above) - user cant set phi BC? 
  if(ISPOT_DOF)then
    do i_dof=1,nndofphi
      idof=idofphi(i_dof)
      do j_dof=1,nnode
        if(bcnodalv(idof,j_dof)/=ZERO)nodalphi(j_dof)=bcnodalv(idof,j_dof)
      enddo
    enddo
  endif

end subroutine apply_nonzero_bc
!===========================================================================

end module bc
!===========================================================================
