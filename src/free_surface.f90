module free_surface
use set_precision
implicit none
character(len=250),private :: myfname=' => free_surface.f90'
character(len=500),private :: errsrc

integer :: nelmt_fs,nnode_fs
integer,allocatable :: iface_fs(:), id_elem_fs(:), fs_elem_id(:) !stores faces, element IDs of fs elems
! We can store only gnum_fs and gnum4_fs can be later extracted from it.
integer,allocatable :: gnum4_fs(:,:)
integer,allocatable :: gnum_fs(:,:)
integer,allocatable :: gnode_fs(:)
! Connectivity renumberd to the free surface nodes.
integer,allocatable :: rgnum_fs(:,:)
contains
!-------------------------------------------------------------------------------

! This routine reads free surface information specified in the free surface file.
! NOTE: This routine must be called after the spectral elements are created.
! REVISION
!   HNG, Oct 08, 2018
subroutine prepare_free_surface(errcode,errtag)
use global,only:ismpi,nproc,nenode,maxngll2d,g_num, &
fsfile,inp_path,part_path,ptail_inp,savedata, nelmt
use element,only:hexface
use math_library,only:i_uniinv

implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: i,i_face,ios
integer :: ielmt,iface
integer :: num(nenode)

integer :: n1,n2
integer :: nsnode_all
integer,allocatable :: inode_order(:),nodelist(:)
logical,allocatable :: isnode(:)
character(len=80) :: fname
character(len=80) :: data_path



errtag=""
errcode=0
! Set data path
if(ismpi.and.nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

fname=trim(data_path)//trim(fsfile)//trim(ptail_inp)


open(unit=11,file=trim(fname),status='old',action='read',iostat=ios)
if (ios /= 0)then
  write(errtag,'(a)')'ERROR: input file "'//trim(fname)//'" cannot be opened!'
  errcode=-1
  return
endif

read(11,*,iostat=ios)nelmt_fs
if(ios/=0.or.nelmt_fs.eq.0)then
  nnode_fs=0
  ! No need to plot free surface files
  savedata%fsplot=.false.
  savedata%fsplot_plane=.false.
  return
endif

allocate(iface_fs(nelmt_fs), id_elem_fs(nelmt_fs), fs_elem_id(nelmt) )
allocate(gnum4_fs(4,nelmt_fs),gnum_fs(maxngll2d,nelmt_fs))
nsnode_all=nelmt_fs*maxngll2d
allocate(nodelist(nsnode_all),inode_order(nsnode_all))


fs_elem_id = 0 

n1=1; n2=maxngll2d
do i_face=1,nelmt_fs
  read(11,*)ielmt,iface

  iface_fs(i_face)=iface
  num=g_num(:,ielmt)
  gnum4_fs(:,i_face)=num(hexface(iface)%gnode)
  gnum_fs(:,i_face)=num(hexface(iface)%node)
  
  id_elem_fs(i_face) = ielmt  ! Input ielmt_fs --> returns global elemnt ID
  fs_elem_id(ielmt) = i_face  ! Input global element ID --> returns ielmt_fs
  nodelist(n1:n2)=num(hexface(iface)%node)
  n1=n2+1; n2=n1+maxngll2d-1
enddo

close(11)

! Renumber connectivity for the surface elements
call i_uniinv(nodelist,inode_order)

nnode_fs=maxval(inode_order)
allocate(isnode(nnode_fs))
isnode=.false.

! Store global node IDs for free surface nodes
allocate(gnode_fs(nnode_fs))
gnode_fs(inode_order(1))=nodelist(1)
isnode(inode_order(1))=.true.
do i=2,nsnode_all
  if(.not.isnode(inode_order(i)))then
     isnode(inode_order(i))=.true.
     gnode_fs(inode_order(i))=nodelist(i)
  endif
enddo
deallocate(isnode,nodelist)
! Store the renumbered connectivity for the free surface elements.
n1=1; n2=maxngll2d
allocate(rgnum_fs(maxngll2d,nelmt_fs))
do i_face=1,nelmt_fs
 rgnum_fs(:,i_face)=inode_order(n1:n2)
 n1=n2+1; n2=n1+maxngll2d-1
enddo

deallocate(inode_order)

return

end subroutine prepare_free_surface
!===============================================================================

! Determine the elevation on the free surface for a given point.
subroutine free_surface_elevation(xp,elevation,ilocated)
use math_constants,only:ZERO
use math_library,only:IsPointInPolygon
use shape_library,only:shape_function_quad4p
use map_location,only:map_point2naturalquad4
use global,only:g_coord
implicit none
real(kind=kreal),intent(in) :: xp(2)
real(kind=kreal),intent(out) :: elevation
integer,intent(out) :: ilocated

integer :: i_face,iface
integer :: niter
real(kind=kreal) :: coord(2,4),vx(4),vy(4),vz(4)
real(kind=kreal) :: located_x(2),xip(2),errx
real(kind=kreal) :: shape_quad4(4)
logical :: isinside

elevation=ZERO
ilocated=0

iface=0
do i_face=1,nelmt_fs
  vx=g_coord(1,gnum4_fs(:,i_face))
  vy=g_coord(2,gnum4_fs(:,i_face))
  coord(1,:)=vx
  coord(2,:)=vy
  vz=g_coord(3,gnum4_fs(:,i_face))
  isinside=IsPointInPolygon(vx,vy,xp(1),xp(2))
  if(isinside)then
    iface=iface+1
    call  map_point2naturalquad4(coord,xp,xip,located_x,niter,errx)
    call shape_function_quad4p(4,xip(1),xip(2),shape_quad4)
    elevation=sum(vz*shape_quad4)
    
    ilocated=1
    return
    ! NOTE: we could run for all surface elements and 
    ! check whether it is actually located only on the single element!
  endif
enddo
return
end subroutine free_surface_elevation
!===============================================================================

subroutine cleanup_free_surface
implicit none
if(allocated(gnum4_fs))deallocate(gnum4_fs)
if(allocated(gnum_fs))deallocate(gnum_fs)
if(allocated(gnode_fs))deallocate(gnode_fs)
if(allocated(rgnum_fs))deallocate(rgnum_fs)
end subroutine cleanup_free_surface
!===============================================================================



subroutine check_surface_normals()
  use global
  use set_precision
  use integration 
  use math_constants
  use element,only:hexface,hexface_sign,hex8_gnode

  implicit none 

  integer                        :: i_elmtfs, i_gll ! loops
  real(kind=kreal)               :: detjac 
  integer                        :: iface           ! face ID for elmt 
  integer                        :: i_elmt, errcode          ! face ID for elmt 
  integer                        :: nfgll           ! ngll on 2D face
  real(kind=kreal), allocatable  :: gw(:), nodalsl(:) ! GLL weights 2D
  real(kind=kreal), allocatable  :: dshape4(:,:,:)
  real(kind=kreal)               :: coord(ndim,4), face_normal(3),dx_dxi(NDIM), dx_deta(NDIM), vertical(3), dot_w_vert
        
  integer :: num4(4), counter


  ! Code
  allocate(gw(maxngll2d))
  allocate(dshape4(2,4,maxngll2d))


  if(myrank.eq.0)then
    write(*,*)'* checking surface normal'
  endif 

  vertical    = zero
  vertical(3) = one

  counter = 0 

  do i_elmtfs = 1, nelmt_fs
    ! Face number (ie between 1 and 6) and get related properties
    iface = iface_fs(i_elmtfs)    
    if(iface==1 .or. iface==3)then
        nfgll             = ngllzx
        gw(1:nfgll)       = gll_weights_zx
        dshape4(:,:,1:nfgll) = dshape_quad4_zx

      elseif(iface==2 .or. iface==4)then
        nfgll             = ngllyz
        gw(1:nfgll)       = gll_weights_yz
        dshape4(:,:,1:nfgll) = dshape_quad4_yz

      elseif(iface==5 .or. iface==6)then
        nfgll             = ngllzx
        gw(1:nfgll)       = gll_weights_xy
        dshape4(:,:,1:nfgll) = dshape_quad4_xy
      else
        return
    endif
    
    num4  = gnum4_fs(:, i_elmtfs)
    coord = g_coord(:,num4)

    do i_gll = 1, nfgll 
        ! Calculate the magnitude of the 2D jacobian 
        dx_dxi  = matmul(coord,dshape4(1,:,i_gll))
        dx_deta = matmul(coord,dshape4(2,:,i_gll))

        ! Calc normal and therefore jac dec (2D) on the fly
        face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
        face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
        face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)


        detjac=sqrt(dot_product(face_normal,face_normal))

        if (detjac.eq.zero)then 
          write(*,*)'ERROR! 2D jacobian magnitude is 0'
          write(*,*)'ielmtfs     = ', i_elmtfs
          write(*,*)'iface       = ', iface
          write(*,*)'i_gll       = ', i_gll
          write(*,*)'nfgll       = ', nfgll
          write(*,*)'dshape4 1   = ', coord,dshape4(1,:,i_gll)
          write(*,*)'dshape4 2   = ', coord,dshape4(2,:,i_gll)
          write(*,*)'num4        = ', num4
          write(*,*)'coord       = ', coord
          write(*,*)'dxi         = ', dx_dxi
          write(*,*)'dx_deta     = ', dx_deta
          write(*,*)'face_normal = ', face_normal
          write(*,*)'--------------------- '
          stop
        endif 

        face_normal=hexface_sign(iface)*face_normal/detjac
        dot_w_vert = dot_product(face_normal, vertical)

        ! Check dot with vertical is positive (or zero if orthogonal)
        if (dot_w_vert.lt.zero)then 
          write(*,*)'Error: negative normals on the free surface (pointing into mesh)' 
          write(*,*)'STOPPING...' 
          stop
        endif 

    enddo 
  enddo 


  if(myrank.eq.0)then
    write(*,*)' ✓ Done'
  endif

  deallocate(gw)
  deallocate(dshape4)
end subroutine check_surface_normals




end module free_surface
!===============================================================================
