module electrical
contains
!-------------------------------------------------------------------------------

! DESCRIPTION
!  This routine computes the load contributed by the electrical moment-tensor
!  source. The moement tensor is read either directly from CMTSOLUTION or
!  computed from the slip information file.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
!  Leah Langer, Princeton University
! REVISION
!  HNG, Feb 19, 2016; HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO
!  - read and implement CMTSOLUTION
!  - check for the sources shared among elements/processors
subroutine electrical_load(errcode,errtag)
use nondimensionpar
use global
use element,only:hex8_gnode
use math_constants
use conversion_constants
use math_library,only:cross_product,determinant,sqdistance,invert,  &
IsPointInHexahedron,norm,vector_rotateZ
use string_library
use shape_library,only : dshape_function_hex8p
use gll_library,only : gll_lagrange3d_point,zwgljd
use map_location
#if (USE_MPI)
use mpi_library
use math_library_mpi
#else
use serial_library
use math_library_serial
#endif
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer,parameter :: nmax_line=100 ! maximum number of lines in the eqsource
!file
integer :: i_elmt,i_line,ielmt,imid,inum,ios,istat,ix
integer :: num(nenode),egdofphi(nedofphi)

real(kind=kreal) :: coord(ndim,8),jac(ndim,ndim),xp(ndim),xip(ndim)
real(kind=kreal) :: detjac

real(kind=kreal) :: proj
real(kind=kreal) :: located_x(ndim),source_x(ndim),source_xi(ndim)
real(kind=kreal) :: deriv(ndim,ngll),eload(nedofphi)

integer :: mdomain

logical :: is_located
integer :: this_src_located,total_src_located
integer,allocatable:: isrc_located(:)
integer :: niter
real(kind=kreal) :: errx,errxd,minerr

integer :: npoint_ecurrent
real(kind=kreal),allocatable :: ecurrent_coord(:,:),ecurrent(:)

integer :: f_elmt,n_felmt

real(kind=kreal),dimension(:,:),allocatable :: dshape_hex8
real(kind=kreal),dimension(:),allocatable :: lagrange_gll
real(kind=kreal),dimension(:,:),allocatable :: dlagrange_gll

real(kind=kreal) :: sload(0:neq)

integer :: i_gll,ielmt_min,inode_min,i_point,i_src
integer :: nelmt_srctry
integer,allocatable :: ielmt_srctry(:)
logical,allocatable :: isnode(:),iselmt(:)
real(kind=kreal) :: minsqdist,sqdist

integer :: ielmt_src
integer :: indmin(1),src_inrank
real(kind=kreal) :: gminerr,gmaxerr
real(kind=kreal) :: gminerr_src,gmaxerr_src
real(kind=kreal) :: all_minerr(1,0:nproc-1)
integer :: ipass_strict,nfail_strict
logical :: isinside

character(len=1) :: tchar
character(len=80) :: token
character(len=80) :: fname
character(len=80) :: data_path
character(len=250) :: pfile

errtag="ERROR: unknown!"
errcode=-1

! set data path
data_path=trim(inp_path)

! Read electrical current
! open electrical current file
fname=trim(data_path)//trim(ecfile)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif

read(11,*) ! Skip a line
read(11,*)npoint_ecurrent
allocate(ecurrent_coord(3,npoint_ecurrent),ecurrent(npoint_ecurrent))
read(11,*) ! Skip a line
do i_point=1,npoint_ecurrent
  read(11,*)ecurrent_coord(:,i_point),ecurrent(i_point)
enddo
close(11)

nsource=npoint_ecurrent

! Get derivatives of shape functions for 8-noded hex.
allocate(dshape_hex8(ngnode,ndim))
allocate(lagrange_gll(ngll),dlagrange_gll(ndim,ngll))

imid=(ngll+1)/2

allocate(isrc_located(nsource))
isrc_located=0
! Compute contribution of the electrical source/s to the elemental loads.
allocate(isnode(nnode),iselmt(nelmt))
gminerr_src=INFTOL
gmaxerr_src=ZERO
nfail_strict=0
source: do i_src=1,nsource
  sload=zero
  source_x=ecurrent_coord(:,i_src)
  ! Find the element which contains this electrical source.
  is_located=.false.
  
  print*,myrank, source_x(1),pmodel_minx,pmodel_maxx
  print*,myrank, source_x(2),pmodel_miny,pmodel_maxy
  print*,myrank, source_x(3),pmodel_minz,pmodel_maxz

  ! Check if the source is within the model range.
  prange:if(source_x(1).lt.NONDIM_L*pmodel_minx .or. source_x(1).gt.NONDIM_L*pmodel_maxx .or. & 
     source_x(2).lt.NONDIM_L*pmodel_miny .or. source_x(2).gt.NONDIM_L*pmodel_maxy .or. & 
     source_x(3).lt.NONDIM_L*pmodel_minz .or. source_x(3).gt.NONDIM_L*pmodel_maxz)then
    nelmt_srctry=0
    ! NOTE: we cannot cycle the loop here otherwise the process hangs for ever
    ! because inbetween there are some statements which require all procesors!
  else
 
    ! First, find the element which is not too far and has the nearest GLL point
    ! from the target point.
    ielmt_min=1
    inode_min=g_num(1,ielmt_min)
    minsqdist=INFTOL
    do i_elmt=1,nelmt
      ! Do not locate sources in transition or infinite elements.
      mdomain=mat_domain(mat_id(i_elmt))
      if(mdomain==ELASTIC_TRINFDOMAIN .or.      &
         mdomain==ELASTIC_INFDOMAIN   .or.      &
         mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
         mdomain==VISCOELASTIC_INFDOMAIN)cycle
      ! Preliminary test. Discard far enough elements.
      num=g_num(:,i_elmt)
      xp=g_coord(:,num(imid))
      if(sqdistance(source_x,xp,ndim).gt.sqmaxsize_elmt)cycle

      do i_gll=1,ngll
        xp=g_coord(:,num(i_gll))
        sqdist=sqdistance(source_x,xp,ndim)
        if(sqdist.lt.minsqdist)then
          minsqdist=sqdist
          ielmt_min=i_elmt
          inode_min=num(i_gll)
        endif
      enddo
    enddo

    ! Count the number of and tag elements that share the GLL point just found.
    isnode=.false.
    isnode(g_num(hex8_gnode,ielmt_min))=.true.
    iselmt=.false.
    iselmt(ielmt_min)=.true.
    nelmt_srctry=1
    do i_elmt=1,nelmt
      if(i_elmt.eq.ielmt_min)cycle
      if(any(isnode(g_num(hex8_gnode,i_elmt))))then
        nelmt_srctry=nelmt_srctry+1
        iselmt(i_elmt)=.true.
      endif
    enddo
    ! Find and assign the list of elements to try.
    allocate(ielmt_srctry(nelmt_srctry))
    ielmt_srctry=-9999
    inum=0
    do i_elmt=1,nelmt
      if(iselmt(i_elmt))then
        inum=inum+1
        ielmt_srctry(inum)=i_elmt
      endif
    enddo
    if(inum.ne.nelmt_srctry)then
      write(*,*)'ERROR: nelmt_srctry mismatch!'
      stop
    endif
  endif prange
  !print*,'try elements:',i_src,myrank,nelmt_srctry,count(iselmt) 

  n_felmt=0
  minerr=INFTOL
  ! Find the actual element that contains the source.
  ! Loop through the list of elements
  ! STRICT TEST: Element must be a proper hexagon, i.e., must have six faces. 
  ipass_strict=0
  element_try1: do i_elmt=1,nelmt_srctry
    ielmt=ielmt_srctry(i_elmt)

    num=g_num(:,ielmt)
    coord=g_coord(:,num(hex8_gnode))
    ! Check if the source is in this element.
    call IsPointInHexahedron(coord,source_x,isinside)
    if(isinside)then
      n_felmt=n_felmt+1
      ! Map source location to natural coordinates.
      call  map_point2naturalhex8(coord,source_x,xip,located_x,niter,errx)
      errxd=errx*DIM_L
      if(n_felmt==1)then
        ! Initialize errors
        minerr=errxd
        ! Set source element 
        ielmt_src=ielmt
      elseif(n_felmt.gt.1)then
        if(errxd.lt.minerr)then
          minerr=errxd 
          ! Reset source element 
          ielmt_src=ielmt
        endif
      endif
      ipass_strict=1
    endif !(isinside)
  enddo element_try1 ! i_elmt
  
  ipass_strict=sumscal(ipass_strict)
  ! FLEXIBLE TEST: Element may be a proper hexagon, i.e., may have more than 
  !                six faces. 
  ! There will be a WARNING in the log file. 
  ! This test will be performed only if the STRICT TEST fails.
  if(ipass_strict.lt.1)then
  ! Source point was not located on the element list
    nfail_strict=nfail_strict+1
    element_try2: do i_elmt=1,nelmt_srctry
      ielmt=ielmt_srctry(i_elmt)

      num=g_num(:,ielmt)
      coord=g_coord(:,num(hex8_gnode))
      n_felmt=n_felmt+1
      ! map source location to natural coordinates
      call  map_point2naturalhex8(coord,source_x,xip,located_x,niter,errx)
      errxd=errx*DIM_L
      if(n_felmt==1)then
        ! initialize errors
        minerr=errxd
        ! Set source element 
        ielmt_src=ielmt
      elseif(n_felmt.gt.1)then
        if(errxd.lt.minerr)then
          minerr=errxd 
          ! Reset source element 
          ielmt_src=ielmt
        endif
      endif
    enddo element_try2 ! i_elmt
  endif !(ipass_strict.lt.1)
  if(allocated(ielmt_srctry))deallocate(ielmt_srctry)
  ! Only the element with the least error was chosen
  n_felmt=sumscal(n_felmt)
  if(n_felmt.lt.1)then
    if(myrank==0)then
      write(logunit,'(a,3(g0.6,1x))')'WARNING: source cannot be located: ',i_src,DIM_L*source_x
      flush(logunit)
    endif
    cycle source
  else 
    isrc_located(i_src)=1
  endif
  gminerr=minscal(minerr)

  if(gminerr.lt.gminerr_src)gminerr_src=gminerr
  if(gminerr.gt.gmaxerr_src)gmaxerr_src=gminerr

  if(myrank==0)then
    if(gminerr.gt.DIM_L*maxsize_elmt)then
      write(logunit,'(a)')'WARNING: this source location may be inaccurate!'
      write(logunit,'(3(g0.6,1x))')xip
      write(logunit,'(2(g0.6,1x))')gminerr,maxsize_elmt
      flush(logunit)
    endif
  endif

  ! Determine the processor this source belongs to.
  all_minerr=allgather_scal(minerr)
  !if(myrank==0)write(*,'(g0)')all_minerr
  indmin=minloc(all_minerr,2)
  ! NOTE: minloc gives the position starting from 1
  src_inrank=indmin(1)-1
  !print*,i_src,src_inrank
  ! Compute source contribution in the appropriate rank
  if(myrank==src_inrank)then
    !print*,'myrank:',myrank,i_src,minerr,gminerr
    f_elmt=ielmt_src
    source_xi=xip
    is_located=.true.
    isrc_located(i_src)=1
    ! compute load
    ! get deriative of shape function at source location
    call dshape_function_hex8p(ngnode,source_xi(1),source_xi(2), &
                               source_xi(3),dshape_hex8)

    ! compute jacobian
    !WAARNING: need to double check jacobian
    jac=matmul(transpose(dshape_hex8),transpose(coord))
    detjac=determinant(jac)
    call invert(jac)

    ! compute GLL lagrange function and their derivative at source location
    call gll_lagrange3d_point(ndim,ngllx,nglly,ngllz,ngll,source_xi(1), &
                              source_xi(2),source_xi(3),       &
                              lagrange_gll,dlagrange_gll)

    egdofphi=gdof_elmt(edofphi,f_elmt)

    ! compute electrical load
    eload=ecurrent(i_src)*lagrange_gll
    sload(egdofphi)=sload(egdofphi)+eload
    !exit element ! this will place the source in only one element
  endif !(myrank==src_inrank)
  
  call sync_process
  ! NOTE: We allow only a single element to have the source point. 
  ! Therefore, no averaging is necessary.
  ! However, this may have to be modified for the instance where
  ! the source point lies in the interface of two different material domains!
  !n_felmt=sumscal(n_felmt)
  !if(n_felmt.gt.1)then
  !  ! add average load per source
  !  load=load+sload/real(n_felmt,kreal)
  !else
    load=load+sload
  !endif
enddo source ! i_src
deallocate(isnode,iselmt)

call sync_process
isrc_located=maxvec(isrc_located,nsource)
where(isrc_located.gt.1)isrc_located=1
total_src_located=sum(isrc_located)
if(myrank==0)then
  write(logunit,'(a,i0,1x,a,i0)')'Total defined sources: ',nsource, &
                           'Total located sources: ',total_src_located
  write(logunit,'(a,i0)')'Total sources fail strict test: ',nfail_strict
  if(nfail_strict.gt.0)then
    write(logunit,'(a,i0)')'NOTE: failed strict test indicates that some of &
    &the elements may NOT be proper hexahedra!'
  endif
  flush(logunit)
endif
if(total_src_located.gt.0)then
  if(myrank==0)then
    write(logunit,'(a,g0.6,a,g0.6)')'source location errors (m) => min: ', &
    gminerr_src, ' max: ',gmaxerr_src
    flush(logunit)
  endif
endif
deallocate(isrc_located)
deallocate(dshape_hex8)
deallocate(lagrange_gll,dlagrange_gll)

deallocate(ecurrent_coord)

errcode=0

return

end subroutine electrical_load
!===============================================================================

end module electrical
!===============================================================================
