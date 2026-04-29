module station
use set_precision
use global,only:myrank
implicit none
integer,allocatable :: station_id(:),station_myrank(:)
integer,allocatable :: station_funit(:)
integer,allocatable :: station_elmt(:)
real(kind=kreal),allocatable :: station_u(:,:)
real(kind=kreal),allocatable :: station_x(:,:)
real(kind=kreal),allocatable :: station_lagrange(:,:)
real(kind=kreal),allocatable :: station_deriv(:,:,:)
logical,allocatable :: station_islocated(:)
contains
!-------------------------------------------------------------------------------

! DESCRIPTION
!  This routine locates the stations and store their coordinates in a natural
!  coordinates system.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
!  Leah Langer, Princeton University
! REVISION
!  HNG, Feb 19, 2016; HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO
subroutine locate_station(errcode,errtag)
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
use elastic,only:compute_cmat_elastic
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

! maximum number of lines in the station file
integer,parameter :: nmax_line=100
integer :: i_elmt,i_line,nline,ielmt,imid,inum,ios,istat,istation
integer :: num(nenode),egdofu(nedofu)

real(kind=kreal) :: coord(ndim,8),jac(ndim,ndim),xp(ndim),xip(ndim)
real(kind=kreal) :: detjac

integer :: i_w,i_l,i_p,i_rec,i0,i1,i2,i3,i4,i0_nextrow
real(kind=kreal) :: located_x(ndim),rec_x(ndim),rec_xi(ndim)
real(kind=kreal) :: deriv(ndim,ngll)

integer :: mdomain

logical :: is_located
integer :: this_rec_located,total_rec_located
integer,allocatable:: irec_located(:)
integer :: niter
real(kind=kreal) :: errx,errxd,minerr,maxerr,minerrg,maxerrg

integer :: r_elmt,n_relmt

integer :: id
real(kind=kreal) :: x(NDIM)

real(kind=kreal),dimension(:,:),allocatable :: dshape_hex8
real(kind=kreal),dimension(:),allocatable :: lagrange_gll
real(kind=kreal),dimension(:,:),allocatable :: dlagrange_gll

integer :: i_gll,ielmt_min,inode_min
integer :: nelmt_rectry
integer,allocatable :: ielmt_rectry(:)
logical,allocatable :: isnode(:),iselmt(:)
real(kind=kreal) :: minsqdist,sqdist

integer :: ielmt_rec
integer :: indmin(1),rec_inrank
real(kind=kreal) :: gminerr,gmaxerr
real(kind=kreal) :: gminerr_rec,gmaxerr_rec
real(kind=kreal) :: all_minerr(1,0:nproc-1)
integer :: ipass_strict,nfail_strict
logical :: isinside

character(len=1) :: tchar
character(len=250) :: line,tag
character(len=80) :: token
character(len=80) :: fname
character(len=80) :: data_path
character(len=250) :: pfile

errtag="ERROR: unknown!"
errcode=-1

! set data path
data_path=trim(inp_path)

! station files
fname=trim(data_path)//trim(stationfile)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif

nstation=0
nline=0
! count the number of stations
do i_line=1,nmax_line
  ! This will read a line and proceed to next line
  read(11,'(a)',iostat=ios)line
  if (ios/=0)exit
  nline=nline+1
  ! check for blank and comment line
  if (isblank(line) .or. iscomment(line,'#'))cycle
  nstation=nstation+1
enddo
close(11)
if(myrank==0)then
  write(logunit,'(a,i0)')'Total number of stations: ',nstation
  flush(logunit)
endif
allocate(station_id(nstation),station_myrank(nstation),station_elmt(nstation), &
station_x(NDIM,nstation))
allocate(station_funit(nstation))
station_funit=-9999
allocate(station_u(ndim,nstation))
allocate(station_islocated(nstation))
station_islocated=.false.
allocate(station_lagrange(ngll,nstation),station_deriv(NDIM,ngll,nstation))
istation=0
! reopen station file
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
do i_line=1,nline
  ! This will read a line and proceed to next line
  read(11,'(a)',iostat=ios)line
  if (ios/=0)exit
  ! check for blank and comment line
  if (isblank(line) .or. iscomment(line,'#'))cycle
  istation=istation+1
  read(line,*)id,x
  station_id(istation)=id
  station_x(:,istation)=x
enddo
close(11)
! Nondimensionalize the station coordinates
station_x=NONDIM_L*station_x
if(myrank==0)then
  ! plot VTK file: stations
  pfile=trim(out_path)//trim(file_head)//'_station'//'.vtk'
  open(100,file=trim(pfile),action='write',status='replace')

  write(100,'(a)')'# vtk DataFile Version 2.0'
  write(100,'(a)')'Unstructured Grid Example'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(100,'(a,i5,a)')'POINTS',nstation,' float'
  do i_p=1,nstation
    write(100,'(3(e14.6,1x))')DIM_L*station_x(:,i_p)
  enddo
  write(100,*)
  write(100,'(a,i5,a,i5)')'CELLS',nstation,' ',2*nstation
  do i_p=1,nstation
    ! VTK format indexing starts from 0
    write(100,'(5(i5,1x))')1,i_p-1
  enddo
  write(100,*)
  write(100,'(a,i5)')'CELL_TYPES',nstation
  do i_p=1,nstation
    write(100,'(i1)')1
  enddo
  close(100)
  flush(100)
endif

! Get derivatives of shape functions for 8-noded hex.
allocate(dshape_hex8(ngnode,ndim))
allocate(lagrange_gll(ngll),dlagrange_gll(ndim,ngll))

imid=(ngll+1)/2

allocate(irec_located(nstation))
irec_located=0
! Compute contribution of the earthquake source/s to the elemental loads.
allocate(isnode(nnode),iselmt(nelmt))
gminerr_rec=INFTOL
gmaxerr_rec=ZERO
nfail_strict=0
station: do i_rec=1,1!nstation
  rec_x=station_x(:,i_rec)
  ! Find the element which contains this earthquake source.
  is_located=.false.
  
  ! Check if the source is within the model range.
  prange:if(rec_x(1).lt.NONDIM_L*pmodel_minx .or. rec_x(1).gt.NONDIM_L*pmodel_maxx .or. & 
     rec_x(2).lt.NONDIM_L*pmodel_miny .or. rec_x(2).gt.NONDIM_L*pmodel_maxy .or. & 
     rec_x(3).lt.NONDIM_L*pmodel_minz .or. rec_x(3).gt.NONDIM_L*pmodel_maxz)then
    nelmt_rectry=0
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
      if(sqdistance(rec_x,xp,ndim).gt.sqmaxsize_elmt)cycle

      do i_gll=1,ngll
        xp=g_coord(:,num(i_gll))
        sqdist=sqdistance(rec_x,xp,ndim)
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
    nelmt_rectry=1
    do i_elmt=1,nelmt
      if(i_elmt.eq.ielmt_min)cycle
      if(any(isnode(g_num(hex8_gnode,i_elmt))))then
        nelmt_rectry=nelmt_rectry+1
        iselmt(i_elmt)=.true.
      endif
    enddo
    ! Find and assign the list of elements to try.
    allocate(ielmt_rectry(nelmt_rectry))
    ielmt_rectry=-9999
    inum=0
    do i_elmt=1,nelmt
      if(iselmt(i_elmt))then
        inum=inum+1
        ielmt_rectry(inum)=i_elmt
      endif
    enddo
    if(inum.ne.nelmt_rectry)then
      write(*,*)'ERROR: nelmt_rectry mismatch!'
      stop
    endif
  endif prange
  n_relmt=0
  minerr=INFTOL
  ! Find the actual element that contains the source.
  ! Loop through the list of elements
  ! STRICT TEST: Element must be a proper hexagon, i.e., must have six faces. 
  ipass_strict=0
  element_try1: do i_elmt=1,nelmt_rectry
    ielmt=ielmt_rectry(i_elmt)

    num=g_num(:,ielmt)
    coord=g_coord(:,num(hex8_gnode))
    ! Check if the source is in this element.
    call IsPointInHexahedron(coord,rec_x,isinside)
    if(isinside)then
      n_relmt=n_relmt+1
      ! Map source location to natural coordinates.
      call  map_point2naturalhex8(coord,rec_x,xip,located_x,niter,errx)
      errxd=errx*DIM_L
      if(n_relmt==1)then
        ! Initialize errors
        minerr=errxd
        ! Set source element 
        ielmt_rec=ielmt
      elseif(n_relmt.gt.1)then
        if(errxd.lt.minerr)then
          minerr=errxd 
          ! Reset source element 
          ielmt_rec=ielmt
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
    element_try2: do i_elmt=1,nelmt_rectry
      ielmt=ielmt_rectry(i_elmt)

      num=g_num(:,ielmt)
      coord=g_coord(:,num(hex8_gnode))
      n_relmt=n_relmt+1
      ! map source location to natural coordinates
      call  map_point2naturalhex8(coord,rec_x,xip,located_x,niter,errx)
      errxd=errx*DIM_L
      if(n_relmt==1)then
        ! initialize errors
        minerr=errxd
        ! Set source element 
        ielmt_rec=ielmt
      elseif(n_relmt.gt.1)then
        if(errxd.lt.minerr)then
          minerr=errxd 
          ! Reset source element 
          ielmt_rec=ielmt
        endif
      endif
    enddo element_try2 ! i_elmt
  endif !(ipass_strict.lt.1)
  if(allocated(ielmt_rectry))deallocate(ielmt_rectry)
  ! Only the element with the least error was chosen
  n_relmt=sumscal(n_relmt)
  if(n_relmt.lt.1)then
    if(myrank==0)then
      write(logunit,'(a,3(g0.6,1x))')'WARNING: station cannot be located: ',i_rec,DIM_L*rec_x
      flush(logunit)
    endif
    cycle station
  else 
    irec_located(i_rec)=1
  endif
  gminerr=minscal(minerr)

  if(gminerr.lt.gminerr_rec)gminerr_rec=gminerr
  if(gminerr.gt.gmaxerr_rec)gmaxerr_rec=gminerr

  if(myrank==0)then
    if(gminerr.gt.DIM_L*maxsize_elmt)then
      write(logunit,'(a)')'WARNING: this station location may be inaccurate!'
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
  rec_inrank=indmin(1)-1
  ! Compute source contribution in the appropriate rank
  if(myrank==rec_inrank)then
    station_myrank(i_rec)=myrank
    station_elmt(i_rec)=ielmt_rec
    station_islocated(i_rec)=.true.
    rec_xi=xip
    is_located=.true.
    irec_located(i_rec)=1
    ! compute load
    ! get deriative of shape function at source location
    call dshape_function_hex8p(ngnode,rec_xi(1),rec_xi(2), &
                               rec_xi(3),dshape_hex8)

    ! compute jacobian
    !WAARNING: need to double check jacobian
    jac=matmul(transpose(dshape_hex8),transpose(coord))
    !jac=matmul(coord,dshape_hex8)
    detjac=determinant(jac)
    call invert(jac)

    ! compute GLL lagrange function and their derivative at source location
    call gll_lagrange3d_point(ndim,ngllx,nglly,ngllz,ngll,rec_xi(1), &
                              rec_xi(2),rec_xi(3),       &
                              lagrange_gll,dlagrange_gll)

    deriv=matmul(jac,dlagrange_gll) ! use der for gll
    station_lagrange(:,i_rec)=lagrange_gll
    station_deriv(:,:,i_rec)=deriv

  endif !(myrank==rec_inrank)
  
  call sync_process
enddo station ! i_rec
deallocate(isnode,iselmt)

call sync_process
irec_located=maxvec(irec_located,nstation)
where(irec_located.gt.1)irec_located=1
total_rec_located=sum(irec_located)
if(myrank==0)then
  write(logunit,'(a,i0,1x,a,i0)')'Total defined stations: ',nstation, &
                           'Total located stations: ',total_rec_located
  write(logunit,'(a,i0)')'Total stations fail strict test: ',nfail_strict
  if(nfail_strict.gt.0)then
    write(logunit,'(a,i0)')'NOTE: failed strict test indicates that some of &
    &the elements may NOT be proper hexahedra!'
  endif
  flush(logunit)
endif
if(total_rec_located.gt.0)then
  if(myrank==0)then
    write(logunit,'(a,g0.6,a,g0.6)')'station location errors (m) => min: ', &
    gminerr_rec, ' max: ',gmaxerr_rec
    flush(logunit)
  endif
endif
deallocate(irec_located)
deallocate(dshape_hex8)
deallocate(lagrange_gll,dlagrange_gll)

deallocate(station_x)

errcode=0

return

end subroutine locate_station
!===============================================================================

subroutine open_station_files
use global,only:out_path,nstation
implicit none
integer :: i_rec
character(len=255) :: fname

do i_rec=1,nstation
 if(station_myrank(i_rec)==myrank)then
    if(station_islocated(i_rec))then
      write(fname,'(a,i0,a)')trim(out_path)//'Station',station_id(i_rec),'.semd'
      open(newunit=station_funit(i_rec),file=trim(fname))
    endif
  endif
enddo

end subroutine open_station_files
!===============================================================================

subroutine close_station_files
use global,only:nstation
implicit none

integer :: i_rec

do i_rec=1,nstation
  if(station_myrank(i_rec)==myrank)then
    if(station_islocated(i_rec))then
      close(station_funit(i_rec))
    endif
  endif
enddo

end subroutine close_station_files
!===============================================================================

subroutine write_station_files(step)
use global,only:nstation
use nondimensionpar,only:DIM_L
implicit none
real(kind=kreal),intent(in) :: step
integer :: i_rec

do i_rec=1,nstation
  if(station_myrank(i_rec)==myrank)then
    if(station_islocated(i_rec))then
      write(station_funit(i_rec),'(e12.6,1x,e12.6,1x,e12.6,1x,e12.6)')step,DIM_L*station_u(:,i_rec)
      flush(station_funit(i_rec))
    endif
  endif
enddo

end subroutine write_station_files
!===============================================================================

subroutine compute_station(errcode,errtag)
use nondimensionpar
use global
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_rec,ielmt
integer :: num(nenode),egdofu(nedofu)

real(kind=kreal) :: ugllx(ngll),uglly(ngll),ugllz(ngll)

errtag="ERROR: unknown!"
errcode=-1

! number of stations
station: do i_rec=1,nstation
  if(station_myrank(i_rec)==myrank)then
    if(station_islocated(i_rec))then
      ielmt=station_elmt(i_rec)
      num=g_num(:,ielmt)
      ugllx=nodalu(1,num)
      uglly=nodalu(2,num)
      ugllz=nodalu(3,num)
 
      station_u(1,i_rec)=dot_product(station_lagrange(:,i_rec),ugllx)
      station_u(2,i_rec)=dot_product(station_lagrange(:,i_rec),uglly)
      station_u(3,i_rec)=dot_product(station_lagrange(:,i_rec),ugllz)
    endif
  endif
enddo station ! i_rec

errcode=0

return

end subroutine compute_station
!===============================================================================

end module station
!===============================================================================
