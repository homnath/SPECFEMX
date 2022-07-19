! This module contains routines to create spectral elements from the eight-noded
! hexahedral elements
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module mesh_spec
use set_precision
private :: rank, get_global,get_global_indirect_addressing, swap_all

contains
!-------------------------------------------------------------------------------

! This subroutine just runs a few commands to keep the driver file clean
subroutine create_spec_elem(tot_nelmt,max_nelmt,min_nelmt, &
                            tot_nnode,max_nnode,min_nnode, &
                            errcode,errtag)

! USES
use global
use output_to_user
#if (USE_MPI)
use mpi_library
use math_library_mpi
#else
use serial_library
use math_library_serial
#endif
use element
implicit none 
! IO variables
character(len=250) :: errtag ! error message
integer :: errcode
integer :: tot_nelmt,max_nelmt,min_nelmt,tot_nnode,max_nnode,min_nnode

! Local variables

! Code: 
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

  return

end subroutine create_spec_elem








! This subroutine convert all hexahedral meshes (8-noded) to spectral elements
! of arbitrary order defined by ngllx, nglly, and ngllz
subroutine hex2spec(ndim,ngnod,nelmt,nnode,ngllx,nglly,ngllz,errcode,errtag)
use global,only : g_coord,g_num, logunit
use shape_library,only : shape_function_hex8
use gll_library,only:gllpx,gllpy,gllpz

implicit none
integer,intent(in) :: ndim,ngnod,ngllx,nglly,ngllz
integer,intent(inout) :: nelmt,nnode
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: ngll
!double precision
real(kind=kreal),parameter :: jacobi_alpha=0.0d0,jacobi_beta=0.0d0,zero=0.0d0
integer :: i,i_elmt,i_gnod,j,k
real(kind=kreal),dimension(:),allocatable :: xstore,ystore,zstore
real(kind=kreal),dimension(ngnod,ngllx,nglly,ngllz) :: shape_hex8

real(kind=kreal) :: xgll,ygll,zgll
real(kind=kreal) :: xmin,xmax

integer :: ipoint,npoint
integer :: istat
integer :: ienode,inode

integer, dimension(:), allocatable :: iglob

errtag="ERROR: unknown!"
errcode=-1

write(logunit,*)'Running hex2spec...'

ngll=ngllx*nglly*ngllz
xmin=minval(g_coord(1,:))
xmax=maxval(g_coord(1,:))
npoint=nelmt*(ngllx*nglly*ngllz)

write(logunit,*)'   ngllx:', ngllx
write(logunit,*)'   nglly:', nglly
write(logunit,*)'   ngllz:', ngllz
write(logunit,*)'   ngll (ngllx * nglly * ngllz):', ngll
write(logunit,*)'   npoint (ngll * npoint): ', npoint


allocate(xstore(npoint),ystore(npoint),zstore(npoint),stat=istat)
if(istat/=0)then
  write(errtag,'(a)')'ERROR: cannot allocate memory!'
  return
else
  write(logunit,*)'   Created xstore, ystore, zstore each of length :', npoint
endif



! get shape function for 8-noded hex
call shape_function_hex8(ngnod,ngllx,nglly,ngllz,gllpx,gllpy,gllpz,shape_hex8)


write(logunit,*)' Filling in xstore, ystore, zstore with coordinates of local gll points:'
! compute coordinates all local gll points
xstore=zero
ystore=zero
zstore=zero
write(logunit,*)''
ipoint=0
do i_elmt=1,nelmt
  write(logunit,*)' element = ', i_elmt
  do k=1,ngllz
    do j=1,nglly
      do i=1,ngllx
        !write(logunit,*)'   i: ', i, 'j: ', j, 'k: ', k

        xgll = zero
        ygll = zero
        zgll = zero
        !write(logunit,*)'   reset xgll, ygll, zgll to 0'

       ! write(logunit,*)'looping from i_gnod = 1 to ngnod: ', ngnod
        do i_gnod=1,ngnod
          xgll = xgll + shape_hex8(i_gnod,i,j,k)*g_coord(1,g_num(i_gnod,i_elmt))
          ygll = ygll + shape_hex8(i_gnod,i,j,k)*g_coord(2,g_num(i_gnod,i_elmt))
          zgll = zgll + shape_hex8(i_gnod,i,j,k)*g_coord(3,g_num(i_gnod,i_elmt))

          !write(logunit,*)'   xgll += ', shape_hex8(i_gnod,i,j,k), ' * ', g_coord(1,g_num(i_gnod,i_elmt))
        enddo

        ipoint=ipoint+1
        xstore(ipoint) = xgll
        ystore(ipoint) = ygll
        zstore(ipoint) = zgll

      enddo
    enddo
  enddo
enddo




write(logunit,*)' Finished looping through all elements'

deallocate(g_coord,g_num) ! no longer need these
allocate(iglob(npoint))

! gets ibool indexing from local (gll points) to global points
call get_global(ndim,xstore,ystore,zstore,iglob,nnode,npoint,xmin,xmax)

write(*,*)'nnode:  ', nnode

write(*,*)'XSTORE: '
do i =1, npoint 
  write(*,*) xstore(i)
enddo 

write(*,*)'YSTORE: '
do i =1, npoint 
  write(*,*) ystore(i)
enddo 

write(*,*)'ZSTORE: '
do i =1, npoint 
  write(*,*) zstore(i)
enddo 

write(*,*)'IGLOB: '
do i =1, npoint 
  write(*,*) iglob(i)
enddo 


write(*,*)'Now using indirect addressing'
! now we got the new number of nodes (nnode)
!- we can create a new indirect addressing to reduce cache misses
call get_global_indirect_addressing(nnode,npoint,iglob)

write(*,*)'IGLOB is now : '
do i =1, npoint 
  write(*,*) iglob(i)
enddo 


allocate(g_coord(3,nnode),g_num(ngll,nelmt)) ! alocate with new number of nodes
ipoint=0
do i_elmt=1,nelmt
  ienode=0
  do k=1,ngllz
    do j=1,nglly
      do i=1,ngllx
        ienode=ienode+1
        ipoint=ipoint+1
        !if(ifseg(ipoint))then
        inode=iglob(ipoint) !ibool(i,j,k,i_elmt) ! iglob(locval(ipoint))
        g_num(ienode,i_elmt)=inode
        g_coord(1,inode)=xstore(ipoint)
        g_coord(2,inode)=ystore(ipoint)
        g_coord(3,inode)=zstore(ipoint)
        !endif
      enddo
    enddo
  enddo
enddo

deallocate(iglob,xstore,ystore,zstore)

write(logunit, *)'Finished hex2spec'

errcode=0
return

end subroutine hex2spec
!============================================

subroutine get_global(ndim,xold,yold,zold,iglob,nnode,npoint,xmin,xmax)
! this routine must be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by paul f. fischer

! leave the sorting subroutines in the same source file to allow for inlining

implicit none
integer :: ndim,npoint,nnode
integer :: iglob(npoint)
integer :: iloc(npoint) ! at first this is the location (indices) of all points, at the end
!it becomes the orderd location (indices) according to the sorted rank
logical :: ifseg(npoint) ! initailly all are false, afterward only the unique
!ponts become true
real(kind=kreal),intent(in) :: xold(npoint),yold(npoint),zold(npoint) !double precision
real(kind=kreal) :: xp(npoint),yp(npoint),zp(npoint) ! double precision. initially these are the
!original coordinates of all points, at the end these become the coordinates of
!ordered points
real(kind=kreal) :: xmin,xmax

integer :: i,j,ier
integer :: nseg,ioff,iseg,ig

integer, dimension(:), allocatable :: ind,ninseg,iwork
real(kind=kreal), dimension(:), allocatable :: work !double precision

! geometry tolerance parameter to calculate number of independent grid points
! small value for double precision and to avoid sensitivity to roundoff
real(kind=kreal) :: smalltol !double precision

xp=xold
yp=yold
zp=zold

! define geometrical tolerance based upon typical size of the model
smalltol = 1.e-10_kreal * abs(xmax - xmin)
!write(*,*)'SMALL TOLERANCE IS:', smalltol


!write(*,*)'Entered get_global(): '

! dynamically allocate arrays
  allocate(ind(npoint), &
          ninseg(npoint), &
          iwork(npoint), &
          work(npoint),stat=ier)
  if( ier /= 0 )then
    write(*,*)'ERROR: error allocating arrays!'
    stop
  endif



! establish initial pointers
  do i=1,npoint
    iloc(i)=i
  enddo

  ifseg=.false.

  nseg=1
  ifseg(1)=.true.
  ninseg(1)=npoint

  !write(*,*)'Sorting: '



  do j=1,ndim
    write(*,*)'j: ',j



! sort within each segment
    ioff=1
    !write(*,*), 'nseg:', nseg

    do iseg=1,nseg

      !write(*,*)'  iseg = ',iseg

      if(j == 1) then
        call rank(xp(ioff), ind, ninseg(iseg))
      else if(j == 2) then
        !write(*,*)'   yp(ioff): ', yp(ioff)
        !write(*,*)'   ind: ', ind 
        !write(*,*)'   ninseg(iseg): ', ninseg(iseg) 
        !write(*,*)'   RUN RANK YP: '

        call rank(yp(ioff),ind,ninseg(iseg))
      else

        call rank(zp(ioff),ind,ninseg(iseg))
      endif

      call swap_all(iloc(ioff),xp(ioff),yp(ioff),zp(ioff),iwork,work,ind,ninseg(iseg))
      !write(*,*)' AFTER SWAP THE VALUES ARE NOW: '
      !write(*,*)'iloc: ', iloc
      !!write(*,*)'xp: ', xp
      !write(*,*)'yp: ', yp
      !write(*,*)'zp: ', zp
      !write(*,*)'i_work: ', iwork
      !!write(*,*)'ind: ', ind
      !write(*,*)'ninseg: ', ninseg
      
      ioff=ioff+ninseg(iseg)
      !write(*,*)'new ioff value is : ', ioff 
    enddo




! check for jumps in current coordinate
! compare the coordinates of the points within a small tolerance
    if(j == 1) then
      do i=2,npoint
        if(abs(xp(i)-xp(i-1)) > smalltol) ifseg(i)=.true.
      enddo
    else if(j == 2) then
      do i=2,npoint
        if(abs(yp(i)-yp(i-1)) > smalltol) ifseg(i)=.true.
      enddo
    else
      do i=2,npoint
        if(abs(zp(i)-zp(i-1)) > smalltol) ifseg(i)=.true.
      enddo
    endif


! count up number of different segments
    nseg=0
    
    do i=1,npoint
      if(ifseg(i)) then
        nseg=nseg+1
        ninseg(nseg)=1
      else
        ninseg(nseg)=ninseg(nseg)+1
      endif
    enddo



  enddo ! j=1,ndim



! assign global node numbers (now sorted lexicographically)
  ig=0
  do i=1,npoint
    if(ifseg(i)) ig=ig+1
    iglob(iloc(i))=ig
  enddo

  nnode=ig

  !write(*,*)'IGLOB: ', iglob

! deallocate arrays
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iwork)
  deallocate(work)

  end subroutine get_global
!===========================================




! sorting routines put in same file to allow for inlining
  subroutine rank(a,ind,n)
!
! use heap sort (numerical recipes)
!
  implicit none

  integer :: n
  real(kind=kreal) :: a(n) !double precision
  integer :: ind(n)

  integer :: i,j,l,ir,indx
  real(kind=kreal) :: q !double precision


  !write(*,*)'inside rank: '
  !write(*,*)'a = ', a
  !write(*,*)'n = ', n
  !write(*,*)'ind = ', ind



  do j=1,n
   ind(j)=j
  enddo
  !write(*,*)'ind array is now: ', ind 


  if (n == 1) return

  l=n/2+1
  ir=n

  !write(*,*)'l = ', l
  !write(*,*)'ir = ', ir


  100 continue

   !write(*,*) 'l=  ', l
   if (l>1) then
       !write(*,*)'l is larger than 1 so subtract 1'
      l=l-1
      indx=ind(l)
      !write(*,*)'indx = ind(l) = ', indx
      q=a(indx)
      !write(*,*)'q = a(indx) = ', q

   else
    !write(*,*)'l is NOT larger than 1'

      indx=ind(ir)
      !write(*,*)'indx = ind(ir) = ', indx 

      q=a(indx)
      !write(*,*)'q = a(indx) = ', q 

      ind(ir)=ind(1)
      !write(*,*)'ind(ir) = ind(1)', ind(ir) 

      ir=ir-1
      !write(*,*)'ir -= 1 so ir =', ir 

      if (ir == 1) then
        !write(*,*)'ir now == 1 so'

         ind(1)=indx
         !write(*,*)'ind(1) = indx and RETURN'
         return
      endif
   endif

  
   i=l
   !write(*,*)'i = l =', i 
   j=l+l
   !write(*,*)'j = l+l =', j 

   !write(*,*)'____ AT 200 pt'
  200    continue
   if (j <= ir) then
    !write(*,*)'j <= ir so enter'

      if (j<ir) then
        !write(*,*)'j < ir so enter'

         if ( a(ind(j))<a(ind(j+1)) ) then 
          !write(*,*)'a(ind(j))<a(ind(j+1)) so enter'
          j=j+1
          !write(*,*) 'j is now += 1   = ', j
         endif 
      endif

      if (q<a(ind(j))) then
        !write(*,*)'q<a(ind(j)) so enter'
         ind(i)=ind(j)
         !write(*,*)'ind(i) = ind(j) = ', ind(i)
         i=j
         j=j+j
         !write(*,*)'i = j = ', i
         !write(*,*)'j = j+j = ', j
      else
        !write(*,*)'NOT TRUE q<a(ind(j)) so enters else'

         j=ir+1
         !write(*,*)'j =ir+1 =  ', j
      endif
      !write(*,*)'Go to 200'
   goto 200
   endif
   ind(i)=indx
   !write(*,*)'ind(i) = indx =', ind(i)
   !write(*,*)''
   !write(*,*)''
   !write(*,*)'Go to 100'
   !write(*,*)''
   !write(*,*)''

  goto 100

end subroutine rank
!===========================================

subroutine swap_all(ia,a,b,c,iw,w,ind,n)
!
! swap arrays ia, a, b and c according to addressing in array ind
!
  implicit none

  integer :: n

  integer :: ind(n)
  integer :: ia(n),iw(n)
  real(kind=kreal) :: a(n),b(n),c(n),w(n) !double precision

  integer :: i


  write(*,*) " INSIDE SWAPPING SUBROUTINR"

  iw(:) = ia(:)
  w(:) = a(:)

  do i=1,n
    ia(i)=iw(ind(i))
    a(i)=w(ind(i))
  enddo

  w(:) = b(:)

  do i=1,n
    b(i)=w(ind(i))
  enddo

  w(:) = c(:)

  do i=1,n
    c(i)=w(ind(i))
  enddo

end subroutine swap_all
!===========================================

subroutine get_global_indirect_addressing(nnode,npoint,ibool)
!- we can create a new indirect addressing to reduce cache misses
! (put into this subroutine but compiler keeps on complaining that it can't vectorize loops...)

implicit none

integer :: nnode,npoint
integer, dimension(npoint) :: ibool

! mask to sort ibool
integer, dimension(nnode) :: mask_ibool
integer, dimension(npoint) :: copy_ibool_ori
integer :: inumber
integer:: i_point

mask_ibool = -1
copy_ibool_ori = ibool
! reduces misses

inumber = 0
do i_point=1,npoint
  write(*,*)'i_point = ', i_point
  write(*,*)'      copy_ibool_ori(i_point) = ', copy_ibool_ori(i_point) 
  write(*,*)'mask_ibool(copy_ibool_ori(i_point)) = ', mask_ibool(copy_ibool_ori(i_point))

  if(mask_ibool(copy_ibool_ori(i_point)) == -1) then
    write(*,*)'mask_ibool is -1 so enter if...'
    inumber = inumber + 1
    write(*,*)'inumber becomes = ', inumber 
    ibool(i_point) = inumber
    mask_ibool(copy_ibool_ori(i_point)) = inumber
    write(*,*)'make  = mask_ibool(copy_ibool_ori(i_point)) = ', inumber 

  else
    write(*,*)'mask_ibool is NOT -1 so else...'
    ! use an existing point created previously
    ibool(i_point) = mask_ibool(copy_ibool_ori(i_point))
    write(*,*)'ibool(i_point) becomes ', ibool(i_point)

  endif
enddo
return
end subroutine get_global_indirect_addressing
!===========================================

end module mesh_spec
