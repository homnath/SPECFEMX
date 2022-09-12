! DESCRIPTION
!  This module conains the routines to set and analyze the degrees of freedoms.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!  HNG, Feb 19, 2016; HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO
!  - implement power law
module dof
implicit none
character(len=250),private :: myfname=' => degrees_of_freedom.f90'
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

! This subroutine sets number of degrees of freedom and IDs of nodal dof  
subroutine initialize_dof()
use global
implicit none
integer :: i_dof,idof, nelmt_fs, ios
character(len=80) :: data_path
character(len=80) :: fname
integer :: errcode
character(len=250) :: errtag



! total number of degrees of freedom per node
nndof    = 0
nedofu   = 0  
nedofphi = 0 
nedof    = 0

! initialise dof IDs
idof    = 0
idofu   = 0
idofphi = 0


! displacement
if(ISDISP_DOF)then
  nndof=nndof+nndofu 
  nedofu=NNDOFU*nenode 
  nedof=nedof+nedofu

  do i_dof=1,nndofu
    
    idof=idof+1
    idofu(i_dof)=idof

  enddo

  allocate(edofu(nedofu))
endif


! gravity

idof=idofu(nndofu)

if(ISPOT_DOF)then
  nndof=nndof+nndofphi
  nedofphi=NNDOFPHI*nenode
  nedof=nedof+nedofphi

  do i_dof=1,nndofphi
    idof=idof+1

    idofphi(i_dof)=idof

  enddo
  allocate(edofphi(nedofphi))
endif



if(ISSL_DOF)then 
  nndof = nndof + 1 

  ! Assert that all DOFs are on: 
  if (nndof/=5)then 
    write(*,*)'ERROR: Trying to calculate SL but not 5 DOFs switched on'
    stop 
  endif

  ! Note here we are setting a possible DOF for theta on every GLL in an 
  ! element...of course it can only actually be on the surface but 
  ! for now it is convenient ... may be a memory issue if mesh is huge
  ! For things like kmat
  nedof = nedof + nenode  
          
endif 




end subroutine initialize_dof
!===============================================================================



! This subroutine sets IDs for the elemental degrees of freedom for
! u and \phi which may be used to map the elemental matrices
subroutine set_element_dof()
use global,only:ISDISP_DOF,ISPOT_DOF,ISSL_DOF, nedofu,nedofphi,ngll,nndofu, &
                edofu,edofphi,edofsl
implicit none
integer :: i,iu(NNDOFU),iphi,j,nu, isl
integer :: iu0, iphi0, isl0


! order ux,uy,uz,\phi
edofu=-9999
if(ISPOT_DOF)edofphi=-9999
if(ISSL_DOF)edofsl=-9999

iu0=0
iphi0=0
isl0=0

iphi=0
isl=0
nu=0
iu=0
do i=1,NGLL

  if(ISDISP_DOF)then
    iu(1)=iu0+1
    nu=nu+1
    edofu(nu)=iu(1)


    do j=2,NNDOFU
      nu=nu+1
      iu(j)=iu(j-1)+1
      edofu(nu)=iu(j)

    enddo

    iu0=iu(NNDOFU) ! this will be overwritten if POT_DOF is present
    iphi0=iu(NNDOFU)
  endif


  if(ISPOT_DOF)then
    iphi=iphi0+1
    edofphi(i)=iphi

    iu0=iphi ! this will be overwritten if DISP_DOF is present
    iphi0=iphi
  endif

enddo



 ! if(ISSL_DOF)then
 !   isl=isl0+1
 !   edofsl(i)=isl

 !   iu0=isl    ! Not sure when this will be overwritetn!
 !   isl0=isl
 ! endif



return
end subroutine set_element_dof
!===============================================================================

! compute mapping of u  to face elemental matrices
subroutine set_face_vecdof(nfgll,ncomp,nxtra,imapuf)                           
implicit none                                                                    
integer,intent(in) :: nfgll,ncomp,nxtra                                          
integer,intent(out) :: imapuf(ncomp*nfgll)                                       
integer :: i,j,iu(ncomp),nu                                                       
                                                                                 
! u must be in the first order, i.e., order ux,uy,uz,phi                                 
nu=0; iu(1)=1                                                                    
do i=1,nfgll                                                                     
  nu=nu+1                                                                        
  imapuf(nu)=iu(1)                                                               
  do j=2,ncomp                                                                   
    nu=nu+1                                                                      
    iu(j)=iu(j-1)+1                                                              
    imapuf(nu)=iu(j)                                                             
  enddo                                                                          
  iu(1)=iu(ncomp)+nxtra+1
enddo
return
end subroutine set_face_vecdof
!===============================================================================                                     
                                                                                 
! compute mapping of scalar, i.e., phi to face elemental matrices
subroutine set_face_scaldof(nfgll,npos,nxtra,imapsf)                           
implicit none                                                                    
integer,intent(in) :: nfgll,npos,nxtra                                           
integer,intent(out) :: imapsf(nfgll)                                             
integer :: i,is,nu                                                               
                                                                                 
! u must be in first order, i.e., order ux,uy,uz,phi                                 
nu=0; is=npos                                                                    
do i=1,nfgll                                                                     
  nu=nu+1                                                                        
  imapsf(nu)=is                                                                  
  is=is+npos+nxtra                                                               
enddo                                                                            
return                                                                           
end subroutine set_face_scaldof                                                
!===============================================================================

! This subroutine activates degrees of freedoms.
subroutine activate_dof(errcode,errtag)
use global
use free_surface
use math_constants,only:ZERO
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: ios, ctr 
integer :: i_elmt,imat,mdomain, i_face
integer :: inodes(ngll)
integer :: numf(maxngll2d)

errtag="ERROR: unknown!"
errcode=-1

! Initialize all DOFs to OFF
gdof=0

! Prescribed DOFs will be made OFF in apply_bc routine.
! Displacement
if(ISDISP_DOF)then
  do i_elmt=1,nelmt
    inodes=g_num(:,i_elmt)
    imat=mat_id(i_elmt)
    mdomain=mat_domain(imat)
    ! elastic domain
    if(mdomain==ELASTIC_DOMAIN)then
      gdof(idofu,inodes)=1
    ! viselastic domain
    elseif(mdomain==VISCOELASTIC_DOMAIN)then
      gdof(idofu,inodes)=1
    ! acoustic domain
    elseif(mdomain==ACOUSTIC_DOMAIN)then
      write(errtag,'(a)')'ERROR: acoustic domain not supported!'
      return
    ! trasnition infinite/infinite domain: only gravity 
    elseif(mdomain==ELASTIC_TRINFDOMAIN .or. &
           mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
           mdomain==ELASTIC_INFDOMAIN .or. &
           mdomain==VISCOELASTIC_INFDOMAIN)then
      if(infbc)then
        if(.not.isempty_blk(imat))then
        ! This will also have displacement DOF
          gdof(idofu,inodes)=1
        endif
      endif
    else
      write(errtag,'(a)')'ERROR: unsupported material domain!'
      return
    endif
  enddo
  !gdof(idofu,:)=1
endif

! Gravity
if(ISPOT_DOF)then
  ! gravity exists everywhere
  gdof(idofphi,:)=1
endif

! Sea Level
if(ISSL_DOF)then
  ! only for nodes on the free surface
  do i_face = 1, nelmt_fs  
    ! Get g_num values of this face 
    numf  = gnum_fs(:,i_face)
    ! Set these nodes for the index idofsl to 1 (activate them) 
    ! Note here that for SL to be solved we need displacement and 
    ! phi to be solved so theta will always be the 5th dof slot 
    gdof(5, numf) = 1
  enddo 
endif

errcode=0

end subroutine activate_dof
!===============================================================================

! This subroutine finalizes the global degrees of freedom IDs.
subroutine finalize_gdof(errcode,errtag)
use global, only:gdof,neq,nndof,nnode,g_num,part_path,proc_str,file_head
use global, only:myrank,nedof, ISSL_DOF, SLlogunit
use free_surface
implicit none
integer,intent(out) :: errcode
character(len=250) :: ofname
character(len=250),intent(out) :: errtag
integer :: i,istat,j, neqsl

errtag="ERROR: unknown!"
errcode=-1
! Compute modified gdof
neq=0
neqsl = 0

write(SLlogunit,*)
write(SLlogunit,*)'Finalising DOFs: '

if (ISSL_DOF)then 
  ! If SL then we want to run it as original version first by ignoring 
  ! the theta, and then tag on the theta DOF after 

  do j=1,ubound(gdof,2)
    do i=1,ubound(gdof,1)-1 ! -1 so that not including theta
      if(gdof(i,j)/=0)then
        neq=neq+1
        gdof(i,j)=neq
      endif
    enddo
  enddo

  write(SLlogunit,*)' - Number of eq. for U, Phi: ', neq

  ! Now index the SL ones 
  do j=1,ubound(gdof,2)
      if(gdof(5,j)/=0)then ! always 5th dof 
        neq=neq+1
        gdof(5,j)=neq
        neqsl = neqsl + 1 
      endif
  enddo


  write(SLlogunit,*)' - Number of eq. for theta : ', neqsl
  write(SLlogunit,*)'   .......................................'
  write(SLlogunit,*)' - Total number of eqns    : ', neq
  write(SLlogunit,*)'   .......................................'
  write(SLlogunit,*)

else
  ! Original version 
  do j=1,ubound(gdof,2)
    do i=1,ubound(gdof,1)
      if(gdof(i,j)/=0)then
        neq=neq+1
        gdof(i,j)=neq
      endif
    enddo
  enddo

endif 




ofname='tmp/'//trim(file_head)//'_gdof'//trim(adjustl(proc_str))
open(unit=22,file=trim(ofname),access='stream',form='unformatted', &
status='replace',action='write',iostat=istat)
if (istat /= 0)then
  write(errtag,'(a)')'ERROR: output file "'//trim(ofname)//'" cannot be opened!'
  stop
endif
write(22)nnode
write(22)neq
write(22)gdof
close(22)
ofname='tmp/'//trim(file_head)//'_gnum'//trim(adjustl(proc_str))
open(unit=22,file=trim(ofname),access='stream',form='unformatted', &
status='replace',action='write',iostat=istat)
if (istat /= 0)then
  write(errtag,'(a)')'ERROR: output file "'//trim(ofname)//'" cannot be opened!'
  stop
endif
write(22)nnode
write(22)g_num
close(22)

! Compute nodal to global
errcode=0
return
end subroutine finalize_gdof
!===============================================================================

subroutine sea_level_dof()
  ! DOF setup for the sea level stuff - requires free surface to be run first.
  use global 
  use free_surface 
  implicit none 

  integer :: ldof, i 

  allocate(edofsl(nnode_fs))

  ! Get last degree of freedom from phi + 1: 
  ldof = edofphi(nedofphi) + 1

  ! DOFs for SL are sequence starting with ldof since after u and phi
  do i = 1, nnode_fs
    edofsl(i) = ldof
    ldof = ldof + 1 
  enddo 


end subroutine sea_level_dof






end module dof
!===============================================================================
