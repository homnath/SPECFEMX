! count_elements.f90
! Last edit: WE Jun 2 2022

module count_elements 

contains 
!_______________________________________________________________________________
subroutine count_elmts( errcode ) 

  ! Subroutine used to loop through the elements and determine how many
  ! are viscoelastic vs elastic
  
  
use global ! uses mat_domain, nelmt mat_id and DOMAIN variables
use math_library_mpi
implicit none

! Local variables: 
integer :: i_elmt, i_mat     ! Loop index for elements    
integer  :: mdomain
integer :: tot_nelmt_elas
integer :: max_nelmt_elas
integer :: min_nelmt_elas
integer :: tot_nelmt_viscoelas
integer :: max_nelmt_viscoelas
integer :: min_nelmt_viscoelas
integer :: errcode


! count elastic and viscoelastic elements
nelmt_elas=0; nelmt_viscoelas=0
do i_elmt=1,nelmt
  mdomain=mat_domain(mat_id(i_elmt))
  ! infinite elements included in elastic domain
  if(mdomain==ELASTIC_DOMAIN .or.  &
    mdomain==ELASTIC_TRINFDOMAIN .or.  &
    mdomain==ELASTIC_INFDOMAIN)then
    nelmt_elas=nelmt_elas+1

  elseif(mdomain==VISCOELASTIC_DOMAIN .or. &
    mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
    mdomain==VISCOELASTIC_INFDOMAIN)then
    nelmt_viscoelas=nelmt_viscoelas+1

  else
    write(logunit,*)'ERROR: unrecognized material domain',mdomain,'!'
    flush(logunit)
    stop
  endif
enddo

if(nelmt/=nelmt_elas+nelmt_viscoelas)then
  write(logunit,*)'ERROR: total number of elements mismatch!'
  flush(logunit)
  stop
endif

tot_nelmt_elas=sumscal(nelmt_elas)
max_nelmt_elas=maxscal(nelmt_elas)
min_nelmt_elas=minscal(nelmt_elas)

tot_nelmt_viscoelas=sumscal(nelmt_viscoelas)
max_nelmt_viscoelas=maxscal(nelmt_viscoelas)
min_nelmt_viscoelas=minscal(nelmt_viscoelas)
if(myrank==0)then
  write(logunit,'(a,i0,1x,a,i0,1x,a,i0)')'elements elastic => total: ', &
  tot_nelmt_elas,' max: ',max_nelmt_elas,' min: ',min_nelmt_elas
  write(logunit,'(a,i0,1x,a,i0,1x,a,i0)')'elements viscoelastic => total: ', &
  tot_nelmt_viscoelas,' max: ',max_nelmt_viscoelas,' min: ',min_nelmt_viscoelas
  flush(logunit)
endif

errcode=0
return
end subroutine count_elmts
!-------------------------------------------------------------------------------

subroutine split_elas_visco_eids()
! Used to save element ID separately for elastic and viscoelastic 
! elements. Ids are stored in eid_viscoelas and eid_elas
use math_library_mpi
use global!,  only: ngll, nmaxwell, nelmt, mat_domain, mat_id, &
           !        ELASTIC_DOMAIN, ELASTIC_TRINFDOMAIN, &
           !        ELASTIC_INFDOMAIN, VISCOELASTIC_DOMAIN,  & 
            !       VISCOELASTIC_TRINFDOMAIN, VISCOELASTIC_INFDOMAIN

implicit none 
  
! Local variables
integer :: ielmt_elas, ielmt_viscoelas, i_elmt, mdomain


! Allocate arrays:
allocate(eid_elas(nelmt_elas),eid_viscoelas(nelmt_viscoelas))
ielmt_elas=0; ielmt_viscoelas=0

! Assign IDs:
do i_elmt=1,nelmt
  mdomain=mat_domain(mat_id(i_elmt))
  ! infinite elements included in elastic domain
  if(mdomain==ELASTIC_DOMAIN .or.  &
    mdomain==ELASTIC_TRINFDOMAIN .or.  &
    mdomain==ELASTIC_INFDOMAIN)then
    ielmt_elas=ielmt_elas+1
    eid_elas(ielmt_elas)=i_elmt
  elseif(mdomain==VISCOELASTIC_DOMAIN .or. &
    mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
    mdomain==VISCOELASTIC_INFDOMAIN)then
    ielmt_viscoelas=ielmt_viscoelas+1
    eid_viscoelas(ielmt_viscoelas)=i_elmt
  endif
enddo

return 
end subroutine split_elas_visco_eIDs
!-------------------------------------------------------------------------------

subroutine store_elemdof_from_nodaldof_global()

  use global, only: gdof_elmt, nedof, nelmt, gdof, g_num
  implicit none

  ! Local variables: 
  integer :: i_elmt

  allocate(gdof_elmt(nedof,nelmt))
  
  gdof_elmt=0 ! initialise values at 0 
  do i_elmt=1,nelmt
    gdof_elmt(:,i_elmt)=reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
  enddo

end subroutine
!-------------------------------------------------------------------------------

!subroutine calc_prestress(strain_elmt, strain_nodal, stress_elmt, &
!  stress_nodal, evpt, bodyload, slipload, &
!  extload, viscoload, ubcload, load, du,  & 
!  dprecon, kmat, storekmat, ksp_iter, & 
!  errcode, errtag)

!use global ! savedata, isplastic
!use math_constants
!use preprocess
!use element
!use elastic
!use output_to_user
!use mpi_library       ! control_error
!use ghost_library_mpi ! assemble_ghosts

!implicit none 
! IO variables
!real(kind=kreal),allocatable :: strain_elmt(:,:,:),                &
!          strain_nodal(:,:),                 &
!          stress_elmt(:,:,:),                & 
!          stress_nodal(:,:),                 &
!          evpt(:,:,:),  bodyload(:),         &
!          slipload(:),  extload(:),          &
!          viscoload(:), ubcload(:), load(:), &
!          du(:), dprecon(:), kmat(:,:),      &
!          storekmat(:,:,:)
!logical :: isgravity, ispseudoeq
!integer :: ksp_iter, errcode
!character(len=250) :: errtag ! error message


! Local variables
!integer :: istat



! compute initial stress assuming elastic domain
!------------------------------
!end subroutine calc_prestress

subroutine calculate_valency()
  use global
  use local
  implicit none 

  ! Local: 
  integer :: i_elmt, ielmt 

  allocate(node_valency(nnode))
  node_valency=0
  do i_elmt=1,nelmt
    ielmt=i_elmt
    num=g_num(:,ielmt)
    node_valency(num)=node_valency(num)+1
  enddo

  return 
end subroutine calculate_valency
!-------------------------------------------------------------------------------

end module count_elements
!===============================================================================
