! subroutine to calculate the prestress in the mesh 
! Last edit - WE 8th June 2022

module prestress

contains 
    

subroutine calculate_prestress(strain_elmt, strain_nodal, &
                                stress_elmt, stress_nodal, & 
                                extload, du, dprecon, storekmat, &
                                isgravity, ispseudoeq, &
                                errcode, errtag, ksp_iter, istat)
! USES 
use global
use preprocess
use math_constants
use output_to_user
#if (USE_MPI)
use mpi_library
#else
use serial_library
#endif
use element

implicit none
! IO Variables
 
real(kind=kreal), allocatable :: strain_elmt(:,:,:), &
                                    strain_nodal(:,:),  &
                                    stress_elmt(:,:,:), & 
                                    stress_nodal(:,:), & 
                                    extload(:), & 
                                    du(:), dprecon(:) , storekmat(:,:,:)
logical :: isgravity, ispseudoeq
integer :: ksp_iter

character(len=250) :: errtag ! error message
integer :: errcode
! Local Variables

integer :: istat


if(savedata%stress.or.isplastic)then
    allocate(stress_elmt(nst,ngll,nelmt),stress_nodal(nst,nnode))
    stress_elmt=ZERO
endif

if(savedata%strain)then
    allocate(strain_elmt(nst,ngll,nelmt),strain_nodal(nst,nnode))
    strain_elmt=ZERO
endif

if(isstress0)then
    if(s0_type==0)then

        ! compute initial stress using SEM itself
        allocate(extload(0:neq),du(0:neq),dprecon(0:neq), &
        storekmat(nedof,nedof,nelmt),stat=istat)
        if (istat/=0)then
            write(logunit,*)'ERROR: cannot allocate memory!'
            flush(logunit)
            stop
        endif

        extload=ZERO; isgravity=.true.; ispseudoeq=.false.
        call stiffness_bodyload(nelmt,neq,hex8_gnode,g_num,gdof_elmt,mat_id,gam_blk, &
        storekmat,dprecon,extload,isgravity,ispseudoeq)

        log_msg = 'complete...' ; call write_ifproc0()
        log_msg = '--------------------------------------------' ; call write_ifproc0()

        ! assemble from ghost partitions
        if(nproc.gt.1)then
            call assemble_ghosts(nndof,neq,dprecon,dprecon)
        endif

        dprecon(1:)=one/dprecon(1:); dprecon(0)=ZERO

        ! compute displacement due to graviy loading to compute initial stress
        du=ZERO
        call ksp_pcg_solver(neq,nelmt,storekmat,du,extload,   &
        dprecon,gdof_elmt,ksp_iter,errcode,errtag)
        call control_error(errcode,errtag,stdout,myrank)

            du(0)=ZERO

        call elastic_stress(nelmt,neq,hex8_gnode,g_num,gdof_elmt,du,stress_elmt)
            deallocate(extload,dprecon,du,storekmat)

    elseif(s0_type==1)then
        ! compute initial stress using simple relation for overburden pressure
            call overburden_stress(nelmt,g_num,z_datum,s0_datum,epk0,stress_elmt)
        else
        write(logunit,*)'ERROR: s0_type:',s0_type,' not supported!'
        flush(logunit)
        stop
    endif
endif

return 


end subroutine calculate_prestress
end module prestress