module sea_level 
    use set_precision
    implicit none 

    contains 





subroutine start_SL_log(errcode, errtag)

use global
#if(USE_MPI)
use mpi_library
#else 
use serial_library
#endif  
    implicit none 

    character(len=250) :: errtag
    integer :: ios, errcode


    ! Create file 
    if(myrank==0)then
        SL_log_file = trim(file_head)//'SL.log'
        open(unit=SLlogunit,file=trim(SL_log_file),status='replace',action='write',iostat=ios)
        if(ios.ne.0)then
            write(errtag,'(a)')'ERROR: cannot open log file: '//trim(SL_log_file)
            call control_error(errcode,errtag,stdout,myrank)
        endif
    endif 

    ! Write example 
    write(SLlogunit, *)' ****** CREATED SL LOG FILE ******'
    write(*, *)' Created SL log file'

end subroutine start_SL_log
        




subroutine prepare_sea_level()
    use global 
    use free_surface

    implicit none 

    write(SLlogunit, *)'Preparing sea level variables...'

    allocate(oceanf(nelmt_fs, maxngll2d)) ! Allocate ocean function 
    


end subroutine prepare_sea_level






subroutine calc_SL_LHS(errcode, errtag)
use global 
use free_surface

use integration,only:dshape_quad4_xy,dshape_quad4_yz,dshape_quad4_zx,          &
                     gll_weights_xy,gll_weights_yz,gll_weights_zx,             &
                     lagrange_gll_xy,lagrange_gll_yz,lagrange_gll_zx
use math_library,only : angle

#if(USE_MPI)
use mpi_library
#else 
use serial_library
#endif 

    implicit none 

    character(len=250) :: errtag
    integer :: ios, errcode

    integer :: abar, bbar, gamma, x, y, i_face, iface, i_gll, j, gid,  &
               numf(maxngll2d)
    real(kind=kreal) :: coord(3, maxngll2d), G_cal, theta_tf, Gj_sum, & 
                        u_tf(3), phi_tf

    integer :: nfdof,nfgll !face nodal dof, face gll pts
    real(kind=kreal),dimension(:,:,:),allocatable :: dshape_quad4
    real(kind=kreal),dimension(:),allocatable :: gll_weights
    real(kind=kreal),dimension(:,:),allocatable :: lagrange_gll
    real(kind=kreal),dimension(:,:,:),allocatable :: dlagrange_gll



    allocate(dshape_quad4(2,4,maxngll2d))
    allocate(gll_weights(maxngll2d),lagrange_gll(maxngll2d,maxngll2d),             &
    dlagrange_gll(2,maxngll2d,maxngll2d))


    ! Calculating R matrices as defined in BF latex doc 

    
    ! Need to loop through each element face combination for the solid
    ! surface 
    do i_face=1, nelmt_fs  

        ! Face number (ie between 1 and 6) and get related properties
        iface = iface_fs(i_face)    
        if(iface==1 .or. iface==3)then
            nfgll=ngllzx
            gll_weights(1:nfgll)=gll_weights_zx
            dshape_quad4(:,:,1:nfgll)=dshape_quad4_zx
          elseif(iface==2 .or. iface==4)then
            nfgll=ngllyz
            gll_weights(1:nfgll)=gll_weights_yz
            dshape_quad4(:,:,1:nfgll)=dshape_quad4_yz
          elseif(iface==5 .or. iface==6)then
            nfgll=ngllzx
            gll_weights(1:nfgll)=gll_weights_xy
            dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
          else
            write(errtag,'(a)')'ERROR: wrong face ID for traction!'
            return
        endif

        ! Get coordinates for the face and global IDs 
        numf  = gnum_fs(:,i_face)
        coord = g_coord(:, numf)


        ! Loop through GLL on the surface: 
        do i_gll = 1,nfgll
            
            ! global ID for this GLL point 
            gid = numf(i_gll) 

            !==================== CALCULATE G ==========================
            ! Initialise test functions in case testing code 
            theta_tf = 1.0_kreal 
            phi_tf = 1.0_kreal 
            u_tf = 1.0 

            ! J sum as part of G_cal calculation 
            !Gj_sum = 0.0 
            !do j =1,3
            !    Gj_sum = Gj_sum + u_tf(j)*g0_local(3*(gid-1) + j) 
            !enddo 

            !G_cal = g0(gid)*theta_tf  + oceanf(gid)*(phi_tf + Gj_sum)
            ! ===============  FINISHED CALC G  ====================
    

            ! Because R and R_cal are diagonal matrices it will be better
            ! to store them as vectors 
    
        enddo 

       

    enddo 
    end subroutine calc_SL_LHS






    subroutine update_ocean_function(u, errcode, errtag)
    ! Routine checks each GLL point on the surface to see if its theta
    ! Value is 0 or not - if it is 0 then the ocean function becomes 0 
    ! if SL is larger than 0 then ocean function is 1.
    use global 
    use free_surface
    implicit none 

    ! IO variables
    character(len=250) :: errtag
    integer :: ios, errcode
    real(kind=kreal), allocatable :: u(:)

    ! Local variables 
    integer          :: i_face, iface, numf(maxngll2d), i_numf
    real(kind=kreal) :: sl

    

    ! Loop through the free surface faces: 
    do i_face=1, nelmt_fs  

        ! Get the global ID of the nodes on this face 
        numf  = gnum_fs(:,i_face)

        ! Loop through nodes on this face
        do i_numf = 1, maxngll2d

            if (u(gdof(idofsl(1), numf(i_numf))).gt.0.0 ) then 
                ! Sea level is not zero - ocean func is 1 
                oceanf(i_face, i_numf) = 1.0_kreal

            elseif(u(gdof(idofsl(1), numf(i_numf)))==0.0 ) then 
                ! Sea level is zero - ocean func is 0 
                oceanf(i_face, i_numf) = 0.0_kreal

            else 
                write(errtag,'(a)')'SEA LEVEL VALUE IS NEGATIVE!! '
                return
            endif 
        enddo 

        
    enddo

    end subroutine update_ocean_function






end module