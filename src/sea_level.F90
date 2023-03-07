module sea_level 
    use global 
    use set_precision
    implicit none 

    contains 

    

! Probably not worth using as it is SLOW
subroutine calc_jacdet2d(i_gll, coord, dshape_quad4, detjac)
    use global 
    use set_precision
    implicit none 
    real(kind=kreal) :: face_normal(NDIM)

    integer :: i_gll

    real(kind=kreal) :: coord(:,:), dshape_quad4(:,:,:), dx_dxi(NDIM), dx_deta(NDIM)
    real(kind=kreal) :: detjac

    dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
    dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))

    ! Normal
    face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
    face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
    face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

    
    detjac=sqrt(dot_product(face_normal,face_normal))
end subroutine  calc_jacdet2d





! ####################### FUNCS FOR GETTING FS DETAILS #########################
subroutine get_fs_details(i_elmt, iface, nfgll, gw, dsq4)
    use global 
    use set_precision
    use integration 
    use free_surface
    implicit none 

    ! IO variables: 
    integer           :: i_elmt, nfgll ,iface 
    real(kind=kreal)  :: gw(:), dsq4(:,:,:)


    ! Face number (ie between 1 and 6) and get related properties
    iface = iface_fs(i_elmt)    
    if(iface==1 .or. iface==3)then
        nfgll             = ngllzx
        gw(1:nfgll)       = gll_weights_zx
        dsq4(:,:,1:nfgll) = dshape_quad4_zx

      elseif(iface==2 .or. iface==4)then
        nfgll             = ngllyz
        gw(1:nfgll)       = gll_weights_yz
        dsq4(:,:,1:nfgll) = dshape_quad4_yz

      elseif(iface==5 .or. iface==6)then
        nfgll             = ngllzx
        gw(1:nfgll)       = gll_weights_xy
        dsq4(:,:,1:nfgll) = dshape_quad4_xy
      else
        !write(errtag,'(a)')'ERROR: wrong face ID for traction!'
        return
    endif
end subroutine get_fs_details




subroutine get_fs_details_nodsq(i_elmt, iface, nfgll, gw)
    use global 
    use set_precision
    use integration 
    use free_surface
    implicit none 

    ! IO variables: 
    integer           :: i_elmt, nfgll ,iface 
    real(kind=kreal)  :: gw(:)


    ! Face number (ie between 1 and 6) and get related properties
    iface = iface_fs(i_elmt)    
    if(iface==1 .or. iface==3)then
        nfgll             = ngllzx
        gw(1:nfgll)       = gll_weights_zx

      elseif(iface==2 .or. iface==4)then
        nfgll             = ngllyz
        gw(1:nfgll)       = gll_weights_yz

      elseif(iface==5 .or. iface==6)then
        nfgll             = ngllzx
        gw(1:nfgll)       = gll_weights_xy
      else
        !write(errtag,'(a)')'ERROR: wrong face ID for traction!'
        return
    endif
end subroutine get_fs_details_nodsq


subroutine get_fs_details_noweights(i_elmt, iface, nfgll)
    use global 
    use set_precision
    use integration 
    use free_surface
    implicit none 

    ! IO variables: 
    integer           :: i_elmt, nfgll ,iface 

    ! Face number (ie between 1 and 6) and get related properties
    iface = iface_fs(i_elmt)    
    if(iface==1 .or. iface==3)then
        nfgll             = ngllzx
      elseif(iface==2 .or. iface==4)then
        nfgll             = ngllyz
      elseif(iface==5 .or. iface==6)then
        nfgll             = ngllzx
      else
        !write(errtag,'(a)')'ERROR: wrong face ID for traction!'
        return
    endif
end subroutine get_fs_details_noweights
! #################### END GETTING FS DETAILS #########################






! ################### LOG AND OUTPUT FUNCTIONS  #######################
subroutine print_SL_read()

use global
#if(USE_MPI)
use mpi_library
#else 
use serial_library
#endif  
implicit none 

! Write summary of SL inputs: 
if(myrank.eq.0)then
    write(*,*)'-----------------------------------------------------'
    write(*,*)'Sea level data read from:  ', trim(slfile)
    write(*,*)'-----------------------------------------------------'
    write(*,*)'SL IS_SL                 :  ', IS_SL
    write(*,*)'SL DOF                   :  ', ISSL_DOF
    write(*,*)'Save original SL         :  ', savedata%sl0
    write(*,*)'Save SL                  :  ', savedata%sl
    write(*,*)'Save ocean func. initial :  ', savedata%oceanf0
    write(*,*)'Save ocean func.         :  ', savedata%oceanf
    write(*,*)'-----------------------------------------------------'
    write(*,*)


    if(IS_CART_SIM)then
        write(*,*)'Simulation type         :  CARTESIAN'
    elseif(IS_GLOB_SIM)then 
        write(*,*) 'ERROR: GLOBAL SIMULATIONS NOT IMPLEMETED YET'
        stop
    else
        write(*,*) 'ERROR: SIMULATION MUST BE GLOBAL OR CARTESIAN'
        stop
    endif 
endif 
end subroutine print_SL_read




subroutine summarise_SL_input(nodalsl)
    ! Summarises sea level data for cartesian sims
use global 
use dimensionless
use free_surface
#if(USE_MPI)
use mpi_library
use math_library_mpi
#else
use serial_library
use math_library_serial
#endif

    implicit none 

    real(kind=kreal) :: nodalsl(:)
    ! Local: 

    ! Collect total number of nodes on free surface and total ocean nodes
    totaloceannodes = sumscal(oceannodes) 


    if(myrank.eq.0)then 
        ! Model setup: 
        if(IS_CART_SIM)then 
            write(*,*)'*  Model setup                        :  Cartesian'
        else 
            write(*,*)'* Model setup                        :  Global'  
        endif 
        write(*,*)
        write(*,'(a,i6)')'* Total SL objects added        : ', nsl_obj 
        write(*,*)
    endif 

    call write_min_max_SL(nodalsl)

    if(myrank.eq.0)then 
        write(*,*)'-------------------------------------------------'
        write(*,*)
    endif 

     
end subroutine summarise_SL_input




subroutine write_min_max_SL(nodalsl)
    ! Outputs the min/max SL. Is also called by summarise_SL_input function
    use global
    use dimensionless
#if(USE_MPI)
use mpi_library
use math_library_mpi
#else
use serial_library
use math_library_serial
#endif
    implicit none

    ! IO variables: 
    real(kind=kreal) :: nodalsl(:)
    ! Local variables:
    real(kind=kreal) :: minnodalsl, maxnodalsl

    ! Minimum/maximum SL values: 
    minnodalsl = minscal(minval(nodalsl))
    maxnodalsl = maxscal(maxval(nodalsl))

    if(myrank.eq.0)then
        write(*,'(a,g0.6)')'  -->  Minimum sea level             : ', minnodalsl
        write(*,'(a,g0.6)')'  -->  Maximum sea level             : ', maxnodalsl
        if(devel_nondim)then
            write(*,'(a,g0.6)')'  -->  Dimensionalised min sea level : ', minnodalsl*DIM_L
            write(*,'(a,g0.6)')'  -->  Dimensionalised max sea level : ', maxnodalsl*DIM_L
        endif      
        write(*,*)   
    endif     
end subroutine


subroutine write_SL_to_ensight(nodalsl, i_step)
    use global 
    use postprocess
    use set_precision
    use dimensionless
    use free_surface
#if(USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif
implicit none 

    real(kind=kreal) :: nodalsl(:)
    integer :: i_step

    ! On the free surface
    if(savedata%fsplot)then
      call write_scalar_to_file_freesurf(nnode_fs, nodalsl*DIM_L, &
      ext='sl',istep=i_step) 
    endif
    
    if(savedata%fsplot_plane)then
      call write_scalar_to_file_freesurf(nnode_fs, nodalsl*DIM_L, &
      ext='sl', istep=i_step,plane=.true.) 
    endif

    if(myrank.eq.0)then
        write(*,'(a,i6)')'  ✓ Saved nodal sea level for step ', i_step
        write(*,*)
    endif 
    end subroutine write_SL_to_ensight



subroutine write_OF_to_ensight(i_step)
    ! Writes the ocean function to ensight 
    use global 
    use postprocess
    use free_surface
    implicit none 

    logical :: save_orig
    integer :: i_step


    ! On the free surface
    if(savedata%fsplot)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalOF, &
        ext='oceanf',istep=i_step) 
    endif
    
    if(savedata%fsplot_plane)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalOF, &
        ext='oceanf', istep=i_step, plane=.true.) 
    endif


    if(myrank.eq.0)then
        write(*,'(a,i6)')'  ✓ Saved ocean function for step ', i_step
        write(*,*)
    endif 
end subroutine write_OF_to_ensight


! ################# END  LOG AND OUTPUT FUNCTIONS  ####################




! #################    INITIAL SETUP FUNCTIONS    #####################
subroutine prepare_sea_level(nodalsl, nodalslrate)
    use global 
    use free_surface
    use set_precision
    use math_constants

    implicit none 

    integer :: istattemp, istat
    real(kind=kreal), allocatable :: nodalsl(:), nodalslrate(:)

    if(myrank.eq.0)then
        write(*,*)'Preparing sea level variables...'
    endif 

    istat = 0

    if(myrank.eq.0)then
        write(*,*)'  + number of unique FS nodes: ', nnode_fs
    endif 

    ! Allocate ocean function and nodal ocean func
    allocate(oceanf(nelmt_fs, maxngll2d), nodalOF(nnode_fs), stat=istattemp) 
    istat=istat+istattemp
    oceanf=ZERO 


    if(myrank.eq.0)then
        write(*,*)'  --> Created ocean function'
    endif 

    ! Allocate nodal sea level 
    allocate(nodalsl(nnode_fs), nodalslrate(nnode_fs), stat=istattemp)
    istat=istat+istattemp

    ! Check allocations 
    if(istat/=0)then
        write(*,*)'ERROR: cannot allocate memory in prepare_sea_level!'
        stop
    else 
        nodalsl = ZERO
        nodalslrate = ZERO

        if(myrank.eq.0)then
            write(*,*)'  --> Initialised nodalsl '
            write(*,*)'  --> Initialised nodalslrate '
        endif
    endif

    if(myrank.eq.0)then
        write(*,*)'  ✓ Prepared sea level.'
        write(*,*)
    endif 

    return 
end subroutine prepare_sea_level




subroutine set_original_sea_level(nodalsl)
    ! Uses
    use set_precision
    use global 
    use integration
    use dimensionless
    use free_surface
    use math_constants
    implicit none 

    ! IO variables
    real(kind=kreal) :: nodalsl(:)
    ! Local variables 
    integer          :: i_elmt, iface, nfgll, i_gll, i_obj
    integer          :: slobjtype
    real(kind=kreal) :: params(4)


    ! Code:
    if(myrank.eq.0)then
        write(*,*)'-----------------------------------------------------'
        write(*,*)'       Setting original water distribution          '
        write(*,*)
    endif


    nodalsl = 0.0_kreal

    do i_obj = 1, nsl_obj
        slobjtype = slobjs(i_obj, 1)

        ! Need to nondimensionalise the parameters: 
        params = slobjs(i_obj,2:5)*NONDIM_L

        if (slobjtype.eq.0) then 
            ! Single point of water - args: nodeid, height          
            call add_sl_gll(INT(params(1)), INT(params(2)), params(3),  INT(params(4)), nodalsl)
        
        elseif (slobjtype.eq.1) then 
        ! Uniform sea level z coordinate over all nodes
            if(IS_CART_SIM)then
                call set_cart_constant_SL0(nodalsl, params(1))
            else 
                write(*,*)'ERROR: UNIFORM SL Z ONLY FOR CARTESIAN SETUPS.'
                stop
            endif 
        else
            ! Invalid entry
            write(*,*)'ERROR: Unknown sea level object type: ', slobjtype
            stop
        endif 
    enddo 

end subroutine set_original_sea_level




subroutine set_cart_constant_SL0(nodalsl, sl_zcoord)
    use global
    use free_surface
    use dimensionless
    use set_precision
    ! IO: 
    real(kind=kreal) :: nodalsl(:), sl_zcoord

    ! Local: 
    integer :: i_elmt,i_gll, nfgll, iface
    real(kind=kreal) :: theta, z_coord

    ! SL0 is the z coordinate of the sea surface. We therefore
    ! need to calculate the theta value for each point based on the 
    ! the z coordinate of the face 

    ! Store value: 
    SL0_constant = sl_zcoord



    do i_elmt=1, nelmt_fs  
        call get_fs_details_noweights(i_elmt, iface, nfgll)

        do i_gll = 1, nfgll 
            ! Get Z coordinates for the face and global IDs 
            z_coord = g_coord(3,  gnum_fs(i_gll,i_elmt))
            ! SL = value - z_coord on surface
            theta   = sl_zcoord - z_coord 
            
            ! Now storing negative values too. 
            nodalsl(rgnum_fs(i_gll, i_elmt)) = theta    
        enddo 
    enddo 
    ! Log output
    if(myrank.eq.0)then 
        write(*,*)' -- Added water at constant Z value'
        write(*,*)'   --> value      : ', sl_zcoord
        if(devel_nondim)then 
            write(*,*)'   --> value (DIM): ', sl_zcoord*DIM_L
        endif 
        write(*,*)
    endif 
end subroutine set_cart_constant_SL0






subroutine add_sl_gll(i_elmtfs, i_gll, height, overwrite_int, nodalsl)
    ! Adds ice in the required location to a single GLL point 
    ! Uses
    use set_precision
    use global 
    use integration
    use free_surface
    use dimensionless
    use math_constants

    ! IO vars: 
    real(kind=kreal) :: height, nodalsl(:) 
    integer :: i_elmtfs, i_gll
    integer :: overwrite_int
    logical :: overwrite
    ! Params should be the faceID ON FS, nodeID, height, overwrite

    ! Process overwrite: 
    if (overwrite_int.eq.0.or.overwrite_int.eq.1) then 
        overwrite = overwrite_int
    else 
        write(*,*)'ERROR: OVERWRITE FLAG FOR SL GLL POINT MUST BE 1/0. Value given: ', overwrite_int
        stop 
    endif 

    if(myrank.eq.0)then
        write(*,*)
        write(*,*)' --  Adding sea level at point '
        write(*,'(a,i6)')'     -->  FS Elmt ID             : ', i_elmtfs
        write(*,'(a,i6)')'     -->  GLL Node (1-maxngll2d) : ', i_gll
        write(*,'(a,g0.6)')'     -->  height                 : ', height
        if(devel_nondim)then 
            write(*,'(a,g0.6)')'     --> non-dimensional height  : ', height*DIM_L
        endif
        write(*,'(a,g0.6)')'     -->  overwrite              : ', overwrite
    endif

    ! Add ice height to nodal point: 
    if (overwrite)then 
     nodalsl(rgnum_fs(i_gll, i_elmtfs)) = height
    else
     nodalsl(rgnum_fs(i_gll, i_elmtfs)) = nodalsl(rgnum_fs(i_gll, i_elmtfs)) + height
    endif 

    if(myrank.eq.0)then
        write(*,*)' ✓ Injected at GLL point'
    endif 

end subroutine add_sl_gll

! ################### END INITIAL SETUP FUNCTIONS  #####################






subroutine update_ocean_function(nodalice, nodalsl, errcode, errtag)
! Routine checks each GLL point on the surface to see if it is part of the ocean set
! see Crawford et al 2018, eqn 31-32.
! Set contains any nodes in which rho_w * SL > rho_i * I 
use global 
use free_surface
#if(USE_MPI)
use mpi_library
#else
use serial_library
#endif

implicit none 

! IO variables
character(len=250) :: errtag
integer :: ios, errcode
real(kind=kreal), allocatable :: nodalice(:), nodalsl(:)


! Local variables 
integer          :: i_elmt, iface, numf(maxngll2d), i_numf, i_node, nfgll
integer          :: i_gll  
real(kind=kreal) :: SL, I


! In this case we can calculate ocean function using nodalsl and nodalice: 
oceannodes =0 
do i_elmt=1, nelmt_fs  

    ! Face number (ie between 1 and 6) and get related properties
    call get_fs_details_noweights(i_elmt, iface, nfgll)

    ! Loop through GLL on the surface: 
    do i_gll = 1, nfgll
        SL =  nodalsl(rgnum_fs(i_gll, i_elmt))
        I     =  nodalice(rgnum_fs(i_gll, i_elmt))

        ! if rho_w SL > rho_ice I then part of ocean set
        if ( (rho_water * SL).GT.(I * rho_ice) ) then 
            oceanf(i_elmt, i_gll)            = 1.0_kreal
            nodalOF(rgnum_fs(i_gll, i_elmt)) = 1.0_kreal
            oceannodes = oceannodes + 1 
        else
            oceanf(i_elmt, i_gll)            = 0.0_kreal
            nodalOF(rgnum_fs(i_gll, i_elmt)) = 0.0_kreal
        endif 
    enddo 
enddo

if(myrank.eq.0)then 
    write(*,*)'  ✓ Updated ocean function'
    write(*,*)
endif 
end subroutine update_ocean_function



subroutine calculate_SL_A_per_proc(nodalsl, nodalu)
    ! Calculates the area covered by ocean (integral of ocean func
    ! over the solid surface)
    use global
    use element
    use mpi
    use set_precision_mpi
    use free_surface
    use integration
    use dimensionless
    use math_constants
#if (USE_MPI)
use mpi_library
#else
use serial_library
#endif

    implicit none 

    integer                        :: i_elmtfs, i_gll ! loops
    real(kind=kreal)               :: detjac2d, ocean_height! 2d jacobian
    integer                        :: iface           ! face ID for elmt 
    integer                        :: i_elmt, errcode          ! face ID for elmt 
    integer                        :: nfgll           ! ngll on 2D face
    real(kind=kreal), allocatable  :: gw(:), nodalsl(:), nodalu(:,:) ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:)
    real(kind=kreal)               :: coord(ndim,4), face_normal(3),& 
                                    dx_dxi(NDIM), dx_deta(NDIM)
    integer :: num4(4)

    ! Code
    allocate(gw(maxngll2d))
    allocate(dshape4(2,4,maxngll2d))

    if(myrank.eq.0)then
        write(*,*)'Calculating ocean area and volume'
    endif 

    ! Store old values
    SLarea_old   = SLarea   
    SLvolume_old = SLvolume

    SLarea   = ZERO 
    SLvolume = ZERO 


    do i_elmtfs = 1, nelmt_fs

        call get_fs_details(i_elmtfs, iface, nfgll, gw, dshape4)
        
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
            
            ! NOTE THAT THE SLAREA is actually the integral of the ocean function, 
            ! not the actual area (which requires projection of the normal into the local vertical)
            detjac2d=sqrt(dot_product(face_normal,face_normal))       
            SLarea   = SLarea + oceanf(i_elmtfs, i_gll)*gw(i_gll)*detjac2d
            ocean_height = nodalsl(rgnum_fs(i_gll, i_elmtfs)) - nodalu(3, rgnum_fs(i_gll, i_elmtfs))
            
            ! Project to the vertical (multiply by 0, 0, 1 for z as vertical): 
            face_normal(1) = zero
            face_normal(2) = zero
            detjac2d=sqrt(dot_product(face_normal,face_normal))       

            SLvolume     = SLvolume +  oceanf(i_elmtfs, i_gll)*gw(i_gll)*detjac2d*ocean_height
        enddo ! i_gll
    enddo   ! i_elmtfs

    ! calculate the change in mass
    SLmasschange = (SLvolume - SLvolume_old)*rho_water

    deallocate(gw)
    deallocate(dshape4)
    
end subroutine calculate_SL_A_per_proc





subroutine update_SL_area(nodalsl, nodalu)
    
use set_precision
use dimensionless
use math_constants
use set_precision_mpi
#if (USE_MPI)
use mpi_library
use mpi

#else
use serial_library
#endif
    real(kind=kreal), allocatable  :: nodalsl(:), nodalu(:,:) 
    integer :: errcode

    ! Use ocean function to calculate area of ocean for each processor     
    call calculate_SL_A_per_proc(nodalsl, nodalu)
    call sync_process()
    
  
    ! Sum up Area over all of the nodes: 
    call MPI_Allreduce(SLarea, totalSLA, 1, MPI_KREAL, MPI_SUM, MPI_COMM_WORLD, errcode) 
    call MPI_Allreduce(SLmasschange, SLsummasschange, 1, MPI_KREAL, MPI_SUM, MPI_COMM_WORLD, errcode) 
    SLarea = totalSLA
    SLmasschange = SLsummasschange

    ! Escape if no water. 
    if(SLarea.le.ZERO)then 
      write(*,*)'ERROR: Volume/Area of ocean = 0 -- NO WATER!!!' 
      write(*,*)'SL Area              : ', SLarea 
      write(*,*)'Dimensional SL Area  : ', SLarea * DIM_L * DIM_L 
      write(*,*)'The assumption is that there is at least some defined ocean basin.' 
      stop 
    endif 

    if (myrank.eq.0)then 
        write(*,'(a,g0.6)')'Integrated ocean function (A)    :   ', SLarea
        write(*,'(a,g0.6)')'Total SL mass change across nodes:   ', SLmasschange
        if(devel_nondim)then
            write(*,*)'    - Dimensionalised values: '
            write(*,'(a,g0.6)')'  --> Integrated ocean function (A): ', SLarea*DIM_L*DIM_L
            write(*,'(a,g0.6)')'  --> Total SL mass change         : ', SLmasschange*DIM_M
          endif 
        write(*,*) 
    endif 

end subroutine update_SL_area




end module