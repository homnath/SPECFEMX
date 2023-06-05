module ice 
    use global 
    use sea_level
    use set_precision
    implicit none 

    contains 


    

! ################### LOG AND OUTPUT FUNCTIONS  #######################

subroutine print_ice_read()
! Opens the ice log file and writes initial info 
use global
#if(USE_MPI)
use mpi_library
#else 
use serial_library
#endif  
implicit none 
    
if(myrank.eq.0)then
    ! Write summary of ICE inputs: 
    write(*,*)'-----------------------------------------------------'
    write(*,*)'Ice data read from:  ', trim(icefile)
    write(*,*)'-----------------------------------------------------'
    write(*,*)'IS_ICE            :  ', IS_ICE
    write(*,*)'Save original ice :  ', savedata%ice0
    write(*,*)'Save ice          :  ', savedata%ice
    write(*,*)'Save ice rate     :  ', savedata%icerate
    write(*,*)'Save ice load     :  ', savedata%iceload
    write(*,*)'-----------------------------------------------------'
    write(*,*)
endif 
end subroutine print_ice_read
        


subroutine write_ice_to_ensight(nodalice, i_step)
    use global 
    use postprocess
    use set_precision
    use free_surface
    use nondimensionpar
    implicit none 
    integer :: i_step
    real(kind=kreal) :: nodalice(:)

 
    if(savedata%fsplot)then
        call write_scalar_to_file_freesurf(nnode_fs,  DIM_L*nodalice, &
        ext='ice',istep=i_step) 
    endif

    if(savedata%fsplot_plane)then
        call write_scalar_to_file_freesurf(nnode_fs,  DIM_L*nodalice, &
        ext='ice', istep=i_step, plane=.true.) 
    endif

    if(myrank.eq.0.and.verbose_save_var)then
        write(*,'(a,i6)')'  ✓ Saved nodal ice for step ', i_step
        write(*,*)
    endif 
end subroutine write_ice_to_ensight


subroutine write_icerate_to_ensight(nodalicerate, i_step)
    use global 
    use postprocess
    use free_surface
    use nondimensionpar
    implicit none 
    real(kind=kreal) :: nodalicerate(:) ! Nodal rate of I values
    integer :: i_step


    ! On the free surface
    if(savedata%fsplot)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalicerate*DIM_L, &
        ext='icerate',istep=i_step) 
    endif
    
    if(savedata%fsplot_plane)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalicerate*DIM_L, &
        ext='nodalice', istep=i_step,plane=.true.) 
    endif


    if(myrank.eq.0.and.verbose_bool)then
        write(*,'(a,i6)')'  ✓ Saved ice rate for step ', i_step
        write(*,*)
    endif 

end subroutine write_icerate_to_ensight
! ################# END  LOG AND OUTPUT FUNCTIONS  ####################




! #################    INITIAL SETUP FUNCTIONS    #####################
subroutine prepare_ice(nodalice, nodalicerate)
    use global 
    use free_surface
    use set_precision
    use math_constants
    use nondimensionpar
    implicit none 

    integer :: istattemp, istat
    real(kind=kreal), allocatable :: nodalice(:), nodalicerate(:)

    if(myrank.eq.0.and.verbose_bool)then
        write(*,*)'Preparing ice variables...'
    endif
    istat = 0

   
    if(IS_ICE)then 
        allocate(nodalicerate(nnode_fs), nodalice(nnode_fs), stat=istattemp)
        nodalicerate = ZERO
        nodalice     = ZERO
        istat=istat+istattemp
    endif

 
    ! Check allocations 
    if(istat/=0)then
        write(*,*)'ERROR: cannot allocate memory in prepare_ice!'
        stop
    else 
        nodalice = ZERO
        ! Output confirmation to log. 
        if(myrank.eq.0.and.verbose_bool)then
            write(*, *)'  --> Created initial ICE (nodal)'
            write(*, *)'  --> Created ICE (nodal)'
        endif 
    endif



    if (savedata%iceload)then 
        ! Arrays for saving iceload 
        allocate(nodal_iceload_u(ndim, nnode_fs), & 
                 nodal_iceload_phi(nnode_fs),     &
                 nodal_iceload_sl(nnode_fs),      & 
                 stat=istattemp)
    endif 


       ! Check allocations 
    if(istat/=0)then
        write(*,*)'ERROR: cannot allocate memory for saving iceload!'
        stop
    else 
        if(myrank.eq.0.and.verbose_bool)then
            write(*, *)'  --> Created iceload_save arrays'
        endif 
    endif


    ! Apply non-dimensionalisaiton: 
    icerate = icerate*NONDIM_L


    if(myrank.eq.0)then
        write(*,*)
        write(*,*)'-  Using ice rate: ', icerate
        write(*,*)'  ✓ Prepared ice '
        write(*,*)
    endif 

    return 
end subroutine prepare_ice


subroutine set_original_ice_level(nodalice)
    ! Takes the user inputted ice data read from ice file and sets initial ice values
    ! Uses
    use set_precision
    use global 
    use integration
    use free_surface
    use math_constants
    use nondimensionpar
#if(USE_MPI)
use math_library_mpi
use mpi_library
use mpi
#else
use math_library_serial
use serial_library
#endif
    real(kind=kreal) :: nodalice(:)
    
    ! Local vars: 
    integer :: i_obj ! loop var
    integer :: iceobjtype
    real(kind=kreal) :: params(4), maxvalue, minvalue

    ! Code: 
    if(myrank.eq.0)then
    write(*,*)'-----------------------------------------------------'
    write(*, *)'        Setting original ice distribution           '
    write(*,*)'-----------------------------------------------------'
    endif 
    
    ! Loop through each ice object user inputted :
    do i_obj = 1, nice_obj
        iceobjtype = iceobjs(i_obj, 1)
        
        ! All of the parameters need to be non_dimensionalised in length 
        params = iceobjs(i_obj,2:5)*NONDIM_L


        if (iceobjtype.eq.1) then 
            ! Single point of ice: args: nodeid, height  
            call add_ice_gll(INT(params(1)), INT(params(2)), params(3), nodalice)
        elseif (iceobjtype.eq.2) then 
            ! Cylinder - args: x, y, rad, height
            call add_ice_cylinder(params, nodalice)
        elseif (iceobjtype.eq.3)then 
            ! Gaussian - args: x, y, amp, sigma
            call add_ice_gaussian(params, nodalice)
        elseif (iceobjtype.eq.4)then 
            ! Rectangle - args: x, y, len, height
            ! Assumes square base: xlen = y len
            call add_ice_sqcuboid(params, nodalice)

        else
            ! Invalid entry
            write(*,*)'ERROR: Unknown ice object type: ', iceobjtype
            stop
        endif 
    enddo 
    

    maxvalue = maxscal(maxval(nodalice))
    minvalue = minscal(minval(nodalice))
    ! update log file with results: 
    if(myrank.eq.0)then 
        write(*,*)
        write(*,*)'✓ Finished setting original ice level '
        write(*,'(a,i6)')'  -->  Number of ice objects added   : ', nice_obj
        write(*,'(a,g0.6)')'  -->  Min ice level                 : ', minvalue
        write(*,'(a,g0.6)')'  -->  Max ice level                 : ', maxvalue
        if(devel_nondim)then
            write(*,'(a,g0.6)')'  -->  Dimensionalised min ice level : ', minvalue*DIM_L
            write(*,'(a,g0.6)')'  -->  Dimensionalised ice level     : ', maxvalue*DIM_L
        endif 
        write(*,*)'-----------------------------------------------------'
        write(*,*)
    endif 

end subroutine set_original_ice_level


subroutine add_ice_gll(i_elmtfs, i_gll, height, nodalice)
    ! Adds ice in the required location to a single GLL point 
    ! Uses
    use set_precision
    use global 
    use integration
    use free_surface
    use math_constants
    use nondimensionpar 

    ! IO vars: 
    real(kind=kreal) :: height, nodalice(:) 
    integer :: i_elmtfs, i_gll

    ! Params should be the faceID, nodeID and the height
    ! Add to log file: 
    if(myrank.eq.0)then
        write(*,*)
        write(*,*)' --  Adding ice at point '
        write(*,'(a,i6)')'      ->  FS Elmt ID                : ', i_elmtfs
        write(*,'(a,i6)')'      ->  GLL Node (1 to maxngll2d) : ', i_gll
        write(*,'(a,g0.6)')'      ->  height                    : ', height
    endif 

    ! Add ice height to nodal point: 
    nodalice(rgnum_fs(i_gll, i_elmtfs)) = height

    if(myrank.eq.0.and.verbose_bool)then
        write(*,*)' ✓ Injected at GLL point'
        write(*,*)
    endif 
end subroutine add_ice_gll




subroutine add_ice_cylinder(params, nodalice)
    ! Adds a cylinder of ice in the required location
    ! Uses
    use set_precision
    use global 
    use nondimensionpar
    use integration
    use free_surface
    use math_constants

    ! IO vars: 
    real(kind=kreal) :: params(4), nodalice(:) ! x, y, rad, height

    ! Local vars: 
    integer :: i_obj,i_elmtfs,i_gll, iface,nfgll ,node_ctr
    integer :: iceobjtype
    integer :: ios, errcode

    real(kind=kreal):: x, y, r, h, dis, x_coord, y_coord, coord(3)


    ! Cylinder params 
    x = params(1)
    y = params(2)
    r = params(3)
    h = params(4)

    ! Add to log file: 
    if(myrank.eq.0)then
        write(*,*)' --  Creating cylinder '
        write(*,'(a,g0.6,1x,g0.6)')'      ->  centre (x,y): ', x, y
        write(*,'(a,g0.6)')'      ->  height      : ', h
        write(*,'(a,g0.6)')'      ->  radius      : ', r

        if(devel_nondim)then
            write(*,*)' --  Dimensionalised values '
            write(*,'(a,g0.6,1x,g0.6)')'      ->  centre (x,y): ', x*DIM_L, y*DIM_L
            write(*,'(a,g0.6)')'      ->  height      : ', h*DIM_L
            write(*,'(a,g0.6)')'      ->  radius      : ', r*DIM_L
        endif 
    endif 
    
    ! Searches for nodes on FS that are within the radius of cylinder
    ! Loop through each face on the free surface
    node_ctr = 0
    do i_elmtfs=1, nelmt_fs  

        ! For each GLL point get the coordinates
        call get_fs_details_noweights(i_elmtfs, iface, nfgll)
        do i_gll = 1, nfgll 
            ! Get X, Y, Z coordinates for the face 
            coord   = g_coord(:,  gnum_fs(i_gll,i_elmtfs))
            x_coord = coord(1)
            y_coord = coord(2)

            ! Cartesian distance to the point: 
            dis = ((x_coord - x)**2 + (y_coord - y)**2)**0.5

            if (dis.LE.r)then 
                ! Add ice height to nodal point: 
                nodalice(rgnum_fs(i_gll, i_elmtfs)) = h
                node_ctr = node_ctr + 1 
            endif 
        enddo 
    enddo 

    if(myrank.eq.0.and.verbose_bool)then
        write(*,*)' ✓ Injected cylinder at ', node_ctr, 'nodal points'
    endif
end subroutine add_ice_cylinder




subroutine add_ice_sqcuboid(params, nodalice)
    ! Adds a square-based cuboid of ice in the required location
    ! Uses
    use set_precision
    use global 
    use nondimensionpar
    use integration
    use free_surface
    use math_constants

    ! IO vars: 
    real(kind=kreal) :: params(4), nodalice(:) ! x, y, len, height

    ! Local vars: 
    integer :: i_obj,i_elmtfs,i_gll, iface,nfgll ,node_ctr
    integer :: iceobjtype
    integer :: ios, errcode

    real(kind=kreal):: x, y, r,r2, h, disx, disy, x_coord, y_coord, coord(3)


    ! Cylinder params 
    x = params(1)
    y = params(2)
    r = params(3)
    h = params(4)

    r2=r/two

    ! Add to log file: 
    if(myrank.eq.0)then
        write(*,*)' --  Creating square-based cuboid '
        write(*,'(a,g0.6,1x,g0.6)')'      ->  centre (x,y): ', x, y
        write(*,'(a,g0.6)')'      ->  height      : ', h
        write(*,'(a,g0.6)')'      ->  sq. length  : ', r

        if(devel_nondim)then
            write(*,*)' --  Dimensionalised values '
            write(*,'(a,g0.6,1x,g0.6)')'      ->  centre (x,y): ', x*DIM_L, y*DIM_L
            write(*,'(a,g0.6)')'      ->  height      : ', h*DIM_L
            write(*,'(a,g0.6)')'      ->  sq. length  : ', r*DIM_L
        endif 
    endif 
    
    ! Searches for nodes on FS that are within the domain of the square based
    ! Loop through each face on the free surface
    node_ctr = 0
    do i_elmtfs=1, nelmt_fs  

        ! For each GLL point get the coordinates
        call get_fs_details_noweights(i_elmtfs, iface, nfgll)
        do i_gll = 1, nfgll 
            ! Get X, Y, Z coordinates for the face 
            coord   = g_coord(:,  gnum_fs(i_gll,i_elmtfs))
            x_coord = coord(1)
            y_coord = coord(2)

            ! Cartesian distance to the point: 
            disx = ABS(x - x_coord) 
            disy = ABS(y - y_coord) 
            

            if (disx.LE.r2.and.disy.LE.r2)then 
                ! Add ice height to nodal point: 
                nodalice(rgnum_fs(i_gll, i_elmtfs)) = h
                node_ctr = node_ctr + 1 
            endif 
        enddo 
    enddo 

    if(myrank.eq.0)then
        write(*,*)' ✓ Injected square-based cuboid at ', node_ctr, 'nodal points'
    endif
end subroutine add_ice_sqcuboid



subroutine add_ice_gaussian(params, nodalice)
    ! Adds a gaussian of ice in the required location
    ! Uses 2D gaussian: (Amp/2 pi sigma^2) * exp(- (x^2 + y^2)/(2sigma^2) ) 
    ! input params are x, y, amp, sigma
    use set_precision
    use global 
    use nondimensionpar
    use integration
    use free_surface
    use math_constants

    ! IO vars: 
    real(kind=kreal) :: params(4), nodalice(:) ! x, y, rad, height

    ! Local vars: 
    integer :: i_obj,i_elmtfs,i_gll, iface,nfgll ,node_ctr
    integer :: iceobjtype
    integer :: ios, errcode

    real(kind=kreal):: x, y, amp, sigma, coeff, dis, x_coord, y_coord, coord(3)


    ! Gaussian params 
    x       = params(1)
    y       = params(2)
    amp     = params(3)*NONDIM_L*NONDIM_L ! To maintain mass/vol of gaussian need to div by L^3 (already divided by L)
    sigma   = params(4)

    ! Add to log file: 
    if(myrank.eq.0)then
        write(*,*)
        write(*,*)' --  Creating Gaussian '
        write(*,'(a,g0.6,1x,g0.6)')'      ->  centre (x,y): ', x, y
        write(*,'(a,g0.6)')'      ->  amplitude   : ', amp
        write(*,'(a,g0.6)')'      ->  sigma       : ', sigma

        if(devel_nondim)then
            write(*,*)' --  Dimensionalised values '
            write(*,'(a,g0.6,1x,g0.6)')'      ->  centre (x,y): ', x*DIM_L, y*DIM_L
            write(*,'(a,g0.6)')'      ->  amplitude   : ', amp*DIM_L*DIM_L*DIM_L
            write(*,'(a,g0.6)')'      ->  sigma       : ', sigma*DIM_L
        endif 
    endif
    
    ! Constant amplitude of gaussian:
    if(sigma.ne.ZERO)then
        coeff = (amp/(2*PI * (sigma**2))) 
    else
        write(*,*)'ERROR: TRYING TO ADD GAUSSIAN WITH SIGMA = 0'
        stop
    endif 


    ! Loop through each face on the free surface
    node_ctr = 0
    do i_elmtfs=1, nelmt_fs  

        ! For each GLL point get the coordinates
        call get_fs_details_noweights(i_elmtfs, iface, nfgll)
        do i_gll = 1, nfgll 
            ! Get X, Y, Z coordinates for the face 
            coord   = g_coord(:,  gnum_fs(i_gll,i_elmtfs))
            x_coord = coord(1)
            y_coord = coord(2)

            ! Add ice height to nodal point: 
            nodalice(rgnum_fs(i_gll, i_elmtfs)) =  coeff* EXP(- ( ( (x_coord-x)**2) + ((y_coord-y)**2))/(2 * sigma**2) ) 
        enddo 
    enddo 

    if(myrank.eq.0)then
        write(*,*)' ✓ Overwritten with gaussian'
        write(*,*)
    endif
    
end subroutine add_ice_gaussian




! ################### END INITIAL SETUP FUNCTIONS  #####################





subroutine set_ice_rate(nodalice, nodalicerate, i_step)
    ! Uses:
    ! Removes a constant amount from anywhere with ice
    ! if less than that amount of ice present, removes all of it 
    use set_precision
    use global 
    use free_surface
    use nondimensionpar
    use math_constants 
#if(USE_MPI)
use math_library_mpi
use mpi_library
use mpi
#else
use math_library_serial
use serial_library
#endif
    ! IO variables: 
    real(kind=kreal) :: nodalice(:), nodalicerate(:),z_coord, current_ice, residual, coord(3), icerateval
    ! Local variables: 
    integer :: iface, nfgll, i_elmtfs, i_gll, gid, i_step
    ! Code


    ! Get icerateval for this timestep: 
    icerateval = icerate(i_step)

    if(myrank.eq.0)then
        write(*,*)'-----------------------------------------------------'
        write(*,*)'             Setting ice rate/change                '
        write(*,*)'-----------------------------------------------------'
        write(*,*)'CAUTION: ONLY IMPLEMENTING FIXED ICE RATE ACROSS REGIONS WITH ICE'
        write(*,*)
        write(*,'(a,g0.6)')'* Using fixed ice rate value        : ', icerateval
        if(devel_nondim)then 
            write(*,*)'  --> dimensionalised ice rate value: ', icerateval*DIM_L
        endif 
        write(*,*)'* If larger than height of ice, just removes all ice.'
        write(*,*)
    endif 

    ! Loop elements
    do i_elmtfs = 1, nelmt_fs  
        ! Get number of GLL 
        call get_fs_details_noweights(i_elmtfs, iface, nfgll)
        ! Loop GLL 
        do i_gll = 1, nfgll 
            ! Get X, Y, Z coordinates for the face 
            coord   = g_coord(:,  gnum_fs(i_gll,i_elmtfs))
            z_coord = coord(3)

            gid = rgnum_fs(i_gll, i_elmtfs)
            current_ice = nodalice(gid)

            ! Only for locations with ice present: 
            if (current_ice.gt.zero) then 
                residual  = (current_ice+icerateval)
                

                if(residual.le.zero)then 
                    !Just remove all of what is present currently
                    nodalicerate(gid) = -current_ice
                else
                    ! Trim a bit off the top!
                    nodalicerate(gid) = icerateval
                endif
            endif ! If some ice exists
        enddo !igll
    enddo !inode_fs

    !if(myrank.eq.0)then
    !    write(*,*)' ✓ Finished setting ice rate/change'
    !    write(*,*)'-----------------------------------------------------'
    !    write(*,*)
    !endif 
end subroutine set_ice_rate




subroutine set_ice_rate_slice(nodalice, nodalicerate, i_step)
    ! Uses:
    ! Removes a slice of ice from the top
    use set_precision
    use global 
    use free_surface
    use nondimensionpar
    use math_constants 
#if(USE_MPI)
use math_library_mpi
use mpi_library
use mpi
#else
use math_library_serial
use serial_library
#endif

    ! IO variables: 
    real(kind=kreal) :: nodalice(:), nodalicerate(:),z_coord, current_ice, residual, threshold, coord(3),maxice, icerateval
    integer :: i_step 
    ! Local variables: 
    integer :: iface, nfgll, i_elmtfs, i_gll, gid
    ! Code

   ! Get maxvalue of current ice: 
    maxice=maxscal(maxval(nodalice))

    ! Get icerate for this timestep
    icerateval = icerate(i_step)

    ! Threshold for slice
    threshold = maxice + icerateval

    if(myrank.eq.0.and.verbose_bool)then
        write(*,*)'-----------------------------------------------------'
        write(*,*)'             Setting ice rate/change                '
        write(*,'(a,g0.6)')'* Using fixed ice rate value   : ', icerateval
        write(*,'(a,g0.6)')'* Threshold                    : ', threshold
        write(*,'(a,g0.6)')'* maxice                       : ', maxice
        if(devel_nondim)then 
            write(*,*)'  Dimensional values...'
            write(*,'(a,g0.6)')'   * Using fixed ice rate value: ', icerateval*DIM_L
            write(*,'(a,g0.6)')'   * Threshold                 : ', threshold*DIM_L
            write(*,'(a,g0.6)')'   * maxice                    : ', maxice*DIM_L
        endif 
    endif 
 

    if (maxice.gt.zero)then 
        ! Loop elements
        do i_elmtfs = 1, nelmt_fs  
            ! Get number of GLL 
            call get_fs_details_noweights(i_elmtfs, iface, nfgll)
            ! Loop GLL 
            do i_gll = 1, nfgll 
                ! Get X, Y, Z coordinates for the face 
                coord   = g_coord(:,  gnum_fs(i_gll,i_elmtfs))
                z_coord = coord(3)

                gid = rgnum_fs(i_gll, i_elmtfs)
                current_ice = nodalice(gid)

                
                if(threshold.le.zero) then 
                    ! Remove what is left
                    nodalicerate(gid) = -current_ice
                elseif (current_ice.gt.threshold) then 
                    ! If within top threshold section, reduce down to threshold value 
                    nodalicerate(gid) = threshold-current_ice
                else 
                    nodalicerate(gid) = zero 
                endif 

            enddo !igll
        enddo !inode_fs
    endif 



    if(myrank.eq.0.and.verbose_bool)then
        write(*,*)' ✓ Finished setting ice rate/change'
        write(*,*)'-----------------------------------------------------'
        write(*,*)
    endif
end subroutine set_ice_rate_slice








subroutine calculate_ice_change_volume(nodalicerate)
    ! Simple integration of icerate change over the FS 

    use global
    use element
    use free_surface
    use integration
    use math_constants

    implicit none 

    real(kind=kreal)               :: nodalicerate(:)

    integer                        :: i_elmtfs, i_gll ! loops
    real(kind=kreal)               :: detjac2d, ocean_height! 2d jacobian
    integer                        :: iface           ! face ID for elmt 
    integer                        :: i_elmt            ! face ID for elmt 
    integer                        :: nfgll             ! ngll on 2D face
    real(kind=kreal), allocatable  :: gw(:), nodalsl(:) ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:)
    real(kind=kreal)               :: coord(ndim,4), face_normal(3),& 
                                      dx_dxi(NDIM), dx_deta(NDIM)
    integer :: num4(4)

    ! Code
    allocate(gw(maxngll2d))
    allocate(dshape4(2,4,maxngll2d))

    if(myrank.eq.0.and.verbose_bool)then 
        write(*,*)'* Calculating change in ice volume'
        write(*,*)' --> WARNING: USING PROJECTION OF AREA TO THE Z direction'
    endif 

    icechangevol = ZERO 

    do i_elmtfs = 1, nelmt_fs

        call get_fs_details(i_elmtfs, iface, nfgll, gw, dshape4)
        num4   = gnum4_fs(:, i_elmtfs)
        coord  = g_coord(:,num4)

        do i_gll = 1, nfgll 
            ! Calculate the magnitude of the 2D jacobian 
            dx_dxi  = matmul(coord,dshape4(1,:,i_gll))
            dx_deta = matmul(coord,dshape4(2,:,i_gll))
            face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
            face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
            face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

            ! Project to the vertical (multiply by 0, 0, 1 for z as vertical): 
            face_normal(1) = zero; 
            face_normal(2) = zero;
            detjac2d=sqrt(dot_product(face_normal,face_normal))       
            
            icechangevol = icechangevol + (gw(i_gll) * detjac2d * nodalicerate(rgnum_fs(i_gll, i_elmtfs)))
        enddo ! i_gll
    enddo   ! i_elmtfs
    

    deallocate(gw)
    deallocate(dshape4)
end subroutine calculate_ice_change_volume


subroutine calc_ice_load(nodalicerate,i_step)
! Uses 
use global 
use postprocess
use element
use set_precision 
use free_surface
use math_constants
use nondimensionpar
#if(USE_MPI)
use math_library_mpi
use mpi_library
use mpi
#else
use math_library_serial
use serial_library
#endif

    ! IO vars
    real(kind=kreal)               :: nodalicerate(:)
    integer :: i_step
    ! Local vars
    integer                        :: i_elmtfs,nodeid,dof, gid_elmt, num4(4), num_FS(maxngll2d), gid_abg, gid_xyg, iface, nfgll, abg,xyg, j,k, i_gll, i
    real(kind=kreal)               :: coord(ndim,4), face_normal(3), dx_dxi(NDIM), dx_deta(NDIM), pi_2d_abg, pi_2d_xyg
    real(kind=kreal)               :: utf_dot_bkgrav, utf_dot_bkgrav_xyg, eps_area, xyg_sum
    real(kind=kreal), allocatable  :: gw(:)           ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:),miniceload,maxiceload
    real(kind=kreal)               :: epsilon, sumepsilon, pi_2d, ctr, val
    integer                        :: fgdof(nndof, maxngll2d), errcode ! face global degrees of freedom

    ! for projection to face normal: 
    real(kind=kreal) :: vertical(3), unit_normal(3), face_normal_len, cos_theta

    ! Debug: 
    integer :: iproc


    ! Code: 
    allocate(gw(maxngll2d))
    allocate(dshape4(2,4,maxngll2d))

    if(myrank.eq.0.and.verbose_bool)then
        write(*,*)
        write(*,*)'-----------------------------------------------------'
        write(*,*)'             Calculating iceload (RHS)               '
        write(*,*)
    endif 

    ! Initialise 
    iceload = zero
    if (savedata%iceload)then 
        ! Zero for saving iceload 
        nodal_iceload_u   = zero
        nodal_iceload_phi = zero
        nodal_iceload_sl  = zero
    endif 


    ! Epsilon/Area
    call calc_iceload_epsilon(epsilon, gw, dshape4, num4, coord, nodalicerate)
    call sync_process()
    ! Sum up Epsilon over all of the nodes and copy to local epsilon
    sumepsilon = sumscal(epsilon)
    epsilon    = sumepsilon

    if(myrank.eq.0.and.verbose_bool)then
        write(*,'(a,g0.8)')'  * ε value                 : ', epsilon
        if(devel_nondim)then
            write(*,'(a,g0.8)')'  * ε value dimensionalised : ', epsilon*DIM_L*DIM_L*DIM_L
        endif 
        write(*,*)
    endif 


    eps_area = epsilon/SLarea



    do i_elmtfs=1,nelmt_fs

        ! GET SOME DETAILS ABOUT THE ELEMENT 
        call get_fs_details(i_elmtfs, iface, nfgll, gw, dshape4) ! face number, number of GLL on face, gauss weights, derivative of shape funcs
        num4     = gnum4_fs(:, i_elmtfs)           ! node IDs for corners 
        coord    = g_coord(:,num4)                 ! node coordinates
        gid_elmt = id_elem_fs(i_elmtfs)            ! global element ID 
        num_FS   = gnum_fs(:, i_elmtfs)            ! Global node IDs of the GLL pts on FS 

        ! DOF IDs for the nodes in matrix (property, gll pt) where property goes from 1 - 5 (ux,uy,uz,phi,theta)
        fgdof(:, 1:nfgll) = reshape(gdof(:, g_num(hexface(iface)%node, gid_elmt)),(/nndof, maxngll2d/))


        do i_gll = 1, nfgll 

            nodeid = rgnum_fs(i_gll, i_elmtfs)

            ! Calculate the magnitude of the 2D jacobian 
            dx_dxi  = matmul(coord,dshape4(1,:,i_gll))
            dx_deta = matmul(coord,dshape4(2,:,i_gll))
            ! Calc normal and therefore jac dec (2D) on the fly
            face_normal(1) = dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
            face_normal(2) = dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
            face_normal(3) = dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
            pi_2d          = gw(i_gll) * sqrt(dot_product(face_normal,face_normal)) ! Weights*jac

            !face_normal_len =  sqrt(dot_product(face_normal,face_normal))

            ! We need the dot product of the vertical with the normal to the surface 
            ! The sign is not relevant becaause the water is always pushing down from the top surface
            !vertical    = zero
            !vertical(3) = one

            ! Normalise the length of the face normal to get unit normal to the free surface
            !unit_normal = face_normal/face_normal_len
            !cos_theta   = ABS(dot_product(unit_normal,vertical) ) 


            ! Calculate the coefficient for phi and u that is shared
            ! Note that u values also need to be multiplied by background gravity
            val =( (ONE - oceanf(i_elmtfs, i_gll)) * nodalicerate(nodeid)) - eps_area*oceanf(i_elmtfs, i_gll) 
            val = val * pi_2d !* cos_theta

            ! Displacement for direction j: + ( (1-OF)*I_dot  - epsilon/A * OF  )* pi * \nabla\Phi_j
            do j = 1, NDIM    
                dof = fgdof(j, i_gll) 
                if (dof.gt.0)then 

                    iceload(dof) = iceload(dof) + (val *  (-grav0_nodal(j, num_FS(i_gll))) )
                    
                    if (savedata%iceload)then 
                        nodal_iceload_u(j,nodeid) = nodal_iceload_u(j, nodeid) + (val *  (-grav0_nodal(j, num_FS(i_gll)) ) )
                    endif 
                endif  
            enddo

            
            ! Phi:  + ( (1-OF)*I_dot  - epsilon/A * OF  )* pi
            dof = fgdof(4, i_gll) 
            if (dof.gt.0)then 

                iceload(dof) = iceload(dof) + val  
                if (savedata%iceload)then 
                    nodal_iceload_phi(nodeid) = nodal_iceload_phi(nodeid) + val
                endif 
            endif 


            ! Theta: + epsilon/A * g * pi2d
            dof = fgdof(5, i_gll)
            if (dof.gt.0)then 


                iceload(dof) = iceload(dof) - (eps_area * pi_2d * ABS(g0_nodal(num_FS(i_gll))) )

                if (savedata%iceload)then 
                    nodal_iceload_sl(nodeid) = nodal_iceload_sl(nodeid) - (eps_area * pi_2d * ABS(g0_nodal(num_FS(i_gll))))
                endif 
            endif



        enddo  ! i_gll 
    enddo  ! i_elmtfs





    ! Multiply whole of the vector by rho_i 
    iceload = -(rho_ice * iceload)

    call sync_process()

    miniceload = minscal(minval(iceload))
    maxiceload = maxscal(maxval(iceload))


    if(myrank.eq.0.and.verbose_bool)then
        write(*,*)
        write(*,*)'✓ Calculated ice load:'
        write(*,'(a, g0.6)')'  -->  Min value of iceload     : ', miniceload
        write(*,'(a, g0.6)')'  -->  Max value of iceload     : ', maxiceload
        if(devel_nondim)then
            write(*,'(a)')' --  Dimensionalised values '
            write(*,'(a, g0.6)')'  -->  Min value of iceload     : ', miniceload * DIM_ICELOAD
            write(*,'(a, g0.6)')'  -->  Max value of iceload     : ', maxiceload * DIM_ICELOAD
        endif 
        write(*,*)
    endif 

    ! Save iceload to file: 
    if (savedata%iceload)then 
        nodal_iceload_u   =   -nodal_iceload_u   * rho_ice * DIM_ICELOAD
        nodal_iceload_phi =   -nodal_iceload_phi * rho_ice * DIM_ICELOAD
        nodal_iceload_sl  =   -nodal_iceload_sl  * rho_ice * DIM_ICELOAD

        call write_iceload_to_ensight(i_step=i_step-1)

        if(i_step.eq.nstep)then 
            ! Need to save at last timestep for ensight files to be 
            ! created properly/viewed in paraview, but these values are 
            ! irrelevant...would act on the (nstep+1)th timestep
            nodal_iceload_u   = zero 
            nodal_iceload_phi = zero
            nodal_iceload_sl  = zero
            call write_iceload_to_ensight(i_step=i_step)
        endif
    endif

    deallocate(gw)
    deallocate(dshape4)

end subroutine calc_ice_load    






subroutine write_iceload_to_ensight(i_step)
    use global
    use postprocess
    use free_surface
    use set_precision
    integer :: i_step

    ! Sea level iceload
    if(savedata%fsplot)then
      call write_scalar_to_file_freesurf(nnode_fs, nodal_iceload_sl,&
      ext='iceload_sl',istep=i_step) 
    endif
    if(savedata%fsplot_plane)then
      call write_scalar_to_file_freesurf(nnode_fs, nodal_iceload_sl,&
      ext='iceload_sl',istep=i_step,plane=.true.) 
    endif

    ! Save phi iceload
    if(savedata%fsplot)then
        call write_scalar_to_file_freesurf(nnode_fs, nodal_iceload_phi,&
        ext='iceload_phi',istep=i_step) 
      endif
    if(savedata%fsplot_plane)then
        call write_scalar_to_file_freesurf(nnode_fs, nodal_iceload_phi,&
        ext='iceload_phi',istep=i_step,plane=.true.) 
    endif

    ! Displacement iceload 
    if(savedata%fsplot)then
      call write_vector_to_file_freesurf(nnode_fs,nodal_iceload_u,&
      ext='iceload_u',istep=i_step)
    endif
    if(savedata%fsplot_plane)then
      call write_vector_to_file_freesurf(nnode_fs,nodal_iceload_u,&
      ext='iceload_u',istep=i_step,plane=.true.)
    endif

    if(myrank.eq.0.and.verbose_bool)then
        write(*,'(a,i6)')'  ✓ Saved ice load for step ', i_step
        write(*,*)
    endif 

end subroutine write_iceload_to_ensight


subroutine calc_iceload_epsilon(epsilon, gw, dshape4, num4, coord, nodalicerate)
    use free_surface
    use global
    use set_precision
    use math_constants
    ! IO: 
    integer                        :: num4(4)
    real(kind=kreal)               :: gw(:), nodalicerate(:)        
    real(kind=kreal)               :: dshape4(:,:,:)
    real(kind=kreal)               :: coord(ndim,4),  epsilon

    ! local vars: 
    integer :: i_elmtfs, iface, nfgll, i_gll
    real(kind=kreal) :: dx_dxi(NDIM), dx_deta(NDIM), pi_2d, val, face_normal(3) 

    ! First calculate the average ice rate load change (epsilon): 
    epsilon = zero
    do i_elmtfs=1,nelmt_fs
        ! Get details for the element: 
        call get_fs_details(i_elmtfs, iface, nfgll, gw, dshape4)
        num4   = gnum4_fs(:, i_elmtfs)
        coord  = g_coord(:,num4)

        ! Loop over GLL nodes of the face
        do i_gll = 1, nfgll 
            ! Calculate the magnitude of the 2D jacobian 
            dx_dxi  = matmul(coord,dshape4(1,:,i_gll))
            dx_deta = matmul(coord,dshape4(2,:,i_gll))

            ! Calc normal and therefore jac dec (2D) on the fly
            face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
            face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
            face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
            pi_2d         = gw(i_gll) * sqrt(dot_product(face_normal,face_normal)) ! Weights*jac

            val = (ONE - oceanf(i_elmtfs, i_gll)) * nodalicerate(rgnum_fs(i_gll, i_elmtfs)) * pi_2d
            epsilon = epsilon + val 
        enddo 
    enddo 

end subroutine calc_iceload_epsilon






subroutine summarise_ice_vol_change()
    use global 
    use nondimensionpar
#if(USE_MPI)
use mpi_library
use math_library_mpi
use ghost_library_mpi
#else
use serial_library
use math_library_serial
#endif

    ! Update the total mass change to date 
    icechangevoltmp = sumscal(icechangevol) 
    icechangevol = icechangevoltmp
    icemasschange_per_ts = icechangevol*rho_ice
    total_ice_mass_change = total_ice_mass_change + icemasschange_per_ts
    if(myrank.eq.0.and.verbose_bool)then
      write(*,'(a,g0.6)')'  --> Change in ice mass to occur: ', icemasschange_per_ts
      write(*,'(a,g0.6)')'  --> Total ice change so far    : ', total_ice_mass_change
      if(devel_nondim)then
        write(*,*)'    - Dimensionalised values: '
        write(*,'(a,g0.6)')'  --> Change in ice mass to occur: ', icemasschange_per_ts * DIM_M
        write(*,'(a,g0.6)')'  --> Total ice change so far    : ', total_ice_mass_change* DIM_M
      endif 
      write(*,*)
    endif   
end subroutine



end module
