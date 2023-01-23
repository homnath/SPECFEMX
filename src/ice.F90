module ice 
    use global 
    use sea_level
    use set_precision
    implicit none 

    contains 


    

! ################### LOG AND OUTPUT FUNCTIONS  #######################

subroutine start_ICE_log(errcode, errtag)
! Opens the ice log file and writes initial info 
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
        ICE_log_file = trim(file_head)//'ICE.log'
        open(unit=ICElogunit,file=trim(ICE_log_file),status='replace',action='write',iostat=ios)
        if(ios.ne.0)then
            write(errtag,'(a)')'ERROR: cannot open log file: '//trim(ICE_log_file)
            call control_error(errcode,errtag,stdout,myrank)
        endif
    endif 
    
    ! Write confirmation of file: 
    write(*, *)' Created ICE log file'
    write(ICElogunit, '(A,/,A)')' ****** CREATED ICE LOG FILE ****** ', ' '
    
    ! Write summary of ICE inputs: 
    write(ICElogunit,*)
    write(ICElogunit,*)'-----------------------------------------------------'
    write(ICElogunit,*)'IS_ICE            :  ', IS_ICE
    write(ICElogunit,*)'Save original ice :  ', savedata%ice0
    write(ICElogunit,*)'Save ice          :  ', savedata%ice
    write(ICElogunit,*)'Save ice_rate     :  ', savedata%icerate
    write(ICElogunit,*)'-----------------------------------------------------'
    write(ICElogunit,*)
    
    
    write(ICElogunit, *)'Ice data read from:  ', trim(icefile)
    if(IS_CART_SIM)then
        call summarise_ICE_input_cart()
    elseif(IS_GLOB_SIM)then 
        write(*,*) 'GLOBAL SIMULATIONS NOT IMPLEMETED YET'
        stop
    else
        write(*,*) 'SIMULATION MUST BE GLOBAL OR CARTESIAN'
        stop
    endif 
end subroutine start_ICE_log
        
subroutine summarise_ICE_input_cart()
    use global 
    implicit none 

    write(ICElogunit,*)
    write(ICElogunit,*)'Model setup          : Cartesian'
end subroutine summarise_ICE_input_cart



subroutine write_ICE0_to_ensight(nodalice)
    use global 
    use postprocess
    use set_precision
    use free_surface
    implicit none 

    real(kind=kreal),allocatable :: nodalice(:)

    write(ICElogunit,*)'Saving the original ICE values'
    

    ! On the free surface
    if(savedata%fsplot)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalice, &
        ext='ice0',istep=0) 
    endif
    
    if(savedata%fsplot_plane)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalice, &
        ext='ice0', istep=0,plane=.true.) 
    endif

    write(ICElogunit,*)'  ✓ Saved original ice level '
    write(ICElogunit,*)
end subroutine write_ICE0_to_ensight


subroutine write_icerate_to_ensight(nodalicerate)
    use global 
    use postprocess
    use free_surface
    implicit none 
    real(kind=kreal) :: nodalicerate(:) ! Nodal rate of I values


    write(ICElogunit,*)'Saving the icerate values'
    

    ! On the free surface
    if(savedata%fsplot)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalicerate, &
        ext='icerate',istep=0) 
    endif
    
    if(savedata%fsplot_plane)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalicerate, &
        ext='nodalice', istep=0,plane=.true.) 
    endif

    write(ICElogunit,*)'  ✓ Saved ice rate level '
    write(ICElogunit,*)
end subroutine write_icerate_to_ensight


! ################# END  LOG AND OUTPUT FUNCTIONS  ####################




! #################    INITIAL SETUP FUNCTIONS    #####################
subroutine prepare_ice(nodalice, nodalicerate)
    use global 
    use free_surface
    use set_precision
    use math_constants

    implicit none 

    integer :: istattemp, istat
    real(kind=kreal), allocatable :: nodalice(:), nodalicerate(:)

    write(ICElogunit, *)'Preparing ice variables...'
    
    istat = 0

   
    if(IS_ICE)then 
        allocate(nodalicerate(nnode_fs), nodalice(nnode_fs), stat=istattemp)
        nodalicerate = ZERO
        nodalice     = ZERO
        istat=istat+istattemp
    endif

 
    ! Check allocations 
    if(istat/=0)then
        write(logunit,*)'ERROR: cannot allocate memory in prepare_ice!'
        flush(logunit)
        stop
    else 
        nodalice = ZERO
        ! Output confirmation to log. 
        write(ICElogunit, *)'  --> Created initial ICE (nodal)'
        write(ICElogunit, *)'  --> Created ICE (nodal)'
    endif

    write(ICElogunit,*)'  ✓ Prepared ice '
    write(ICElogunit,*)

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

    real(kind=kreal) :: nodalice(:)
    
    ! Local vars: 
    integer :: i_obj ! loop var
    integer :: iceobjtype
    real(kind=kreal) :: params(4)

    ! Code: 

    ! Loop through each ice object user inputted :
    do i_obj = 1, nice_obj
        
        iceobjtype = iceobjs(i_obj, 1)

        if (iceobjtype.eq.1) then 
            ! Single point of ice 
            write(*,*)'ERROR: Single ice point not implemented yet!'
            stop
        elseif (iceobjtype.eq.2) then 
            ! Cylinder - args: x, y, rad, height
            params = iceobjs(i_obj,2:5)
            call add_ice_cylinder(params, nodalice)
        else
            ! Invalid entry
            write(*,*)'ERROR: Unknown ice object type: ', iceobjtype
            stop
        endif 
    enddo 


    ! update log file with results: 
    write(ICElogunit,*)''
    write(ICElogunit,*)'* Finished setting original ice level '
    write(ICElogunit,*)'  -->  Number of ice objects added : ', nice_obj
    write(ICElogunit,*)'  -->  Min ice level               : ', minval(nodalice)
    write(ICElogunit,*)'  -->  Max ice level               : ', maxval(nodalice)
end subroutine set_original_ice_level


subroutine add_ice_cylinder(params, nodalice)
    ! Adds a cylinder of ice in the required location
    ! Uses
    use set_precision
    use global 
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
    write(ICElogunit,*)
    write(ICElogunit,*)' --  Creating cylinder '
    write(ICElogunit,*)'      ->  centre (x,y): ', x, y
    write(ICElogunit,*)'      ->  height      : ', h
    write(ICElogunit,*)'      ->  radius      : ', r
    write(ICElogunit,*)'NEED TO CHECK DIMENSIONS (NONDIM) of cylinder coords.'

    
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

    write(ICElogunit,*)' ✓ Injected cylinder at ', node_ctr, 'nodal points'
    flush(ICElogunit)
end subroutine add_ice_cylinder

! ################### END INITIAL SETUP FUNCTIONS  #####################





subroutine set_ice_rate(nodalice, nodalicerate)
    ! Uses:
    use set_precision
    use global 
    use free_surface
    use math_constants 
    ! IO variables: 
    real(kind=kreal) :: nodalice(:), nodalicerate(:)  
    ! Local variables: 
    integer :: iface, nfgll, i_elmtfs, i_gll, gid
    ! Code

    write(ICElogunit, *)'CAUTION: ONLY IMPLEMENTING FIXED ICE RATE ACROSS REGIONS WITH ICE'
    write(ICElogunit, *)'Using fixed ice rate value:', icerateval

    ! Loop elements
    do i_elmtfs = 1, nelmt_fs  
        ! Get number of GLL 
        call get_fs_details_noweights(i_elmtfs, iface, nfgll)
        ! Loop GLL 
        do i_gll = 1, nfgll 
            gid = rgnum_fs(i_gll, i_elmtfs)
            ! Apply constant value to anywhere with ice: 
            if (nodalice(gid).gt.zero) then 
                nodalicerate(gid) = icerateval
            endif 
        enddo 
    enddo 


end subroutine set_ice_rate



subroutine calc_ice_load(iceload, nodalice)
! Uses 
use global 
use set_precision 
use free_surface
use math_constants
#if(USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif

    ! IO vars
    real(kind=kreal) iceload(:), nodalice(:)

    ! Local vars
    integer          ::  i_elmt, num4(4), gid_abg, gid_xyg, iface, nfgll, abg,xyg, j,k
    real(kind=kreal) :: coord(ndim,4), face_normal(3), dx_dxi(NDIM), dx_deta(NDIM), pi_2d_abg, pi_2d_xyg
    real(kind=kreal) :: utf_dot_bkgrav, utf_dot_bkgrav_xyg, area_inv, xyg_sum
    real(kind=kreal), allocatable  :: gw(:)           ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:)

    ! Code: 
    allocate(gw(maxngll2d))
    allocate(dshape4(2,4,maxngll2d))


    ! Inverse area
    area_inv = ONE/SLarea

    ! Loop for each FS element
    do i_elmt = 1, nelmt_fs

        ! Get details for the free surface element (e.g. which face etc)
        call get_fs_details(i_elmt, iface, nfgll, gw, dshape4)
        num4   = gnum4_fs(:, i_elmt)
        coord  = g_coord(:,num4)


        ! CALCULATE XYG SUM THAT IS USED FOR EACH ABG, but doesnt depend on ABG: 
        ! I think this can sit outside the ABG loop
        xyg_sum = ZERO  
        do xyg = 1, nfgll 
            gid_xyg = gnum_fs(xyg, i_elmt)

            ! Get Jacobian_2d x weights for XYG 
            dx_dxi  = matmul(coord,dshape4(1,:,xyg))
            dx_deta = matmul(coord,dshape4(2,:,xyg))
            face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
            face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
            face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
            pi_2d_xyg     =gw(xyg) * sqrt(dot_product(face_normal,face_normal)) ! Weights*jacw

            ! Sum of u_tf dot grad Phi at XYG
            utf_dot_bkgrav_xyg = ZERO 
            do k = 1, NDIM 
                utf_dot_bkgrav_xyg = utf_dot_bkgrav_xyg + (u_tf(k) * grav0_nodal(k, gid_xyg))
            enddo 

            ! Sum up! 
            xyg_sum = xyg_sum + pi_2d_xyg * (g0_nodal(gid_xyg) * theta_tf  + & 
                                            oceanf(i_elmt, xyg) * (phi_tf + utf_dot_bkgrav_xyg))
        enddo 



        ! Loop through each GLL on the face 
        do abg = 1, nfgll 
            gid_abg = gnum_fs(abg, i_elmt)


            ! Get Jacobian_2d x weights for ABG 
            dx_dxi  = matmul(coord,dshape4(1,:,abg))
            dx_deta = matmul(coord,dshape4(2,:,abg))
            face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
            face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
            face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
            pi_2d_abg     = gw(abg) * sqrt(dot_product(face_normal,face_normal)) ! Weights*jacw

            ! Sum over u_tf \cdot background grav @ abg nodes: 
            utf_dot_bkgrav = ZERO 
            do j = 1, NDIM 
                utf_dot_bkgrav = utf_dot_bkgrav + (u_tf(j) * grav0_nodal(j, gid_abg))
            enddo 

            
            iceload(gid_abg) = iceload(gid_abg) + ( (1 -  oceanf(i_elmt, abg)) * nodalice(rgnum_fs(abg, i_elmt)) * pi_2d_abg * & 
                               (phi_tf + utf_dot_bkgrav - area_inv * xyg_sum) ) 
        enddo ! abg



        ! Print some stats about iceload: 
        write(ICElogunit,*)' --------------------------------------------- '
        write(ICElogunit,*)' Calculated ice load: '
        write(ICElogunit,*)' Min value          : ', minscal(minval(iceload))
        write(ICElogunit,*)' Max value          : ', maxscal(maxval(iceload))


    enddo !nelmt_fs 



    deallocate(gw)
    deallocate(dshape4)

end subroutine calc_ice_load    

end module