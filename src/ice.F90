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
    write(ICElogunit,*)'Save ice rate     :  ', savedata%icerate
    write(ICElogunit,*)'Save ice load     :  ', savedata%iceload
    write(ICElogunit,*)'-----------------------------------------------------'
    write(ICElogunit,*)
    
    
    write(ICElogunit, *)'Ice data read from:  ', trim(icefile)
    if(IS_CART_SIM)then
        call summarise_ICE_input_cart()
    elseif(IS_GLOB_SIM)then 
        write(*,*) 'ERROR: GLOBAL SIMULATIONS NOT IMPLEMETED YET'
        stop
    else
        write(*,*) 'ERROR: SIMULATION MUST BE GLOBAL OR CARTESIAN'
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

    write(ICElogunit, *)
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
        params = iceobjs(i_obj,2:5)

        if (iceobjtype.eq.1) then 
            ! Single point of ice: args: nodeid, height  
            call add_ice_gll(INT(params(1)), INT(params(2)), params(3), nodalice)
        elseif (iceobjtype.eq.2) then 
            ! Cylinder - args: x, y, rad, height
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


subroutine add_ice_gll(i_elmtfs, i_gll, height, nodalice)
    ! Adds ice in the required location to a single GLL point 
    ! Uses
    use set_precision
    use global 
    use integration
    use free_surface
    use math_constants

    ! IO vars: 
    real(kind=kreal) :: height, nodalice(:) 
    integer :: i_elmtfs, i_gll

    ! Params should be the faceID, nodeID and the height
    ! Add to log file: 
    write(ICElogunit,*)
    write(ICElogunit,*)' --  Adding ice at point '
    write(ICElogunit,*)'      ->  FS Elmt ID             : ', i_elmtfs
    write(ICElogunit,*)'      ->  GLL Node (1-maxngll2d) : ', i_gll
    write(ICElogunit,*)'      ->  height                 : ', height


    ! Add ice height to nodal point: 
    nodalice(rgnum_fs(i_gll, i_elmtfs)) = height

    write(ICElogunit,*)' ✓ Injected at GLL point'
    flush(ICElogunit)
end subroutine add_ice_gll




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

    write(ICElogunit, *)
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



subroutine calc_ice_load(iceload, nodalicerate,nodalu, i_step)
! Uses 
use global 
use postprocess
use element
use set_precision 
use free_surface
use math_constants

#if(USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif

    ! IO vars
    real(kind=kreal),intent(inout) :: iceload(:)
    real(kind=kreal)               :: nodalicerate(:), nodalu(:,:)
    integer :: i_step
    ! Local vars
    integer                        :: i_elmtfs,nodeid,dof, gid_elmt, num4(4), num_FS(maxngll2d), gid_abg, gid_xyg, iface, nfgll, abg,xyg, j,k, i_gll, i
    real(kind=kreal)               :: coord(ndim,4), face_normal(3), dx_dxi(NDIM), dx_deta(NDIM), pi_2d_abg, pi_2d_xyg
    real(kind=kreal)               :: utf_dot_bkgrav, utf_dot_bkgrav_xyg, area_inv, xyg_sum
    real(kind=kreal), allocatable  :: gw(:)           ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:)
    real(kind=kreal)               :: epsilon, pi_2d, ctr, val
    integer                        :: fgdof(nndof, maxngll2d) ! face global degrees of freedom

    real(kind=kreal), allocatable  :: nodal_iceload_u(:,:), nodal_iceload_phi(:), nodal_iceload_sl(:) ! dim nnode_fs

    ! Code: 
    allocate(gw(maxngll2d))
    allocate(dshape4(2,4,maxngll2d))


    write(ICElogunit,*)'Calculating iceload...'
    iceload = zero


    if (savedata%iceload)then 
        write(ICElogunit,*) ' Saving nodal_icelog values.'
        allocate(nodal_iceload_u(ndim, nnode_fs), nodal_iceload_phi(nnode_fs), nodal_iceload_sl(nnode_fs))
        nodal_iceload_u   = zero
        nodal_iceload_phi = zero
        nodal_iceload_sl  = zero
    endif 


    call calc_iceload_epsilon(epsilon, gw, dshape4, num4, coord, face_normal, nodalicerate)

    ! Epsilon/Area
    area_inv = epsilon/SLarea

    ! Now calculate the actual iceload  - TODO there may be a way to combine the two loops
    ! and calculate epsilon, then apply it at the end, which would be more efficient
    ! but leave that til later.  
    do i_elmtfs=1,nelmt_fs

        ! GET SOME DETAILS ABOUT THE ELEMENT 
        call get_fs_details(i_elmtfs, iface, nfgll, gw, dshape4) ! face number, number of GLL on face, gauss weights, derivative of shape funcs
        num4     = gnum4_fs(:, i_elmtfs)           ! node IDs for corners 
        coord    = g_coord(:,num4)                 ! node coordinates
        gid_elmt = id_elem_fs(i_elmtfs)            ! global element ID 
        num_FS   = gnum_fs(:, i_elmtfs)            ! Global node IDs of the GLL pts on FS 

        ! DOF IDs for the nodes in matrix (property, gll pt) where property goes from 1 - 5 (ux,uy,uz,phi,theta)
        fgdof(:, 1:nfgll) = reshape(gdof(:, g_num(hexface(iface)%node, gid_elmt)),(/nndof, maxngll2d/))
        
        !write(ICElogunit,*)'  FACE_FS : ', i_elmtfs
        !do i_gll = 1, maxngll2d 
        !    write(ICElogunit,*)'     *  ', fgdof(:, i_gll)
        !enddo 
        !write(ICElogunit,*) 


        write(ICElogunit,*)' gid_elmt : ', gid_elmt
        write(ICElogunit,*)' face     : ', iface
        write(ICElogunit,*)'  fgdof : '
        do i_gll = 1, maxngll2d 
            write(ICElogunit,*)'     *  ', fgdof(:, i_gll)
        enddo 


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

            ! Calculate the coefficient for phi and u that is shared
            ! Note that u values also need to be multiplied by background gravity
            val = ( (ONE - oceanf(i_elmtfs, i_gll)) * nodalicerate(nodeid)) - area_inv*oceanf(i_elmtfs, i_gll) 
            val = val * pi_2d


            ! Displacement for direction j: + ( (1-OF)*I_dot  - epsilon/A * OF  )* pi * \nabla\Phi_j
            do j = 1, NDIM    
                dof = fgdof(j, i_gll) +1  
                if (dof.gt.0)then 
                    iceload(dof) = iceload(dof) + ( val *  grav0_nodal(j, num_FS(i_gll)) )

                    if (savedata%iceload)then 
                        nodal_iceload_u(j,nodeid) = nodal_iceload_u(j, nodeid) + (val *  grav0_nodal(j, num_FS(i_gll)) )
                    endif 
                endif  
            enddo

            ! Phi:  + ( (1-OF)*I_dot  - epsilon/A * OF  )* pi
            dof = fgdof(4, i_gll) + 1
            if (dof.gt.0)then 
                iceload(dof) = iceload(dof) + val 

                if (savedata%iceload)then 
                    nodal_iceload_phi(nodeid) = nodal_iceload_phi(nodeid) + val
                endif 
            endif 

            ! Theta: - epsilon/A * g * pi2d
            dof = fgdof(5, i_gll) + 1
            if (dof.gt.0)then 
                iceload(dof) = iceload(dof) - (area_inv * pi_2d * g0_nodal(num_FS(i_gll)))

                if (savedata%iceload)then 
                    nodal_iceload_sl(nodeid) = nodal_iceload_sl(nodeid) - (area_inv * pi_2d * g0_nodal(num_FS(i_gll)))
                endif 
            endif


        enddo  ! i_gll 
    enddo  ! i_elmtfs


    ! Multiply whole of the vector by - rho_i 
    iceload = - (rho_ice * iceload)

    ! Save iceload to file: 
    if (savedata%iceload)then 
        nodal_iceload_u   =  - nodal_iceload_u   * rho_ice
        nodal_iceload_phi =  - nodal_iceload_phi * rho_ice
        nodal_iceload_sl  =  - nodal_iceload_sl  * rho_ice

        call write_iceload_to_ensight(nodal_iceload_sl, nodal_iceload_phi, nodal_iceload_u, nodalu, i_step=0)
      
        write(ICElogunit,*)' Min value of nodal_iceload_u: ', minval(nodal_iceload_u)
        write(ICElogunit,*)' Max value of nodal_iceload_u: ', maxval(nodal_iceload_u)
        write(ICElogunit,*)
        write(ICElogunit,*)' Min value of nodal_iceload_phi: ', minval(nodal_iceload_phi)
        write(ICElogunit,*)' Max value of nodal_iceload_phi: ', maxval(nodal_iceload_phi)
        write(ICElogunit,*)
        write(ICElogunit,*)' Min value of nodal_iceload_sl: ', minval(nodal_iceload_sl)
        write(ICElogunit,*)' Max value of nodal_iceload_sl: ', maxval(nodal_iceload_sl)
        write(ICElogunit,*)'--------------------------------------------------------'
    endif
    write(ICElogunit,*)' Min value of nodalicerate: ', minval(nodalicerate)
    write(ICElogunit,*)' Max value of nodalicerate: ', maxval(nodalicerate)
    write(ICElogunit,*)
    write(ICElogunit,*)' Calculated ice load '
    write(ICElogunit,*)'  -->  Min value of iceload     : ', minscal(minval(iceload))
    write(ICElogunit,*)'  -->  Max value of iceload     : ', maxscal(maxval(iceload))

    deallocate(gw)
    deallocate(dshape4)

end subroutine calc_ice_load    






subroutine write_iceload_to_ensight(nodal_iceload_sl, nodal_iceload_phi, nodal_iceload_u, nodalu, i_step)
    use global
    use postprocess
    use free_surface
    use set_precision

    real(kind=kreal) :: nodal_iceload_sl(:), nodal_iceload_phi(:), nodal_iceload_u(:,:)
    real(kind=kreal) :: nodalu(:,:)
    integer :: i_step

        ! Save sea level iceload
    if(savedata%fsplot)then
        write(ICElogunit,*)'Saving iceload_sl'
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

    ! Displacement
    write(ICElogunit,*)'Saving iceload_u - currently using nodalu as the vector for whole mesh (instead of free surface) as a proxy...not real.'
    call write_vector_to_file(nnode,nodalu,ext='iceload_u',istep=i_step) 
    ! On the free surface
    if(savedata%fsplot)then
      call write_vector_to_file_freesurf(nnode_fs,nodal_iceload_u,&
      ext='iceload_u',istep=i_step)
    endif
    if(savedata%fsplot_plane)then
      call write_vector_to_file_freesurf(nnode_fs,nodal_iceload_u,&
      ext='iceload_u',istep=i_step,plane=.true.)
    endif
end subroutine write_iceload_to_ensight


subroutine calc_iceload_epsilon(epsilon, gw, dshape4, num4, coord, face_normal, nodalicerate)
    use free_surface
    use global
    use set_precision
    use math_constants
    ! IO: 
    integer                        :: num4(4)
    real(kind=kreal)               :: gw(:), nodalicerate(:)        
    real(kind=kreal)               :: dshape4(:,:,:)
    real(kind=kreal)               :: coord(ndim,4), face_normal(3), epsilon

    ! local vars: 
    integer :: i_elmtfs, iface, nfgll, i_gll
    real(kind=kreal) :: dx_dxi(NDIM), dx_deta(NDIM), pi_2d, val 

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
    write(ICElogunit,*)' --------------------------------------------- '
    write(ICElogunit,*)' ε value            : ', epsilon
    write(ICElogunit,*)' --------------------------------------------- '
end subroutine calc_iceload_epsilon




end module