! Last edited - W Eaton Aug 2022 
! Bilinear form matrices for rate-dependent viscoelastic loading 
! see Al-Attar and Tromp 2014

! Note that matrices are calculated via explicit summations rather than 
! matrix multiplication as done by HNG

module bilinear_form
    use set_precision
    implicit none 
    
    contains 
    
    
    
        real(kind=kreal) function kronecker(x,y)
        !Kronecker delta 
        implicit none
        integer, intent(in) :: x, y
        if (x==y)then 
            kronecker = 1.0
        else
            kronecker = 0.0
        endif     
        end function kronecker
    
    
    
    subroutine calculate_jacobian()
        ! Stores the jacobian for each element, its inverse and the 
        ! determinant to avoid multiple calculations
        ! This is fine as long as the mapping is time-invariant
    
        use global,       only: weJACINV, weDETJAC, weJAC, NDIM, nelmt, & 
                                g_coord, ngll, ngnode, g_num, nenode
        use element,      only: hex8_gnode
        use math_library, only: determinant, invert
        use integration,  only: dshape_hex8
        use set_precision
    
        implicit none 
    
        ! Local variables: 
        integer :: i_elem, igll  
        integer :: num(nenode)
        real(kind=kreal) :: coord(ngnode,NDIM), jac(NDIM,NDIM)      
    
    
        ! CODE: 
        allocate(weJAC(NDIM, NDIM, ngll, nelmt),   & 
        weJACINV(NDIM, NDIM, ngll, nelmt),&
        weDETJAC(ngll, nelmt))
    
    
        do i_elem = 1, nelmt
            ! Get global ID of nodes in element in question and coords.
            num   = g_num(:,i_elem)
            coord = transpose(g_coord(:, num(hex8_gnode)))
    
            ! For each GLL point 
            do igll = 1, ngll
                jac = matmul(dshape_hex8(:,:,igll),coord)
                weJAC(:,:,igll,i_elem)    = jac(:,:)
                weDETJAC(igll,i_elem)     = determinant(jac)
                call invert(jac)
                weJACINV(:,:,igll,i_elem) = jac
            enddo 
        enddo 
    
    
    
    
    
    end subroutine calculate_jacobian
    
    
    
    
    
    
    
    ! ______________________________________________________________________
    subroutine BF_poissons_T1(test_BF)
        ! Calculates the coefficient matrix for bilinear form Term 1
        ! due to poissons 
        ! USES 
        use set_precision
        use global, only: nelmt,ngll, g_num, ngnode, NDIM, nenode, g_coord,&
                            weJACINV, weDETJAC, nnode, Pmatglob
        use integration, only: dshape_hex8, dlagrange_gll, gll_weights
        use element, only: hex8_gnode
        use math_library,only:determinant,invert,issymmetric
        use math_constants
        implicit none 
    
        ! Input/Output variables: 
        logical :: test_BF 
    
        ! Local variables: 
        real(kind=kreal) ::  coords(3)
        real(kind=kreal) :: jacw_bars, i_sum, t1, t2, quad_sum, tf
                
        integer :: i_elem, abg, stv, bars, j,i,q, g_ind1, g_ind2,num(nenode)
    
        ! ---------------------------- CODE --------------------------------
        ! Initialise
        allocate(Pmatglob(nnode, nnode))
        Pmatglob = 0.0
        
        do i_elem = 1, nelmt ! loop elements
            num = g_num(:,i_elem)          ! global IDs for assembly
            do abg = 1, ngll               ! loop abg summation
                do stv = 1, ngll           ! loop stv summation
                    quad_sum = 0 
                    do bars = 1, ngll      ! loop quadrature summation
            
                        ! Calculate jacobian with respect to bars: 
                        jacw_bars = gll_weights(bars)* weDETJAC(bars,i_elem)
    
    
                        i_sum = 0 
                        do i = 1,3
    
                            t1 = 0 
                            do j = 1,3
                                t1 = t1 + (weJACINV(i,j, stv, i_elem) & 
                                            * dlagrange_gll(j,bars,abg))
                            enddo 
    
                            t2 = 0 
                            do q = 1,3
                                t2 = t2 + (weJACINV(i,q, stv, i_elem) & 
                                            * dlagrange_gll(q,bars,stv))
                            enddo 
    
                            i_sum = i_sum + (t1*t2)
                        enddo
                        quad_sum = quad_sum + (i_sum * jacw_bars)
                    enddo   
    
    
                    ! Our elemental matrix P(ielem, stv, abg) is then 
                    ! put directly into the assembly: 
                    ! We need the mapping for this 
                    g_ind1 = num(stv)
                    g_ind2 = num(abg)
    
                    
     
                    tf = ONE          ! VALUE IF NOT TESTING
                    if (test_BF) then ! FLAG FOR TESTING ASSEMBLY
                    ! Get test function value for tf(stv, ielem)  
                    !  --> requires coordinates 
                    ! test function = x + y + 2z for Test Case 1
                        coords =  g_coord(:,g_ind1)
                        tf     =  coords(1) + coords(2) + TWO*coords(3)
                    endif 
                    
                    ! Add to global matrix 
                    Pmatglob(g_ind1, g_ind2) = Pmatglob(g_ind1, g_ind2) &
                                               +  (quad_sum*tf)  
                enddo                       
            enddo  
        enddo 
    
    end subroutine BF_poissons_T1
    ! ______________________________________________________________________
    
    
    ! ______________________________________________________________________
    subroutine BF_bulk_T1(test_BF)
        ! Calculates the global matrix for 2nd term in BF (incl. kappa)
        use global, only: Kmatglob, bulkmod_elmt, nenode, ngll, nelmt, & 
                          nnode, g_num, weJACINV, weDETJAC, g_coord
        use integration, only: dlagrange_gll, gll_weights
        use set_precision
        use math_constants
        implicit none 
    
        ! IO variables: 
        logical :: test_BF
        ! Local variables: 
        integer :: i_elem, i, j, abg, stn, bars, gi1, gi2, q
        integer :: num(nenode)
        real(kind=kreal) :: quad_sum, Li, Lj, k, jacw, tf(3), coords(3)
        !real(kind=kreal) :: Kmatlocal(nelmt, 3*nnode, 3*nnode) !temp 
        real(kind=kreal) :: kappa(ngll) !temp 
        ! ---------------------------- CODE --------------------------------
        ! Initialise matrix
        allocate(Kmatglob(3*nnode, 3*nnode))
        Kmatglob = 0.0
    
        ! Loop through elements: 
        do i_elem = 1, nelmt 
            num    = g_num(:,i_elem)          ! global IDs for assembly
            kappa  = bulkmod_elmt(:,i_elem)   ! Bulk modulus for element
    
            do abg = 1,ngll 
                do i = 1,3
                    do stn = 1,ngll 
                        do j = 1,3 
                            
                            quad_sum = 0.0
                            do bars = 1,ngll 
    
                                Li = 0.0 ; Lj = 0.0 ! initialise
                                do q=1,3
                                Li = Li + weJACINV(i,q, bars, i_elem) * & 
                                          dlagrange_gll(q, bars, abg)
                                Lj = Lj + weJACINV(j,q, bars, i_elem) * & 
                                          dlagrange_gll(q, bars, stn)
                                enddo 
    
                                jacw = weDETJAC(bars, i_elem) * & 
                                       gll_weights(bars)
                    
                                ! If testing then kappa becomes xy/z: 
                                if (test_BF) then 
                                    ! Get coordinates of bars GLL point
                                    coords =  g_coord(:,num(bars))
                                    k = coords(2)*coords(1)/coords(3)
                                else
                                    k = kappa(bars)
                                endif 
    
    
                                quad_sum=quad_sum + (jacw*k*Li*Lj)
                            enddo ! bars
    
                            ! Global matrix indices
                            gi1 = 3*(num(stn)-1) + j 
                            gi2 = 3*(num(abg)-1) + i
    
                            tf = ONE ! Test vector - 1 unless testing. 
                            if (test_BF) then
                                ! Running test on this assembly so need the 
                                ! tf vector to be (4y, 2yz, 3xz)
                                ! where x,y,z are for GLL point stn 
                                coords  =  g_coord(:,num(stn))
                                tf(1)   =  FOUR  * coords(2)
                                tf(2)   =  TWO   * coords(2) * coords(3)
                                tf(3)   =  THREE * coords(1) * coords(3)
                            endif 
    
                            Kmatglob(gi1,gi2) = Kmatglob(gi1, gi2)  & 
                                                + (tf(j)*quad_sum)
    
                        enddo ! j
                    enddo ! stn
                enddo ! i
            enddo ! abg 
        enddo  ! elmt
    
    end subroutine BF_bulk_T1
    ! ______________________________________________________________________
    
    
    subroutine BF_backgrav_1(test_BF)
        ! Uses
        use global
        use math_constants
        use set_precision
        use integration
    
        implicit none 
        ! IO variables
        logical :: test_BF
        ! Local variables: 
        integer :: i_elem, i, j, abg, stn, i_bars, k , gi1, gi2 
        integer :: num(nenode)
    
        real(kind=kreal) :: quad_sum, rho, jacw, ksum, tf_stnj(3), coords(3), & 
                            g0_local(3), g0
                            
        ! Code: 
    
        allocate(B1glob(3*nnode, 3*nnode))
        B1glob = ZERO
    
    
        ! Loop through elements: 
        do i_elem = 1, nelmt
            num   = g_num(:,i_elem)          ! global IDs for assembly
            
            do abg=1,ngll 
                do i=1,3
                    do stn=1,ngll
    
                        ! If testing the function: 
                        if (test_BF) then 
                            ! Calc test function
                            coords  =  g_coord(:,num(stn))
                            tf_stnj(1) = 2*coords(1)
                            tf_stnj(2) = 3*coords(2)
                            tf_stnj(3) = coords(3)
                        
                            rho   = coords(1) + coords(2)
    
                            ! Background gravity - USES coords of ABG GLL point
                            coords  =  g_coord(:,num(abg))
                            g0_local(1) = 1.0
                            g0_local(2) = 0
                            g0_local(3) = 2.0*coords(3)
    
                        else
                            tf_stnj(:)  = 1.0
                            g0_local(:) = 1.0
                            rho         = massdens_elmt(stn, i_elem)
                            write(*,*)'Using background grav = 1'
                        endif 
    
    
                        do j=1,3
                            
                            g0    = g0_local(i)
                            jacw  = weDETJAC(stn, i_elem) * gll_weights(stn)
    
                            ksum = 0 
                            do k=1,3
                                ksum = ksum + (weJACINV(j,k,stn, i_elem) &
                                              * dlagrange_gll(k,stn,abg))
                            enddo 
    
    
                            ! Global matrix indices
                            gi1 = 3*(num(stn)-1) + j 
                            gi2 = 3*(num(abg)-1) + i
    
                            ! Add to global matrix: 
                            B1glob(gi1, gi2) = B1glob(gi1, gi2) +  tf_stnj(j)*(rho*g0*jacw*ksum)
    
    
                        enddo !j
                    enddo !stn
                enddo !i
            enddo !abg
        enddo !ielem
    
    end subroutine BF_backgrav_1
    
    
    
    
    
    
    subroutine BF_backgrav_2(test_BF)
        ! Uses
        use global
        use math_constants
        use set_precision
        use integration
    
        implicit none 
        ! IO variables
        logical :: test_BF
        ! Local variables: 
        integer :: i_elem, i, j, abg, stn, i_bars, k , gi1, gi2 
        integer :: num(nenode)
    
        real(kind=kreal) :: quad_sum, rho, jacw, ksum, tf_stnj(3), coords(3), & 
                            g0_local(3), g0
                            
        ! Code: 
    
        allocate(B2glob(3*nnode, 3*nnode))
        B2glob = ZERO
    
    
        ! Loop through elements: 
        do i_elem = 1, nelmt
            num   = g_num(:,i_elem)          ! global IDs for assembly
            
            do abg=1,ngll 
                do i=1,3
                    do stn=1,ngll
    
                        ! If testing the function: 
                        if (test_BF) then 
                            coords  = g_coord(:,num(abg))
                            rho     = coords(1) + coords(2)
    
                            ! Background gravity - USES coords of STN GLL point
                            coords      =  g_coord(:,num(stn))
                            g0_local(1) = 1.0
                            g0_local(2) = 0
                            g0_local(3) = 2.0*coords(3)
                            ! Calc test function
                            tf_stnj(1)  = 2*coords(1)
                            tf_stnj(2)  = 3*coords(2)
                            tf_stnj(3)  = coords(3)
    
                        else
                            tf_stnj(:)  = 1.0
                            g0_local(:) = 1.0
                            rho         = massdens_elmt(abg, i_elem)
                            write(*,*)'Using background grav = 1'
                        endif 
    
    
                        do j=1,3
                            g0    = g0_local(j)
                            jacw  = weDETJAC(abg, i_elem) * gll_weights(abg)
    
                            ksum = 0 
                            do k=1,3
                                ksum = ksum + (weJACINV(i,k,abg, i_elem) &
                                               * dlagrange_gll(k,abg,stn))
                            enddo 
    
    
                            ! Global matrix indices
                            gi1 = 3*(num(stn)-1) + j
                            gi2 = 3*(num(abg)-1) + i
    
                            ! Add to global matrix: 
                            B2glob(gi1, gi2) = B2glob(gi1, gi2) + (tf_stnj(j)*rho*g0*jacw*ksum)
    
                        enddo !j
                    enddo !stn
                enddo !i
            enddo !abg
        enddo !ielem
    
    end subroutine BF_backgrav_2
    
    
    
    subroutine BF_backgrav_3(test_BF)
        use global 
        use integration 
        use set_precision
        use math_constants
        implicit none 
    
        ! IO variables
        logical :: test_BF
    
        ! Local variables: 
        integer :: i_elem, num(nenode), abg, i, stn, j, bars, k, row, col 
        real(kind=kreal) :: coords(3), g0_abg(3), g0_stn(3), tf_u(3), & 
                            quad_sum, rho, t1, t2, jacw, kd1, kd2
    
    
        ! CODE: 
    
    
        allocate(B3(4*nnode, 4*nnode))
        B3 = 0.0
    
        do i_elem = 1, nelmt 
            num = g_num(:,i_elem)
            do abg = 1, ngll 
                do i=1,3
                    do stn = 1, ngll 
    
                        ! Test function
                        if (test_BF) then 
                            ! STN coords
                            coords       = g_coord(:, num(stn))
                            tf_u(1)      = 2*coords(1)
                            tf_u(2)      = 3*coords(2) 
                            tf_u(3)      =   coords(3)
                            
                            ! Background gravity
                            g0_stn(1)    = 1
                            g0_stn(2)    = 0
                            g0_stn(3)    = 2*coords(3)
    
                            ! ABG coords
                            coords       = g_coord(:, num(abg))
                            g0_abg(1)    = 1
                            g0_abg(2)    = 0
                            g0_abg(3)    = 2*coords(3)
                        else 
                            tf_u(:)    = 1.0
                            g0_abg(:)  = 1.0
                            g0_stn(:)  = 1.0
                        endif 
    
    
                        do j=1,3
            
                            quad_sum = 0.0 
                            do bars=1,ngll
                            
                                jacw = weDETJAC(bars,i_elem)*gll_weights(bars)
                                
                                if (test_BF) then
                                    coords = g_coord(:, num(bars))
                                    rho    = coords(1) + coords(2)
                                else
                                    rho = massdens_elmt(bars, i_elem)
                                endif
                            
                                ! Internal sums: 
                                t1 = 0.0
                                t2 = 0.0 
    
                                do k=1,3
                                    t1 = t1 + (weJACINV(i,k,bars,i_elem) * &
                                               dlagrange_gll(k,bars,abg))
    
                                    t2 = t2 + (weJACINV(j,k,bars,i_elem) * &
                                               dlagrange_gll(k,bars,stn))           
                                enddo 
    
                              
                                 kd1 =  kronecker(bars,stn)
                                 kd2 =  kronecker(bars,abg)
    
    
                                quad_sum = quad_sum + & 
                                           (jacw*rho * & 
                                           (  (g0_stn(j)*t1*kd1) & 
                                            + (g0_abg(i)*t2*kd2)) )
                            
                            enddo ! bars
                          
                            
                            ! Assembly: 
                            row = (num(stn) -1)*3 + j    ! Row index
                            col = (num(abg) -1)*3 + i    ! Column index
    
                            B3(row,col) = B3(row,col) + (quad_sum*tf_u(j))
    
                        enddo ! j
                    enddo ! stn
                enddo ! i
            enddo !abg
        enddo !elem
    
    
    end subroutine BF_backgrav_3 
    
    
    
    
    
    subroutine BF_disp_grav_coupling(test_BF)
        ! Calculates the 6th Term of the BF with the coupling between 
        ! gravity and displacement 
    
        use global 
        use integration 
        use set_precision
        use math_constants
        implicit none 
    
        ! IO variables: 
        logical :: test_BF
    
        ! Local variables: 
        integer :: i_elem, u_dim,  abg, stn, i, j, rd, cd, rp, cp
        integer :: num(nenode)
    
        real(kind=kreal) :: phi_sum, tf_u_i_stn(3), rho_phi, jacw_phi, L_phi, & 
                            tf_phi_stn, rho_disp, jacw_disp, L_disp, & 
                            disp_sum, coords(3)
        ! Code: 
    
    
    
        ! Initialise the coupling matrix 
        allocate(DCG_global(nnode*4, nnode*4))
        DCG_global = 0.0 
    
        u_dim = 3*nnode
    
    
        ! Default test function values 
        tf_u_i_stn(:) = 1.0      ! Disp. test function
        tf_phi_stn    = 1.0      ! Phi test function
    
        do i_elem = 1, nelmt
            num = g_num(:,i_elem)
    
            do abg = 1, ngll 
                do stn = 1, ngll 
    
                    phi_sum = 0.0
                    rho_phi   = massdens_elmt(stn,i_elem)
                    rho_disp  = massdens_elmt(abg,i_elem)
    
                    ! Test functions if testing: 
                    if (test_BF) then 
                        tf_u_i_stn(1) = 4
                        tf_u_i_stn(2) = 1
                        tf_u_i_stn(3) = 2
    
                        ! STN coords
                        coords        = g_coord(:, g_num(stn, i_elem))
                        tf_phi_stn    = 2*coords(1)*coords(2)
                        rho_phi       = coords(1) + coords(2) - coords(3)
                        ! ABG coords: 
                        coords        = g_coord(:, g_num(abg, i_elem))
                        rho_disp      = coords(1) + coords(2) - coords(3)
                    endif 
    
                    do i=1,3
                        do j=1,3
    
                            ! Gravity component:
    
                            jacw_phi = weDETJAC(stn,i_elem)*gll_weights(stn)
                            L_phi    = ( weJACINV(i,j,stn,i_elem) & 
                                         * dlagrange_gll(j,stn,abg))
                            phi_sum = (phi_sum + &
                                       (tf_u_i_stn(i)*rho_phi*jacw_phi*L_phi) )            
    
                            ! Displacement component: 
    
                            jacw_disp = weDETJAC(abg,i_elem)*gll_weights(abg)
                            L_disp    = weJACINV(i,j,abg,i_elem) * &
                                        dlagrange_gll(j,abg,stn)
                            disp_sum  = tf_phi_stn*rho_disp*jacw_disp*L_disp
    
               
                            ! Assemble the displacement part of matrix: 
                            rd = (num(stn) -1)*3 + j    ! Row index
                            cd = (num(abg) -1)*3 + i    ! Column index
    
                            DCG_global(rd, cd) = DCG_global(rd, cd)+disp_sum 
    
                        enddo 
                    enddo 
    
                    ! Assemble the gravity part of matrix:
                    rp = num(stn) + u_dim
                    cp = num(abg) + u_dim
                    DCG_global(rp, cp) = DCG_global(rp, cp)+phi_sum 
    
    
                   
                enddo ! stn
            enddo ! abg
        enddo ! elem
    
    
    end subroutine BF_disp_grav_coupling
    
    
    
    ! NEED TO ADD: STRAIN DEVIATOR TERM 
    
    subroutine BF_strain_deviator(test_BF)
        ! USES 
        use global 
        use integration 
        use set_precision
        use math_constants
        implicit none 
    
        ! IO variables: 
        logical :: test_BF
    
        ! Local variables: 
        integer          :: num(ngll)
        integer          :: i_elem, abg, lam, stn, gam, bars, i, j, q, g, &
                            row, col 
        real(kind=kreal) :: sum, Ajq, Air, Apg, Bjq, Bir, Bpg, Agsum, Bgsum,&
                            A_ij_gam, B_ij_lam, mu, jacw, tf(3), coords(3)
    
        ! Allocate the strain deviator matrix, D: 
        allocate(Dmat(4*nnode, 4*nnode))
    
    
    
        ! Begin the plethora of loops!! 
        do i_elem = 1, nelmt
    
            ! Get the IDs for the element 
            num = g_num(:,i_elem)
    
            do abg = 1, ngll
                do lam = 1, 3
                    do stn = 1, ngll
                        do gam = 1, 3
    
                            sum = 0.0 ! initialise sum 
                            do bars = 1, ngll 
                                do i = 1,3
                                    do j = 1,3
    
                                        ! Mu and integration weightings
                                        if (test_BF) then 
                                            coords = g_coord(:, g_num(bars, i_elem))
                                            mu = coords(1)+coords(2)+coords(3) 
                                            !write(*,*)'mu: ', mu
                                        else
                                            mu = shearmod_elmt(bars, i_elem)
                                        endif
    
                                        jacw = weDETJAC(bars,i_elem) * &
                                               gll_weights(bars)
    
                                        ! Now calculate A_ijgamma, B_ijgamma
                                        Ajq = 0.0; Air = 0.0; Apg = 0.0;
                                        Bjq = 0.0; Bir = 0.0; Bpg = 0.0;
                                        
                                        do q = 1,3
                                            Ajq = Ajq + weJACINV(j,q,bars,i_elem)*dlagrange_gll(q,bars,stn)
                                            Bjq = Bjq + weJACINV(j,q,bars,i_elem)*dlagrange_gll(q,bars,abg)
    
                                            Air = Air + weJACINV(i,q,bars,i_elem)*dlagrange_gll(q,bars,stn)
                                            Bir = Bir + weJACINV(i,q,bars,i_elem)*dlagrange_gll(q,bars,abg)
    
                                            Agsum = 0.0 ; Bgsum = 0.0;
    
                                            do g = 1, 3
                                                Agsum = Agsum + weJACINV(q,g,bars,i_elem)*dlagrange_gll(g,bars,stn)
                                                Bgsum = Bgsum + weJACINV(q,g,bars,i_elem)*dlagrange_gll(g,bars,abg)
                                            enddo !g 
    
                                            Apg = Apg + (Agsum*kronecker(q,gam) * TWO_THIRD * kronecker(i,j))
                                            Bpg = Bpg + (Bgsum*kronecker(q,lam) * TWO_THIRD * kronecker(i,j))
                                        enddo ! q 
    
    
                                        A_ij_gam = (Ajq*kronecker(gam,i)) + (Air*kronecker(gam,j)) - Apg
                                        B_ij_lam = (Bjq*kronecker(lam,i)) + (Bir*kronecker(lam,j)) - Bpg
    
                                        sum = sum + mu*jacw*A_ij_gam*B_ij_lam
    
    
                                    enddo ! j
                                enddo ! i
                            enddo !bars
                            
    
                            ! ASSEMBLEY 
                            row = (num(stn) -1)*3 + gam    ! Row index
                            col = (num(abg) -1)*3 + lam    ! Column index
    
    
                            ! Test function if testing:
                            if (test_BF) then 
                                ! Calculate test function 
                                coords = g_coord(:, g_num(stn, i_elem))
    
                                tf(1) = coords(3) 
                                tf(2) = coords(2)
                                tf(3) = coords(1)
                            else
                                tf = 1.0 ! 3x1 array
                            endif
                            
                            !write(*,*)'tf: ', tf(gam)
                            Dmat(row, col) = Dmat(row, col) + HALF*sum*tf(gam) 
                                              
                        enddo !gam 
                    enddo !stn 
                enddo !lam 
            enddo !abg 
    
    
        enddo !i_elem 
    
        
    
    
    
    end subroutine BF_strain_deviator
    
    
    
    
    
    
    
    
        !! ____________________ TESTS ______________________________________
    
        subroutine BF_poissons_T1_test_eval()
        ! Calculates a global phi and completes the integration for 
        ! Test case 1 of BF_poisson 
        use global
        use set_precision
        use dimensionless
        implicit none 
    
    
        real(kind=kreal) :: Phi_global(nnode), vector(nnode), coords(3)
        real(kind=kreal) :: total_sum, analytical, x1, x2, y1, y2, z1, &
                            z2, diff 
        integer :: i, igll, num, ielem
    
    
    
        Phi_global = 0.0
    
        ! Get the phi values at each nodal point. These are given by 
        ! xy + zx + y^2 
        do ielem=1,nelmt
            do igll = 1, ngll
                num = g_num(igll, ielem)
                coords = g_coord(:, num)
    
                Phi_global(num) =  coords(1)*coords(2)  & 
                                    + coords(3)*coords(1)  &
                                    + coords(2)*coords(2)
            enddo 
        enddo 
    
    
        ! We now have the global coefficient matrix Pmatglob and 
        ! the Phi values: 
        vector(:) = matmul(Pmatglob, Phi_global)
        total_sum = 0 
        do i = 1, nnode
            total_sum = total_sum + vector(i)
        enddo 
    
        ! Now calculating the analytical solution for Test 1: 
        ! Coordinates are init. in global and calculated in driver code
        x1 = model_minx*NONDIM_L;  x2 = model_maxx*NONDIM_L
        y1 = model_miny*NONDIM_L;  y2 = model_maxy*NONDIM_L
        z1 = model_minz*NONDIM_L;  z2 = model_maxz*NONDIM_L
    
        analytical =  1.5 * ((x2**2) - (x1**2)) * (y2-y1) * (z2-z1) &
                    + 1.5 * (x2-x1) * ((y2**2) - (y1**2)) * (z2-z1) &
                    + 0.5 * (x2-x1) * (y2-y1) * ((z2**2) - (z1**2))
        diff = total_sum - analytical 
    
            write(*,*)'precision           :  ', BA_TEST_PRECISION
            write(*,*)'Analytical integral :  ', analytical
            write(*,*)'SEM Integral        : ', total_sum
            write(*,*)'difference          : ', diff
            write(*,*)'percentage error    : ', diff*100/analytical,' %'
        if (diff.gt.BA_TEST_PRECISION) then 
            write(*,*)'!!!!! Test failed: BF_poissons_T1_test_eval !!!!!!'
            stop
        else
            write(*,*)'PASSED: BF_poissons_T1_test_eval'
        endif 
    
    end subroutine BF_poissons_T1_test_eval
    
    
    
    
    
    
        subroutine BF_bulk_T2_test_eval()
            ! Calculates a global Bulk_k and completes the integration for 
            ! Test case 1 of BF_bulk
            use global
            use set_precision
            use dimensionless
            use math_constants
            implicit none 
        
        
            real(kind=kreal) :: U_global(nnode*3), vector(nnode*3), coords(3)
            real(kind=kreal) :: total_sum, analytical, x1, x2, y1, y2, z1, &
                                z2, diff, xx, xxx, yy, yyy, & 
                                zz, zzz
            integer :: i, igll, num, ielem, ind 
        
        
        
            U_global = 0.0
        
            ! Get the disp values at each nodal point. These are given by 
            ! (xyz, yz^2, 2x+2y )
            do ielem=1,nelmt
                do igll = 1, ngll
                    num = g_num(igll, ielem)
                    coords = g_coord(:, num)
    
                    ind = (num-1)*3 + 1 
        
                    U_global(ind  ) =  coords(1) * coords(2) * coords(3)
                    U_global(ind+1) =  coords(2) * (coords(3)**2)
                    U_global(ind+2) =  TWO*(coords(1) + coords(2)) 
    
                enddo 
            enddo 
        
        
            ! We now have the global coefficient matrix Kmatglob and 
            ! the U values: 
            vector(:) = matmul(Kmatglob, U_global)
            total_sum = 0 
            do i = 1, nnode*3 
                total_sum = total_sum + vector(i)
            enddo 
        
        
            ! Now calculating the analytical solution for Test 1: 
        
            ! Coordinates are init. in global and calculated in driver code
            
            x1 = model_minx*NONDIM_L;  x2 = model_maxx*NONDIM_L
            y1 = model_miny*NONDIM_L;  y2 = model_maxy*NONDIM_L
            z1 = model_minz*NONDIM_L;  z2 = model_maxz*NONDIM_L
            
            ! Calc. squares and cubes of each term 
            xx = (x2**2) - (x1**2)  ;  xxx = (x2**3) - (x1**3)
            yy = (y2**2) - (y1**2)  ;  yyy = (y2**3) - (y1**3)
            zz = (z2**2) - (z1**2)  ;  zzz = (z2**3) - (z1**3)
    
            
            
            analytical = ONE_SIXTH*xx*yyy*zz + ONE_THIRD*xxx*yyy*(z2-z1) + & 
                         ONE_SIXTH*xx*yy*zzz + ONE_QUARTER*xxx*yy*zz
        
        
            diff = total_sum - analytical 
        
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*)'precision           :  ', BA_TEST_PRECISION
            write(*,*)'Analytical integral :  ', analytical
            write(*,*)'SEM Integral        : ', total_sum
            write(*,*)'difference          : ', diff
            write(*,*)'percentage error    : ', diff*100/analytical,' %'
            if (diff.gt.BA_TEST_PRECISION) then 
                write(*,*)'!!!!! Test failed: BF_bulk_T2_test_eval !!!!!!'
                stop
            else
                write(*,*)'PASSED: BF_bulk_T2_test_eval'
            endif 
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
    
        
        end subroutine BF_bulk_T2_test_eval
    
    
    
    
        subroutine BF_bg1_T3_test_eval()
            ! Calculates a global background grav T3 and completes the integration for 
            ! Test case 1 of BF_backgroundgrav1 
            use global
            use set_precision
            use dimensionless
            use math_constants
            implicit none 
        
        
            real(kind=kreal) :: U_global(nnode*3), vector(nnode*3), coords(3)
            real(kind=kreal) :: total_sum, analytical, x1, x2, y1, y2, z1, &
                                z2, diff, XN(4), YN(4), ZN(4)
            integer :: i, igll, num, ielem, ind 
        
            U_global = 0.0
            ! Get the velocity values at each nodal point. These are given by 
            ! (2x+y, 3xyz, y^2z)
            do ielem=1,nelmt
                do igll = 1, ngll
                    num    = g_num(igll, ielem)
                    coords = g_coord(:, num)
    
                    ind = (num-1)*3 + 1 
        
                    U_global(ind  ) =  2*coords(1) + coords(2) 
                    U_global(ind+1) =  3*coords(1)*coords(2)*coords(3)
                    U_global(ind+2) =  coords(2)*coords(2)*coords(3)
                enddo 
            enddo 
        
        
            ! We now have the global coefficient matrix Kmatglob and 
            ! the U values: 
            vector(:) = matmul(B1glob, U_global)
            total_sum = 0 
            do i = 1, nnode*3 
                total_sum = total_sum + vector(i)
            enddo 
        
        
            ! Now calculating the analytical solution for Test 1: 
            ! Coordinates are init. in global and calculated in driver code
            x1 = model_minx*NONDIM_L;  x2 = model_maxx*NONDIM_L
            y1 = model_miny*NONDIM_L;  y2 = model_maxy*NONDIM_L
            z1 = model_minz*NONDIM_L;  z2 = model_maxz*NONDIM_L
            
            ! Calc. squares and cubes of each term 
            do i=1,4
                XN(i)  = (x2**i) - (x1**i)
                YN(i)  = (y2**i) - (y1**i)
                ZN(i)  = (z2**i) - (z1**i)
    
            enddo 
    
            
            analytical = (FOUR_THIRD     * XN(3) * YN(1) * ZN(1)) + & 
                         (SEVEN_QUARTERS * XN(2) * YN(2) * ZN(1)) + &
                         (EIGHT_NINTHS   * XN(2) * YN(3) * ZN(3)) + &
                         (                 XN(1) * YN(3) * ZN(1)) + &
                         (FOUR_THIRD     * XN(1) * YN(4) * ZN(3))
        
        
            diff = total_sum - analytical 
        
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*)'precision           :  ', BA_TEST_PRECISION
            write(*,*)'Analytical integral :  ', analytical
            write(*,*)'SEM Integral        : ', total_sum
            write(*,*)'difference          : ', diff
            write(*,*)'percentage error    : ', diff*100/analytical,' %'
            if (diff.gt.BA_TEST_PRECISION) then 
                write(*,*)'!!!!! Test failed: BF_bg1_T3_test_eval !!!!!!'
                stop
            else
                write(*,*)'PASSED: BF_bg1_T3_test_eval'
            endif 
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
    
        
        end subroutine BF_bg1_T3_test_eval
    
    
    
    
    
    
    
        subroutine BF_bg2_T4_test_eval()
            ! Calculates a global backgrav B2 matrix and tests 
            use global
            use set_precision
            use dimensionless
            use math_constants
            implicit none 
        
        
            real(kind=kreal) :: U_global(nnode*3), vector(nnode*3), coords(3)
            real(kind=kreal) :: total_sum, analytical, x1, x2, y1, y2, z1, &
                                z2, diff, XN(4), YN(4), ZN(4)
            integer :: i, igll, num, ielem, ind 
        
            U_global = 0.0
    
            ! Get the velocity values at each nodal point. These are given by 
            ! (2x+y, 3xyz, y^2z)
            do ielem=1,nelmt
                do igll = 1, ngll
                    num    = g_num(igll, ielem)
                    coords = g_coord(:, num)
    
                    ind = (num-1)*3 + 1 
        
                    U_global(ind  ) =  2*coords(1) + coords(2) 
                    U_global(ind+1) =  3*coords(1)*coords(2)*coords(3)
                    U_global(ind+2) =  coords(2)*coords(2)*coords(3)
                enddo 
            enddo 
        
        
            ! We now have the global coefficient matrix Kmatglob and 
            ! the U values: 
            vector(:) = matmul(B2glob, U_global)
            total_sum = 0 
            do i = 1, nnode*3 
                total_sum = total_sum + vector(i)
            enddo 
        
        
            ! Now calculating the analytical solution for Test 1: 
            ! Coordinates are init. in global and calculated in driver code
            x1 = model_minx*NONDIM_L;  x2 = model_maxx*NONDIM_L
            y1 = model_miny*NONDIM_L;  y2 = model_maxy*NONDIM_L
            z1 = model_minz*NONDIM_L;  z2 = model_maxz*NONDIM_L
            
            ! Calc. squares and cubes of each term 
            do i=1,4
                XN(i)  = (x2**i) - (x1**i)
                YN(i)  = (y2**i) - (y1**i)
                ZN(i)  = (z2**i) - (z1**i)
            enddo 
    
    
            
            analytical = ((FOUR_THIRD  * XN(3) * YN(1) * ZN(1)) + & 
                          (TWO_NINTHS  * XN(2) * YN(3) * ZN(3)) + &
                          (THREE_HALFS * XN(2) * YN(2) * ZN(1)) + &
                          (ONE_THIRD   * XN(1) * YN(4) * ZN(3)) + &
                          (TWO_THIRD   * XN(1) * YN(3) * ZN(1)))
    
        
            diff = total_sum - analytical 
        
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*)'precision           :  ', BA_TEST_PRECISION
            write(*,*)'Analytical integral :  ', analytical
            write(*,*)'SEM Integral        : ', total_sum
            write(*,*)'difference          : ', diff
            write(*,*)'percentage error    : ', diff*100/analytical,' %'
            if (diff.gt.BA_TEST_PRECISION) then 
                write(*,*)'!!!!! Test failed: BF_bg2_T4_test_eval !!!!!!'
                stop
            else
                write(*,*)'PASSED: BF_bg2_T4_test_eval'
            endif 
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
    
        
        end subroutine BF_bg2_T4_test_eval
    
    
    
    
    
        subroutine BF_bg3_T5_test_eval()
            ! Calculates a global backgrav B3 matrix and tests 
            use global
            use set_precision
            use dimensionless
            use math_constants
            implicit none 
        
        
            real(kind=kreal) :: vector(nnode*4), coords(3), DOF_global(nnode*4)
            real(kind=kreal) :: total_sum, analytical, x1, x2, y1, y2, z1, &
                                z2, diff, X(4), Y(4), Z(4)
            integer :: i, igll, num, ielem, ind, j
    
            DOF_global = 0.0
    
    
            do ielem=1,nelmt
                do igll = 1, ngll
                    
                    num    = g_num(igll, ielem)
                    coords = g_coord(:, num)
    
                    ind = (num-1)*3
        
                    DOF_global(ind+1) =  2*coords(1) + coords(2) 
                    DOF_global(ind+2) =  3*coords(1)*coords(2)*coords(3)
                    DOF_global(ind+3) =  coords(2)*coords(2)*coords(3)
    
                enddo 
            enddo 
    
    
    
      
    
    
        
            ! We now have the global coefficient matrix Kmatglob and 
            ! the U values: 
            vector(:) = matmul(B3, DOF_global)
            total_sum = 0 
            do i = 1, nnode*4
                total_sum = total_sum + vector(i)
            enddo 
        
        
            ! Now calculating the analytical solution for Test 1: 
            ! Coordinates are init. in global and calculated in driver code
            x1 = model_minx*NONDIM_L;  x2 = model_maxx*NONDIM_L
            y1 = model_miny*NONDIM_L;  y2 = model_maxy*NONDIM_L
            z1 = model_minz*NONDIM_L;  z2 = model_maxz*NONDIM_L
            
            ! Calc. squares and cubes of each term 
            do i=1,4
                X(i)  = (x2**i) - (x1**i)
                Y(i)  = (y2**i) - (y1**i)
                Z(i)  = (z2**i) - (z1**i)
            enddo 
    
    
            analytical = ( (EIGHT*X(3)*Y(3)*Z(1) + 192.0_kreal*X(3)*Y(1)*Z(1) + &
                            NINE*X(2)*Y(4)*Z(1)  + 198.0_kreal*X(2)*Y(2)*Z(1) + &
                            72.0_kreal*X(1)*Y(3)*Z(1) )/36.0_kreal  + & 
                           (THREE*X(4)*Y(1)*Z(2) + TWO*X(3)*Y(2)*Z(2))/FOUR& 
                           + ((FOUR*X(3)*Y(1)*Z(4)  + &
                           THREE*X(2)*Y(2)*Z(4)))/EIGHT + &
                           ( 14.0_kreal*X(2)*Y(3)*Z(3)  + 12.0_kreal*X(2)*Y(1)*Z(3)  + &
                           21.0_kreal*X(1)*Y(4)*Z(3)  + 12.0_kreal*X(1)*Y(2)*Z(3) )/18.0_kreal)
    
            diff = total_sum - analytical 
        
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*)'precision           :  ', BA_TEST_PRECISION
            write(*,*)'Analytical integral :  ', analytical
            write(*,*)'SEM Integral        : ', total_sum
            write(*,*)'difference          : ', diff
            write(*,*)'percentage error    : ', diff*100/analytical,' %'
            if (diff.gt.BA_TEST_PRECISION) then 
                write(*,*)'!!!!! Test failed: BF_bg3_T5_test_eval !!!!!!'
                stop
            else
                write(*,*)'PASSED: BF_bg3_T5_test_eval'
            endif 
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
    
    
        end subroutine BF_bg3_T5_test_eval 
    
    
    
    
        subroutine BF_coupled_T6_test_eval()
            ! Calculates a grav disp coupling matrix and tests 
            use global
            use set_precision
            use dimensionless
            use math_constants
            implicit none 
        
        
            real(kind=kreal) :: DOF_global(nnode*4), vector(nnode*4), coords(3)
            real(kind=kreal) :: total_sum, analytical, x1, x2, y1, y2, z1, &
                                z2, diff, X(4), Y(4), Z(4)
            integer :: i, igll, num, ielem, ind , dim_u
        
            DOF_global = 0.0
    
            dim_u = 3*nnode
            
            do ielem=1,nelmt
                do igll = 1, ngll
                    num    = g_num(igll, ielem)
                    coords = g_coord(:, num)
    
                    ind = (num-1)*3 
        
                    ! U dot 
                    DOF_global(ind+1) =  coords(1)
                    DOF_global(ind+2) =  3*coords(2)*coords(3)
                    DOF_global(ind+3) =  coords(3)*coords(3)
    
                    ! Phi dot
                    DOF_global(num + dim_u) = coords(1)+coords(2)+coords(3)  
                enddo 
            enddo 
        
      
        
            ! We now have the global coefficient matrix DGC and 
            ! the U+Phi (DOF) values: 
            vector(:) = matmul(DCG_global, DOF_global)
            total_sum = 0 
            do i = 1, nnode*4
                total_sum = total_sum + vector(i)
            enddo 
        
        
            ! Now calculating the analytical solution for Test 1: 
            ! Coordinates are init. in global and calculated in driver code
            x1 = model_minx*NONDIM_L;  x2 = model_maxx*NONDIM_L
            y1 = model_miny*NONDIM_L;  y2 = model_maxy*NONDIM_L
            z1 = model_minz*NONDIM_L;  z2 = model_maxz*NONDIM_L
            
            ! Calc. squares and cubes of each term 
            do i=1,4
                X(i)  = (x2**i) - (x1**i)
                Y(i)  = (y2**i) - (y1**i)
                Z(i)  = (z2**i) - (z1**i)
            enddo 
    
    
            
            analytical = (HALF        * X(3)*Y(2)*Z(2)   +   &
                          ONE_THIRD   * X(3)*Y(2)*Z(1)   +   &
                          HALF        * X(2)*Y(3)*Z(2)   +   &
                          ONE_THIRD   * X(2)*Y(3)*Z(1)   -   &
                          HALF        * X(2)*Y(2)*Z(3)   -   &
                          ONE_QUARTER * X(2)*Y(2)*Z(2)   +   & 
                          (SEVEN_TWO  *   (X(2)*Y(1)*Z(1)+   &
                                           X(1)*Y(2)*Z(1)-   &
                                           X(1)*Y(1)*Z(2)))  &
                          )
    
        
            diff = total_sum - analytical 
        
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*)'precision           :  ', BA_TEST_PRECISION
            write(*,*)'Analytical integral :  ', analytical
            write(*,*)'SEM Integral        : ', total_sum
            write(*,*)'difference          : ', diff
            write(*,*)'percentage error    : ', diff*100/analytical,' %'
            if (diff.gt.BA_TEST_PRECISION) then 
                write(*,*)'!!!!! Test failed: BF_coupled_T6_test_eval !!!!!!'
                stop
            else
                write(*,*)'PASSED: BF_coupled_T6_test_eval'
            endif 
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
    
        
        end subroutine BF_coupled_T6_test_eval
    
    
    
    
    
    
        subroutine BF_deviator_test_eval()
            ! Calculates a strain deviator matrix and tests 
            use global
            use set_precision
            use dimensionless
            use math_constants
            implicit none 
        
        
            real(kind=kreal) :: DOF_global(nnode*4), vector(nnode*4), coords(3)
            real(kind=kreal) :: total_sum, analytical, x1, x2, y1, y2, z1, &
                                z2, diff, X(4), Y(4), Z(4)
            integer :: i, igll, num, ielem, ind , dim_u
        
            DOF_global = 0.0
    
            dim_u = 3*nnode
            
            do ielem=1,nelmt
                do igll = 1, ngll
                    num    = g_num(igll, ielem)
                    coords = g_coord(:, num)
    
                    ind = (num-1)*3 
        
                    ! U dot 
                    DOF_global(ind+1) =  coords(1)*coords(2)
                    DOF_global(ind+2) =  2*coords(2)
                    DOF_global(ind+3) =  -coords(3)
                enddo 
            enddo 
        
      
        
            ! We now have the global Dmat and the U (DOF) values: 
            vector(:) = matmul(Dmat, DOF_global)
            total_sum = 0 
            do i = 1, nnode*4
                total_sum = total_sum + vector(i)
            enddo 
        
        
            ! Now calculating the analytical solution for Test 1: 
            ! Coordinates are init. in global and calculated in driver code
            x1 = model_minx*NONDIM_L;  x2 = model_maxx*NONDIM_L
            y1 = model_miny*NONDIM_L;  y2 = model_maxy*NONDIM_L
            z1 = model_minz*NONDIM_L;  z2 = model_maxz*NONDIM_L
            
            ! Calc. squares and cubes of each term 
            do i=1,4
                X(i)  = (x2**i) - (x1**i)
                Y(i)  = (y2**i) - (y1**i)
                Z(i)  = (z2**i) - (z1**i)
            enddo 
    
    
            
            analytical = TWO_THIRD*( FIVE_HALF   * X(2)*Y(1)*Z(1)   +      &
                                     FIVE_HALF   * X(1)*Y(2)*Z(1)   +      &
                                     FIVE_HALF   * X(1)*Y(1)*Z(2)   -      &
                                     ONE_QUARTER * X(2)*Y(2)*Z(1)   -      &
                                     ONE_THIRD   * X(1)*Y(3)*Z(1)   -      & 
                                     ONE_QUARTER * X(1)*Y(2)*Z(2)          & 
                                   )
    
        
            diff = total_sum - analytical 
        
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*)'precision           :  ', BA_TEST_PRECISION
            write(*,*)'Analytical integral :  ', analytical
            write(*,*)'SEM Integral        : ', total_sum
            write(*,*)'difference          : ', diff
            write(*,*)'percentage error    : ', diff*100/analytical,' %'
            if (diff.gt.BA_TEST_PRECISION) then 
                write(*,*)'!!!!! Test failed: BF_STRAIN_DEVIATOR !!!!!!'
                stop
            else
                write(*,*)'PASSED: BF_STRAIN_DEVIATOR'
            endif 
            write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
    
        
        end subroutine BF_deviator_test_eval
    
    
    
    
        end module bilinear_form