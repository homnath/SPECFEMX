module bilinear_form
    implicit none 

    contains 



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






    subroutine BA_poissons_T1()
        ! Calculates the coefficient matrix for bilinear form Term 1
        ! due to poissons 

        ! USES 
        use set_precision
        use global, only: nelmt, ngll, g_num, ngnode, NDIM, nenode, g_coord, &
                          weJACINV, weDETJAC, weJAC
        use integration, only: dshape_hex8, dlagrange_gll, gll_weights
        use element, only: hex8_gnode
        use math_library,only:determinant,invert,issymmetric

        implicit none 

        ! Local variables: 
        real(kind=kreal) :: Pois(nelmt,ngll,ngll)
        real(kind=kreal) :: jacw_bars, i_sum, &
                            t1, t2, quad_sum
               
        integer :: i_elem, abg, stv, bars, j,i,q 

        ! Initialise
        Pois = 0.0
        
        do i_elem = 1, nelmt ! loop elements
            do abg = 1, ngll               ! loop abg summation
                do stv = 1, ngll           ! loop stv summation
                    quad_sum = 0 
                    do bars = 1, ngll      ! loop quadrature summation
            
                        ! Calculate jacobian with respect to bars: 
                        jacw_bars   = gll_weights(bars) * weDETJAC(bars, i_elem)


                        i_sum = 0 
                        do i = 1,3

                            t1 = 0 
                            do j = 1,3
                                t1 = t1 + (weJACINV(j,i, stv, i_elem) * dlagrange_gll(j,bars,abg))
                            enddo 

                            t2 = 0 
                            do q = 1,3
                                t2 = t2 + (weJACINV(q,i, stv, i_elem) * dlagrange_gll(q,bars,stv))
                            enddo 

                            i_sum = i_sum + (t1*t2)
                        enddo
                        quad_sum = quad_sum + (i_sum * jacw_bars)
                    enddo   

                    Pois(i_elem, stv, abg) = quad_sum   
                enddo                       
            enddo  
        enddo 


        write(*,*)"POISSON:"
        write(*,*)Pois(1,:,:)


        
    end subroutine BA_poissons_T1

    end module bilinear_form