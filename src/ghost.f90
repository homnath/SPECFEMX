! Ghost.f90 
! last edited: 3 Jun 2022 by WE 

module ghost 
    contains 

    subroutine modify_ghost_gdof(num, egdof, egdofu, coord, deriv, jac, bmat, &
        eld, eload, bload, vload, nodalu, nodalphi, nodalg, nodalB )
        ! Uses: 
        use global ! uses nproc, nenode
        use ghost_library_mpi ! uses prepare_ghost_gdof
        
        implicit none 
        ! I/O variables: 
        integer          , allocatable :: num(:), egdof(:), egdofu(:)

        real(kind=kreal) , allocatable :: coord(:,:), deriv(:,:), jac(:,:), bmat(:,:), &
                                          eld(:), eload(:),bload(:), vload(:), nodalu(:,:), &
                                          nodalphi(:), nodalg(:,:), nodalB(:,:)

        ! Local variables
        integer :: istat 


        !-------------------------------------
        ! modify ghost GDOFs
        if(nproc.gt.1)then
            call prepare_ghost_gdof()
        endif
        allocate(num(nenode),coord(ngnode,ndim),jac(ndim,ndim),deriv(ndim,nenode),     &
        bmat(nst,nedofu),eld(nedofu),bload(nedofu),vload(nedofu),eload(nedofu),        &
        nodalu(nndofu,nnode),egdof(nedof),egdofu(nedofu),stat=istat)
        if (istat/=0)then
            write(logunit,*)'ERROR: cannot allocate memory!'
            flush(logunit)
            stop
        endif
        if(ISPOT_DOF)then
            allocate(nodalphi(nnode),nodalg(ndim,nnode),nodalB(ndim,nnode),stat=istat)
            if(istat/=0)then
            write(logunit,*)'ERROR: cannot allocate memory!'
            flush(logunit)
            stop
            endif
        endif
        !--------------------------------
        
        return 
        
    end subroutine modify_ghost_gdof

end module ghost 

