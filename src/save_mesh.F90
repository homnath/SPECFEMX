module save_mesh

contains 


subroutine write_original_mesh()
    ! USES
    use global 
    use output_to_user

    
    implicit none 
    ! IO variables
    ! Local variables
    ! Code: 


    ! write original meshes
    log_msg = trim('writing original mesh...') ;   call write_ifproc0()
    if(infbc)then
    ! classify finite/infinite elements for multiblock data plot
    npart=3
    allocate(ipart(npart),spart(npart))
    ! Note: we will NOT write all parts in a single file. Therefore, ipart must be
    ! 1 for all.
    ipart=(/ (1,i=1,npart) /)
    spart(1)='finite_domain'
    spart(2)='trinfinite_domain'
    spart(3)='infinite_domain'
    ! classify (in)finite elements
    !call classify_finite_infinite_elements(0)

    ! finite region
    ! write Ensight gold .case file
    case_file=trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.case'
    geo_file=trim(file_head)//'_original'//trim(ptail)//'.geo'
    add_tag='_original'
    call write_ensight_casefile(case_file,geo_file,add_tag,errcode,errtag)
    call control_error(errcode,errtag,stdout,myrank)
    ! write Ensight gold .geo file
    geo_file=trim(out_path)//trim(geo_file)
    call write_ensight_geo_part1(geo_file,ensight_hex8,ipart,spart,1, &
    nelmt_finite,nnode_finite,node_finite,nnode,real(g_coord),g_num_finite)

    ! trinfinite region
    ! write Ensight gold .case file
    trinfcase_file=trim(out_path)//trim(file_head)//'_original_trinf'//trim(ptail)//'.case'
    trinfgeo_file=trim(file_head)//'_original_trinf'//trim(ptail)//'.geo'
    add_tag='_original_trinf'
    call write_ensight_casefile(trinfcase_file,trinfgeo_file,add_tag,errcode,errtag)
    call control_error(errcode,errtag,stdout,myrank)
    ! write Ensight gold .geo file
    trinfgeo_file=trim(out_path)//trim(trinfgeo_file)
    call write_ensight_geo_part1(trinfgeo_file,ensight_hex8,ipart,spart,2, &
    nelmt_trinfinite,nnode_trinfinite,node_trinfinite,nnode,real(g_coord),g_num_trinfinite)

    ! infinite region
    ! write Ensight gold .case file
    infcase_file=trim(out_path)//trim(file_head)//'_original_inf'//trim(ptail)//'.case'
    infgeo_file=trim(file_head)//'_original_inf'//trim(ptail)//'.geo'
    add_tag='_original_inf'
    call write_ensight_casefile(infcase_file,infgeo_file,add_tag,errcode,errtag)
    call control_error(errcode,errtag,stdout,myrank)
    ! write Ensight gold .geo file
    infgeo_file=trim(out_path)//trim(infgeo_file)
    call write_ensight_geo_part1(infgeo_file,ensight_hex8,ipart,spart,3, &
    nelmt_infinite,nnode_infinite,node_infinite,nnode,real(g_coord),g_num_infinite)

    else
    npart=1
    allocate(ipart(npart),spart(npart))
    ipart=(/ (i,i=1,npart) /)
    spart(1)='finite_domain'
    ! write Ensight gold .case file
    case_file=trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.case'
    geo_file=trim(file_head)//'_original'//trim(ptail)//'.geo'
    add_tag='_original'
    call write_ensight_casefile(case_file,geo_file,add_tag,errcode,errtag)
    call control_error(errcode,errtag,stdout,myrank)

    ! write Ensight gold .geo file
    geo_file=trim(out_path)//trim(geo_file)
    call write_ensight_geo(geo_file,ensight_hex8,ipart,spart,nelmt,nnode,    &
    real(g_coord),g_num)
    endif


    ! write cell model for original mesh
    if(savedata%model_cell)then
    call write_model_cell(errcode,errtag)
    endif

    log_msg = trim('complete!') ;   call write_ifproc0()

end subroutine write_original_mesh



















subroutine save_mesh_ensight(infcase_file,infgeo_file,trinfcase_file, &
    trinfgeo_file,isgeo_change,add_tag, twidth, &
    fscase_file,fsgeo_file, fspcase_file,fspgeo_file, &
    ns,fi,fs,ts, errcode, errtag, ext,format_str, &
    case_file,geo_file, ipart, spart,spart_fs, buffer, node_hex8, gnum_hex8, &
    gnum_quad4,node_quad4)

! USES 
use global 
use write_ensight
use element
use free_surface
use postprocess
#if(USE_MPI)
use mpi_library
#else
use serial_library
#endif

implicit none 
! IO variables
integer :: twidth, ns,fi,fs,ts
integer,allocatable :: ipart(:)
character(len=80),allocatable :: spart(:) ! this must be 80 characters long
character(len=80) :: spart_fs(1) ! this must be 80 characters long
integer :: gnum_quad4(4),node_quad4(4)

character(len=80) :: buffer ! this must be 80 characters long
integer :: node_hex8(8), gnum_hex8(8)
character(len=250) :: infcase_file,infgeo_file,trinfcase_file,trinfgeo_file
character(len=250) :: fscase_file,fsgeo_file
character(len=250) :: fspcase_file,fspgeo_file
character(len=250) :: case_file,geo_file
logical :: isgeo_change
character(len=60) :: add_tag
character(len=20) :: ext,format_str

character(len=250) :: errtag ! error message
integer :: errcode
! Local variables
integer :: iounit,iounit_inf,iounit_fs,i,j,k, i_elmt

! Code: 
if(savedata%infinite)then
infcase_file=trim(out_path)//trim(file_head)//'_inf'//trim(ptail)//'.case'
if(nexcav==0)then
infgeo_file=trim(file_head)//'_inf'//trim(ptail)//'.geo'
else
isgeo_change=.true.
infgeo_file=trim(file_head)//'_inf'//'_step'//wild_char(1:twidth)//trim(ptail)//'.geo'
endif

add_tag='_inf'
call write_ensight_casefile_long(infcase_file,infgeo_file,add_tag,isgeo_change, &
ts,ns,fs,fi,twidth,errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)
endif

if(savedata%fsplot)then
fscase_file=trim(out_path)//trim(file_head)//'_free_surface'//trim(ptail)//'.case'
if(nexcav==0)then
fsgeo_file=trim(file_head)//'_free_surface'//trim(ptail)//'.geo'
else
isgeo_change=.true.
fsgeo_file=trim(file_head)//'_free_surface'//'_step'//wild_char(1:twidth)//trim(ptail)//'.geo'
endif

add_tag='_free_surface'
call write_ensight_casefile_long(fscase_file,fsgeo_file,add_tag,isgeo_change, &
ts,ns,fs,fi,twidth,errcode,errtag,freesurf=.true.,isplane=.false.)
call control_error(errcode,errtag,stdout,myrank)
endif

if(savedata%fsplot_plane)then
fspcase_file=trim(out_path)//trim(file_head)//'_free_surface_plane'//trim(ptail)//'.case'
if(nexcav==0)then
fspgeo_file=trim(file_head)//'_free_surface_plane'//trim(ptail)//'.geo'
else
isgeo_change=.true.
fspgeo_file=trim(file_head)//'_free_surface_plane'//'_step'//wild_char(1:twidth)//trim(ptail)//'.geo'
endif

add_tag='_free_surface_plane'
call write_ensight_casefile_long(fspcase_file,fspgeo_file,add_tag,isgeo_change, &
ts,ns,fs,fi,twidth,errcode,errtag,freesurf=.true.,isplane=.TRUE.)
call control_error(errcode,errtag,stdout,myrank)
endif




! Format string
write(tstep_sformat,*)twidth
tstep_sformat='i'//trim(adjustl(tstep_sformat))//'.'// &
trim(adjustl(tstep_sformat))
format_str='(a,'//trim(tstep_sformat)//',a)'
! write geo file for inital stage (original)
! open Ensight Gold geo file to store mesh data
if(isgeo_change)then
write(geo_file,fmt=format_str)trim(out_path)//trim(file_head)// &
'_step',0,trim(ptail)//'.geo'
else
write(geo_file,'(a)')trim(out_path)//trim(file_head)//trim(ptail)//'.geo'
endif
if(savedata%infinite)then
if(isgeo_change)then
write(infgeo_file,fmt=format_str)trim(out_path)//trim(file_head)//'_inf'// &
'_step',0,trim(ptail)//'.geo'
else
write(infgeo_file,'(a)')trim(out_path)//trim(file_head)//'_inf'//trim(ptail)//'.geo'
endif
endif

if(infbc)then
! write .geo file 
call write_ensight_geocoord_part1(geo_file,ipart,spart,1, &
nnode_finite,node_finite,nnode,real(g_coord),iounit)
if(savedata%infinite)then
call write_ensight_geocoord_part1(infgeo_file,ipart,spart,2, &
nnode_infinite,node_infinite,nnode,real(g_coord),iounit_inf)
endif
else
call write_ensight_geocoord(geo_file,ipart,spart,nnode,real(g_coord),iounit)
endif



! Dimensionalize coordinates after they are written.
! Writes element information.
buffer=ensight_hex8
if(infbc)then
write(iounit)buffer
write(iounit)nelmt_finite*(ngllx-1)*(nglly-1)*(ngllz-1)

! do not substract 1 for ensight file
do i_elmt=1,nelmt_finite
do k=1,ngllz-1
do j=1,nglly-1
do i=1,ngllx-1
! corner nodes in a sequential numbering
node_hex8(1)=(k-1)*ngllxy+(j-1)*ngllx+i
node_hex8(2)=node_hex8(1)+1

node_hex8(3)=node_hex8(1)+ngllx
node_hex8(4)=node_hex8(3)+1

node_hex8(5)=node_hex8(1)+ngllxy
node_hex8(6)=node_hex8(5)+1

node_hex8(7)=node_hex8(5)+ngllx
node_hex8(8)=node_hex8(7)+1
! map to exodus/cubit numbering and write
gnum_hex8=g_num_finite(node_hex8(map2exodus_hex8),i_elmt)
write(iounit)gnum_hex8
enddo
enddo
enddo
enddo
close(iounit)
! deallocate variables
deallocate(g_num_finite)!,g_num_infinite,node_finite,node_infinite)
! infinite region
if(savedata%infinite)then
write(iounit_inf)buffer
write(iounit_inf)nelmt_infinite*(ngllx-1)*(nglly-1)*(ngllz-1)

! do not substract 1 for ensight file
do i_elmt=1,nelmt_infinite
do k=1,ngllz-1
do j=1,nglly-1
do i=1,ngllx-1
! corner nodes in a sequential numbering
node_hex8(1)=(k-1)*ngllxy+(j-1)*ngllx+i
node_hex8(2)=node_hex8(1)+1

node_hex8(3)=node_hex8(1)+ngllx
node_hex8(4)=node_hex8(3)+1

node_hex8(5)=node_hex8(1)+ngllxy
node_hex8(6)=node_hex8(5)+1

node_hex8(7)=node_hex8(5)+ngllx
node_hex8(8)=node_hex8(7)+1
! map to exodus/cubit numbering and write
gnum_hex8=g_num_infinite(node_hex8(map2exodus_hex8),i_elmt)
write(iounit_inf)gnum_hex8
enddo
enddo
enddo
enddo
close(iounit_inf)
! deallocate variables
deallocate(g_num_infinite)!,g_num_infinite,node_finite,node_infinite)

endif
else
write(iounit)buffer
write(iounit)nelmt*(ngllx-1)*(nglly-1)*(ngllz-1)

! do not substract 1 for ensight file
do i_elmt=1,nelmt
do k=1,ngllz-1
do j=1,nglly-1
do i=1,ngllx-1
! corner nodes in a sequential numbering
node_hex8(1)=(k-1)*ngllxy+(j-1)*ngllx+i
node_hex8(2)=node_hex8(1)+1

node_hex8(3)=node_hex8(1)+ngllx
node_hex8(4)=node_hex8(3)+1

node_hex8(5)=node_hex8(1)+ngllxy
node_hex8(6)=node_hex8(5)+1

node_hex8(7)=node_hex8(5)+ngllx
node_hex8(8)=node_hex8(7)+1
! map to exodus/cubit numbering and write
gnum_hex8=g_num(node_hex8(map2exodus_hex8),i_elmt)
write(iounit)gnum_hex8
enddo
enddo
enddo
enddo
close(iounit)
endif

! Write GEO file for the free surface
if(savedata%fsplot)then
if(isgeo_change)then
write(fsgeo_file,fmt=format_str)trim(out_path)//trim(file_head)//'_free_surface'// &
'_step',0,trim(ptail)//'.geo'
else
write(fsgeo_file,'(a)')trim(out_path)//trim(file_head)//'_free_surface'//trim(ptail)//'.geo'
endif
spart_fs(1)='free_surface'
! write .geo file 
call write_ensight_geocoord_part1(fsgeo_file,ipart,spart_fs,1, &
nnode_fs,gnode_fs,nnode,real(g_coord),iounit_fs)



! Writes element information.
buffer=ensight_quad4
write(iounit_fs)buffer
! WARNING: statement/segment below assumes that ngllx=nglly=ngllz.
! It must be modified for unequal GLL points along different axes.
write(iounit_fs)nelmt_fs*(ngllx-1)*(nglly-1)



! Do not substract 1 for ensight file
do i_elmt=1,nelmt_fs
do j=1,nglly-1
do i=1,ngllx-1
! Corner nodes in a sequential numbering
node_quad4(1)=(j-1)*ngllx+i
node_quad4(2)=node_quad4(1)+1

node_quad4(3)=node_quad4(1)+ngllx
node_quad4(4)=node_quad4(3)+1

! Map to exodus/cubit numbering and write
gnum_quad4=rgnum_fs(node_quad4(map2exodus_quad4),i_elmt)
write(iounit_fs)gnum_quad4
enddo
enddo
enddo
close(iounit_fs)
endif

if(savedata%fsplot_plane)then
if(isgeo_change)then
write(fspgeo_file,fmt=format_str)trim(out_path)//trim(file_head)//'_free_surface_plane'// &
'_step',0,trim(ptail)//'.geo'
else
write(fspgeo_file,'(a)')trim(out_path)//trim(file_head)//'_free_surface_plane'//trim(ptail)//'.geo'
endif
spart_fs(1)='free_surface'
! write .geo file 
call write_ensight_geocoord_plane_part1(fspgeo_file,ipart,spart_fs,1,3, &
nnode_fs,gnode_fs,nnode,real(g_coord),iounit_fs)

! Writes element information.
buffer=ensight_quad4
write(iounit_fs)buffer
! WARNING: statement/segment below assumes that ngllx=nglly=ngllz.
! It must be modified for unequal GLL points along different axes.
write(iounit_fs)nelmt_fs*(ngllx-1)*(nglly-1)

! Do not substract 1 for ensight file
do i_elmt=1,nelmt_fs
do j=1,nglly-1
do i=1,ngllx-1
! Corner nodes in a sequential numbering
node_quad4(1)=(j-1)*ngllx+i
node_quad4(2)=node_quad4(1)+1

node_quad4(3)=node_quad4(1)+ngllx
node_quad4(4)=node_quad4(3)+1

! Map to exodus/cubit numbering and write
gnum_quad4=rgnum_fs(node_quad4(map2exodus_quad4),i_elmt)
write(iounit_fs)gnum_quad4
enddo
enddo
enddo
close(iounit_fs)


! Write Z-coordinate file
call write_scalar_to_file_freesurf(nnode_fs,g_coord(3,gnode_fs),ext='z', &
plane=.true.) 
endif


return


end subroutine save_mesh_ensight 



end module save_mesh