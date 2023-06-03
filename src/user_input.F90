module user_input 

contains
!_______________________________________________________________________________

subroutine process_user_input(cmd, tdate, ttime, tzone, path, &
                              ext, format_str, errcode, errtag, &
                              ismesh_only,arg1,arg2,inp_fname,prog, &
                              cpu_tstart)

! USES
use global
use package_version
use shared
use string_library
use input 
#if(USE_MPI)
use mpi_library
#else 
use serial_library
#endif 
implicit none 
! IO variables 
logical :: ismesh_only
character(len=250) :: errtag ! error message
integer :: errcode
character(len=250) :: arg1,arg2,inp_fname,prog
character(len=150) :: path
character(len=20) :: ext,format_str
character(len=250) :: cmd ! command line
character(len=8) :: tdate ! date
character(len=10) :: ttime ! time
character(len=5) :: tzone ! time zone
real(kind=kreal) :: cpu_tstart
! Local variables 
integer ::ios

! Code 

call get_command_argument(0, prog)
if (command_argument_count() <= 0) then
  errcode=-1
  errtag='ERROR: no input file!'
  call control_error(errcode,errtag,stdout,myrank)
endif

call get_command_argument(1, arg1)
if(trim(arg1)==('--help'))then
  if(myrank==0)then
    write(stdout,'(a)')'Usage: '//trim(prog)//' [Options] [input_file]'
    write(stdout,'(a)')'Options:'
    write(stdout,'(a)')'    --help        : Display this information.'
    write(stdout,'(a)')'    --version     : Display version information.'
  endif
  !call sync_process
  call close_process()
elseif(trim(arg1)==('--version'))then
  if(myrank==0)then
    write(stdout,'(a)')trim(packname)//' '//trim(packver)//' '//trim(packtype)
    write(stdout,'(a)')'This is free software; see the source for copying '
    write(stdout,'(a)')'conditions.  There is NO warranty; not even for '
    write(stdout,'(a)')'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
  endif
  !call sync_process
  call close_process()
endif
ismesh_only=.false.
call get_command_argument(2, arg2)
if(trim(arg2)==('--mesh-only'))then
  ismesh_only=.true.
endif

! get input file name
call get_command_argument(1, inp_fname)

! parse input file name for future use
call parse_file(inp_fname,path,file_head,ext)

! open log file
! log file must be opened after parsing the input file name since it uses
! file_head!
! log file must be stored in the current directory since out_path isn't known
! at this stage!
if(myrank==0)then
  log_file = trim(file_head)//'.log'
  open(unit=logunit,file=trim(log_file),status='replace',action='write',iostat=ios)
  if(ios.ne.0)then
    write(errtag,'(a)')'ERROR: cannot open log file: '//trim(log_file)
    call control_error(errcode,errtag,stdout,myrank)
  endif
  write(logunit,'(a)')'--------------------------------------------'
  write(logunit,'(a)')'Result summary produced by '//&
  trim(packname)//' '//trim(packver)//' '//trim(packtype)
  write(logunit,'(a)')'--------------------------------------------'
  write(logunit,'(a)')'DATE(CCYYMMDD) TIME(HHMMSS.SSS) ZONE(+-HHMM)'
  call date_and_time(tdate,ttime,tzone)
  write(logunit,'(3x,a,7x,a,7x,a)')tdate,ttime,tzone
  call get_command(cmd)
  write(logunit,'(a)')trim(cmd)
  flush(logunit)
endif

! starting timer
call cpu_time(cpu_tstart)

! get processor tag
ptail=proc_tag()
if(ismpi.and.nproc.gt.1)then
  ptail_inp=trim(ptail)
else
  ptail_inp=''
endif

proc_str=''
if(ismpi.and.nproc.gt.1)then
  write(format_str,*)ceiling(log10(real(nproc)+1.))
  format_str='(i'//trim(adjustl(format_str))//')'
  write(proc_str,fmt=format_str)myrank
endif

! read input data
call read_input(inp_fname,errcode,errtag)
call sync_process()
call control_error(errcode,errtag,stdout,myrank)

! check method
if (trim(method)/='sem')then
  write(errtag,'(a)')'ERROR: wrong input for SPECFEM3D!'
  call control_error(errcode,errtag,stdout,myrank)
else
  !write(logunit, '(a)')'Correct method: sem'
endif

end subroutine process_user_input
!-------------------------------------------------------------------------------

end module user_input
!===============================================================================
