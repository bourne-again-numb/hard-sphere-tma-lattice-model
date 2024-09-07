program PT_SiO4_TAA
  use variables
  use input_SiO4_SDA
  use MPI
  implicit none
  
  integer :: opstatus,count,pid
  logical :: restart_file_exists
  logical,parameter :: restart_flag = .false.
  integer :: init_error,error
  integer :: comm_rank,comm_size,i
  character(len=200) :: proc_num,ifile,pwd
           
  !initialize MPI
  call MPI_INIT(init_error)
  !find out the rank of the processor
  call MPI_COMM_RANK(MPI_COMM_WORLD,comm_rank,error)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,comm_size,error)
  print*,'comm_size',comm_size
  
  print*,'---------------------------------'
  print*,'Rank of the processor:',comm_rank
  print*,'---------------------------------'
  
  !initialize random seed and start the system clock
  call system_clock(count)
  pid = getpid();pid=pid**2
  seed = - mod(1.0d0*count*pid,1.0d4)
  !seed = -3
  call cpu_time(start_time)
  
  write(proc_num,'(I2.2)')comm_rank
  call getcwd(pwd)
  ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('system_state.cfg')
  inquire(file=ifile,exist=restart_file_exists)
  
  !create output files
  call create_files_append(comm_rank)
  !read all the vaules for the variables
  call input
  
  !initialize the temperature grid
  call init_temp_grid(comm_rank)
  
  !setup the BCC lattice and its neighbor list
  call lattice
  !setup the counters for MC steps and snapshot count
  initial_step = 1
  snapshot = 0
  
  !Print simulation input conditions
  write(*,1000),'-----MC simulation starts-----'
  write(*,2000),'System size:',lx,ly,lz
  write(*,3000),'Concentration(SI/SN/SDA/W)//pH:',xi(1),xi(2),xi(3),xi(4),pH
  write(*,4000),'Molecules(SI/SN/SDA/W):',ni(1),ni(2),ni(3),ni(4)
  write(*,6000),'Temperature:',tstar
  write(*,10000),'SDA Properties/neighbors,radius(nm):',sda_size,sda_dia*a
  write(*,5000),'Ring penalties(3-rings/4rings):',pen3,pen4
  write(*,1000),'Interaction strengths:'
  write(*,6000),'O-SDA:',eOSDA
  !write(*,6000),'SN-SN:',eSNSN
  write(*,11000),'random seed',seed
  write(*,*)
  write(1008,1000)'MC simulation starts';flush(1008)
  write(1008,2000)'System size:',lx,ly,lz;flush(1008)
  write(1008,3000),'Concentration(SI/SN/SDA/W)//pH:',xi(1),xi(2),xi(3),xi(4),pH;flush(1008)
  write(1008,4000)'Molecules(SI/SN/SDA/W):',ni(1),ni(2),ni(3),ni(4);flush(1008)
  write(1008,6000)'Temperature:',tstar;flush(1008)
  write(1008,10000),'SDA Properties/neighbors,radius(nm):',sda_size,sda_dia*a;flush(1008)
  write(1008,5000)'Ring penalties(3-rings/4rings):',pen3,pen4;flush(1008)
  write(1008,1000)'Interaction strengths:';flush(1008)
  write(1008,6000),'O-SDA:',eOSDA;flush(1008)
  !write(1008,6000),'SN-SN:',eSNSN;flush(1008)
  write(1008,11000),'random seed',seed;flush(1008)
  write(1008,*);flush(1008)
    
  !INITIAL CONFIGURATION OF THE LATTICE
  if(initial_config_flag.eq.0)then
     call random_configuration(comm_rank)!random configuration
     !call read_state(comm_rank)
  elseif(initial_config_flag.eq.1)then
     if((xi(2)+xi(4)).eq.1.00)then
        call initial_cristo_config!beta-crystabolite configuration
     else
        print*,'BUG:: incorrect concen for crystobalite configurations'
        call freeandclose
        stop
     end if
  else
     print*,'Error: incorrect value for initial_configuration_flag :: main program'
     call freeandclose
     stop
  end if
  
  !start the PT MC algorithm
  call mc_algo(comm_rank,comm_size)
  
  call cpu_time(end_time)
  runtime = (end_time-start_time)/real(3600)
  !Print simulation summary
  ! write out the move summary
  write(*,1000)'===SN move summary==='
  write(*,*)'Jumps: success/attempts',SN_jumps_success,SN_jumps_attempt
  write(*,*)'Rotations: sucess/attempts:',SN_rotations_success,SN_rotations_attempt
  write(*,1000)'===SDA move summary==='
  !write(*,*)'Jumps: success/attempts',SDA_jumps_success,SDA_jumps_attempt
  write(*,*)'Translations: success/attempts',SDA_translations_success,SDA_translations_attempt
  write(*,1000)'===SN-SDA swap summary==='
  write(*,*)'SN-SDA swap: success/attempts',SN_SDA_swaps_success,SN_SDA_swaps_attempt
  ! write out the time summary
  write(*,8000),'Program runtime(hrs):',runtime
  write(1008,8000)'Program runtime(hrs):',runtime
  write(*,*);write(*,*)

  call MPI_BARRIER(MPI_COMM_WORLD,error)
  !output the final state of the system
  call output_state(comm_rank)
  !deallocate all the allocated memory
  call freeandclose
  
  call MPI_FINALIZE(error)
  
1000 format(x,A)
2000 format(x,A,i8,i8,i8)
3000 format(x,A,2x,f6.4,2x,f6.4,2x,f6.4,2x,f6.4,2x,f10.2)
4000 format(x,A,x,i8,x,i8,x,i8,x,i8)
5000 format(x,A,2x,f6.4,2x,f6.4)
6000 format(x,A,x,f7.4)
7000 format(x,A,i8,A,i8)
8000 format(A,es10.2)
9000 format(x,i10,2x,f6.4,2x,f6.4,2x,f6.4,2x,f6.4,2x,f7.4,2x,f7.4,2x,&
          f7.4,2x,f7.4,2x,f10.4,2x,f7.4)
10000 format(x,A,x,I8,5x,f6.4)
11000 format(x,A,x,i8)
contains
  !---------------------------------------------------------
  !--- INIT_TEMP_GRID: Initialized the temperature grid
  !---------------------------------------------------------
  subroutine init_temp_grid(comm_rank)
    use variables
    use MPI
    implicit none

    integer,intent(in) :: comm_rank
    
    !allocate the temperature matrix for every processor
    ntemp = comm_size  
    allocate(proc_temp(0:(ntemp-1)))

    !!setting up a linear scale in temperature
    do i = 0,ntemp-1
       proc_temp(i) = (ftemp-tstar)/real(ntemp-1)*real(i) + tstar
    end do

    !! additionally we can set up an exponential temperature grid where
    !! points at lower temperature are close to each other and at high temperature
    !! are farther apart T(i) = a*e^{b*i} 'a' and 'b' can be evaluated from
    !! T(0) = tstar T(ntemp-1) = ftemp
    !! To use this section, uncomment the next next few lines of the code
    !do i = 0,ntemp-1 !loop over points in temperature grid from 0->ntemp-1
    !   proc_temp(i) = tstar*(ftemp/tstar)**(real(i)/real(ntemp-1))
    !end do !enddo: loop over points in temperature grid from 1->ntemp-2

    !assigning every processor its temperature
    tstar = proc_temp(comm_rank)

    return
  end subroutine init_temp_grid
  !---------------------------------------------------------
  !--- CREATE_FILES: opens up the output files
  !---------------------------------------------------------
  subroutine create_files_append(comm_rank)
    implicit none
    
    integer,intent(in) :: comm_rank
    character(len=200) :: proc_num,proc_direc,ifile,pwd
    
    write(proc_num,'(I2.2)')comm_rank
    proc_direc = trim('temp_')//'_'//trim(proc_num)
    call getcwd(pwd)
    
    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('q.c_mcsteps.csv')
    open(1001,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'ERROR: opening Qn distribution file'
       call freeandclose
       stop
    end if

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('rings_dist.csv')
    open(1002,file=ifile,status='unknown',iostat=opstatus,access='append')    
    if(opstatus.ne.0)then
       print*,'ERROR: opening Qn distribution file'
       call freeandclose
       stop
    end if

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('lattice_energy.csv')    
    open(1003,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'ERROR: opening lattice_energy.csv'
       call freeandclose
       stop
    end if
    
    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('vmd.xyz')    
    open(1006,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'Error: opening vmd.xyz file'
       call freeandclose
       stop
    end if

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('cluster_stat.csv')        
    open(1007,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'Error: opening cluster_stat.csv file'
       call freeandclose
       stop
    end if

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('screen_output.log')    
    open(1008,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'Error: opening screen_output.log file'
       call freeandclose
       stop
    end if

    !ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('temperature_profile.csv')
    !open(1009,file=ifile,status='unknown',access='append',iostat=opstatus)
    !if(opstatus.ne.0)then
    !   print*,'Error: opening temperature_profile.csv file'
    !   call freeandclose
    !   stop
    !end if

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('cluster_order.csv')            
    open(1012,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'Error: opening cluster_order file'
       call freeandclose
       stop
    end if

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('summary.csv')                
    open(1014,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'Error: opening summary file'
       call freeandclose
       stop
    end if

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('qn.c_mcsteps.csv')                    
    open(1015,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'Error: opening qn.c_mcsteps file'
       call freeandclose
       stop
    end if
    
    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('particle_dimension.csv')                    
    open(1018,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'Error: opening size_range.csv file'
       call freeandclose
       stop
    end if    

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('pddf.csv')                    
    open(1019,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'Error: opening pddf.csv file'
       call freeandclose
       stop
    end if   

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('molecules_vmd.xyz')                    
    open(1021,file=ifile,status='unknown',iostat=opstatus,access='append')
    if(opstatus.ne.0)then
       print*,'Error: opening molecules_vmd.xyz file'
       call freeandclose
       stop
    end if

    !open(5000,file='exchange_info.csv',status='unknown',iostat=opstatus,access='sequential')
    open(5000,file='exchange_info.csv',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening exchange_info.csv file'
       call freeandclose
       stop
    end if

  end subroutine create_files_append
  !---------------------------------------------------------
  !--- LATTICE: CONSTRUSTS THE NEIGHBOR LIST
  !---------------------------------------------------------
  subroutine lattice
    use variables
    implicit none
    
    integer :: ix, iy, iz, itag, icheck, isite
    integer :: ix1,ix2,iy1,iy2,iz1,iz2
    integer :: ixs1,ixs2,iys1,iys2,izs1,izs2
    
    !initialize the coordinates array and the tag to site transformation array
    isite=0
    do ix=1,2*lx
       do iy=1,2*ly
          do iz=1,2*lz
             icheck=mod(ix,2)+mod(iy,2)+mod(iz,2)
             if(icheck.eq.0.or.icheck.eq.3) then
                isite = isite + 1
                itag = (ix-1)*2*ly*2*lz + (iy-1)*2*lz + iz
                rx(isite) = ix
                ry(isite) = iy
                rz(isite) = iz
                tag_site(itag) = isite
             end if
          end do
       end do
    end do
            
    !Potential Neighbor List Errors
    if(isite.ne.nsites) then
       print*,"Error:nlist(site.neq.sites)::subroutine lattice",isite,nsites
       call freeandclose
       stop
    end if
    
    !assign the first,second and third nearest neighbors
    isite=0
    do ix = 1,2*lx
       do iy = 1,2*ly
          do iz = 1,2*lz
             icheck = mod(ix,2) + mod(iy,2) + mod(iz,2)
             if(icheck.eq.0.or.icheck.eq.3) then
                isite = isite + 1
                itag = (ix-1)*2*ly*2*lz + (iy-1)*2*lz + iz

                ! first neighbor periodic boundary conditions
                ix1 = ix - 1
                if(ix1.lt.1) ix1 = 2*lx
                ix2 = ix + 1
                if(ix2.gt.2*lx) ix2 = 1
                iy1 = iy - 1
                if(iy1.lt.1) iy1 = 2*ly
                iy2 = iy + 1
                if(iy2.gt.2*ly) iy2 = 1
                iz1 = iz - 1
                if(iz1.lt.1) iz1 = 2*lz
                iz2 = iz + 1
                if(iz2.gt.2*lz) iz2 = 1

                ! first neighbors: within a dist. of 'a/(2)^0.5' frm a site
                n1list(1,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy1-1)*2*lz + iz2)
                n1list(2,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy1-1)*2*lz + iz2)
                n1list(3,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy2-1)*2*lz + iz2)
                n1list(4,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy2-1)*2*lz + iz2)
                n1list(5,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy1-1)*2*lz + iz1)
                n1list(6,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy1-1)*2*lz + iz1)
                n1list(7,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy2-1)*2*lz + iz1)
                n1list(8,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy2-1)*2*lz + iz1)
                
                ixs1 = ix-2
                ixs2 = ix+2
                iys1 = iy-2
                iys2 = iy+2
                izs1 = iz-2
                izs2 = iz+2
                !second and third neighbor periodic boundary conditions
                if(ixs1.lt.1) ixs1 = ixs1 + 2*lx
                if(ixs2.gt.2*lx) ixs2 = ixs2 - 2*lx
                if(iys1.lt.1) iys1 = iys1 + 2*ly
                if(iys2.gt.2*ly) iys2 = iys2 - 2*ly
                if(izs1.lt.1) izs1 = izs1 + 2*lz
                if(izs2.gt.2*lz) izs2 = izs2 - 2*lz
                !second neighbor: within a dist. of 'a' from a site
                n2list(1,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iy-1)*2*lz + iz)
                n2list(2,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iy-1)*2*lz + iz)
                n2list(3,isite) = tag_site((ix-1)*2*ly*2*lz + (iys1-1)*2*lz + iz)
                n2list(4,isite) = tag_site((ix-1)*2*ly*2*lz + (iys2-1)*2*lz + iz)
                n2list(5,isite) = tag_site((ix-1)*2*ly*2*lz + (iy-1)*2*lz + izs1)
                n2list(6,isite) = tag_site((ix-1)*2*ly*2*lz + (iy-1)*2*lz + izs2)
                
                !third neigbors within a distance of '2^0.5*a' from a site
                n3list(1,isite) = tag_site((ix-1)*2*ly*2*lz + (iys1-1)*2*lz + izs2)
                n3list(2,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys1-1)*2*lz + iz)
                n3list(3,isite) = tag_site((ix-1)*2*ly*2*lz + (iys1-1)*2*lz + izs1)
                n3list(4,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys1-1)*2*lz + iz)
                n3list(5,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iy-1)*2*lz + izs2)
                n3list(6,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iy-1)*2*lz + izs2)
                n3list(7,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iy-1)*2*lz + izs1)
                n3list(8,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iy-1)*2*lz + izs1)
                n3list(9,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys2-1)*2*lz + iz)
                n3list(10,isite) = tag_site((ix-1)*2*ly*2*lz + (iys2-1)*2*lz + izs2)
                n3list(11,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys2-1)*2*lz + iz)
                n3list(12,isite) = tag_site((ix-1)*2*ly*2*lz + (iys2-1)*2*lz + izs1)
                
                !fourth neighbors: within a distance of '3^0.5*a' from a site
                n4list(1,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys1-1)*2*lz + izs2)
                n4list(2,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys1-1)*2*lz + izs2)
                n4list(3,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys2-1)*2*lz + izs2)
                n4list(4,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys2-1)*2*lz + izs2)
                n4list(5,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys1-1)*2*lz + izs1)
                n4list(6,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys1-1)*2*lz + izs1)
                n4list(7,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys2-1)*2*lz + izs1)
                n4list(8,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys2-1)*2*lz + izs1)
            end if
          end do
       end do
    end do
    
    !Potential Neighbor List Errors
    if(isite.ne.nsites) then
       print*,"Error:nlist(site.neq.sites)::subroutine lattice",isite,nsites
       call freeandclose
       stop
    end if

    Si_O(1,1) = 1
    Si_O(2,1) = 3
    Si_O(3,1) = 6
    Si_O(4,1) = 8
    
    Si_O(1,2) = 2
    Si_O(2,2) = 4
    Si_O(3,2) = 5 
    Si_O(4,2) = 7
    
    Si_O(1,3) = 3
    Si_O(2,3) = 8
    Si_O(3,3) = 1
    Si_O(4,3) = 6

    Si_O(1,4) = 4
    Si_O(2,4) = 2
    Si_O(3,4) = 5
    Si_O(4,4) = 7
    
    Si_O(1,5) = 5
    Si_O(2,5) = 2
    Si_O(3,5) = 4
    Si_O(4,5) = 7
    
    Si_O(1,6) = 6
    Si_O(2,6) = 1
    Si_O(3,6) = 3
    Si_O(4,6) = 8
    
    Si_O(1,7) = 7
    Si_O(2,7) = 2
    Si_O(3,7) = 4
    Si_O(4,7) = 5
    
    Si_O(1,8) = 8
    Si_O(2,8) = 3
    Si_O(3,8) = 1
    Si_O(4,8) = 6
    
    Si_notO(1,1) = 2
    Si_notO(2,1) = 4
    Si_notO(3,1) = 5
    Si_notO(4,1) = 7
    
    Si_notO(1,2) = 3
    Si_notO(2,2) = 8
    Si_notO(3,2) = 1
    Si_notO(4,2) = 6
    
    Si_notO(1,3) = 4
    Si_notO(2,3) = 2
    Si_notO(3,3) = 5
    Si_notO(4,3) = 7
    
    Si_notO(1,4) = 1
    Si_notO(2,4) = 3
    Si_notO(3,4) = 6
    Si_notO(4,4) = 8
    
    Si_notO(1,5) = 6
    Si_notO(2,5) = 8
    Si_notO(3,5) = 1
    Si_notO(4,5) = 3
    
    Si_notO(1,6) = 7
    Si_notO(2,6) = 5
    Si_notO(3,6) = 2
    Si_notO(4,6) = 4
    
    Si_notO(1,7) = 8
    Si_notO(2,7) = 6
    Si_notO(3,7) = 1
    Si_notO(4,7) = 3
    
    Si_notO(1,8) = 2
    Si_notO(2,8) = 5
    Si_notO(3,8) = 4
    Si_notO(4,8) = 7
    
    call sda_size_allocation

    return
  end subroutine lattice
  !---------------------------------------------------------
  !--- SDA_SIZE_ALLOCATION: ALLOCATES SITES OCCUPIED BY SDA ON LATTICE
  !--------------------------------------------------------- 
  subroutine sda_size_allocation
    use variables
    use MPI
    implicit none

    integer :: isite, jsite, ksite, index
    integer :: i,j,k,jneigh,kneigh,prev_neighs
    real :: ijdistance,temp,jdist,kdist, half_body_diag
    integer,dimension(:),allocatable :: temp1
    real,dimension(:),allocatable :: temp_dist
    
    allocate(temp_dist(nsites))   
     allocate(temp1(nsites))    
    nlist = int(0);neighbors=int(0)

    !now we setup the neighbor list for every site on the lattice
    !this is done to accomodate an arbitrary size of the SDA molecule
    do isite = 1,nsites ! loop through all the sites
       temp_dist = real(0);temp1 = int(0)
       
       do jsite = 1,nsites ! loop through all the sub-sites of isite
          call distance(isite,jsite,ijdistance)
          temp_dist(jsite) = ijdistance
          temp1(jsite) = jsite
       end do ! loop through all the sub-sites of isite
       
       !sort the temp_dist array and temp1 array accordingly
       !call QsortC(temp_dist)
       do i = 1,nsites
          do j = 1,nsites-1
             if(temp_dist(j).gt.temp_dist(j+1))then
                temp=temp_dist(j)
                temp_dist(j) = temp_dist(j+1)
                temp_dist(j+1) = temp
                temp=temp1(j)
                temp1(j) = temp1(j+1)
                temp1(j+1) = temp
             end if
          end do
       end do
       
       !the ith entry in the neighbor array gives the total neighbors of the 
       !ith neighbor of a lattice site
       !Since, the neighbor array is common to all sites, so we can set 
       !it up only for the first site and use it everywhere else
       if(isite.eq.1)then! if isite.eq.1
          ! do jsite = 1,nsites
          !    write(100000,*)temp_dist(jsite),temp1(jsite)
          ! end do

          !now we setup the neighbor array          
          j = 2
          jneigh = 0;kneigh = 1
          do while(j.ne.nsites-1) ! loop over all the elements on the lattice
             jdist = temp_dist(j)
             kdist = temp_dist(j+1)
             
             if(jdist.eq.kdist)then
                kneigh = kneigh + 1
             else
                jneigh = jneigh + 1
                if(jneigh.gt.sda_size+1)exit
                neighbors(jneigh) = kneigh
                kneigh = 1
             end if
             
             j = j+1
          end do ! enddo: loop over all the elements on the lattice
          ! do jneigh = 1,sda_size+1
          !     write(200000,*)neighbors(jneigh)
          ! end do

          !safety check for SDA
          !the diameter of SDA should be less than half the body diagonal of the box
          !if(sda_dia.ge.())          
          prev_neighs = 1
          do jneigh = 1,sda_size
             prev_neighs = prev_neighs + neighbors(jneigh)
          end do
          sda_dia = temp_dist(prev_neighs + 1)
          half_body_diag = real(0.5)*sqrt(real( lx**2 + ly**2 + lz**2  ))
          if(sda_dia.ge.half_body_diag)then
             print*,'ERROR :: SDA is greater than hald the body diagonal'
             print*,'SDA diameter, half body diagonal',sda_dia,half_body_diag
             call freeandclose
             stop
          end if
          
       end if! endif: isite.eq.1
       
       !------------------
       !now after creating the neighbors array
       !we construct the nlist array for every isite
       do j = 1,sda_size + 1 ! loop over all the neighbors
          jneigh = neighbors(j)

          if(j.eq.1)then !for the first neighbor
             do k = 1,jneigh
                ksite = temp1(k+1)
                nlist(j,k,isite) = ksite
                !write(300000,*)temp_dist(k+1),nlist(j,k,isite)
             end do
          else !for second and higher neighbors

             !find the sum of all previous neighbors < j
             prev_neighs = 0
             do k = 1,j-1
                prev_neighs = prev_neighs + neighbors(k)
             end do
             do k = 1,jneigh !now loop through the entire temp1 array
                ksite = temp1(k + prev_neighs + 1)
                nlist(j,k,isite) = ksite
                !write(300000,*)temp_dist(k + prev_neighs + 1),nlist(j,k,isite)
             end do !enddo: now loop through the entire temp1 array
             
          end if !for the first neighbor
          
          
       end do !enddo: loop over all the neighbors
       !------------------
       
    end do ! loop through all the sites
    
    deallocate(temp_dist)
    deallocate(temp1)
    return
  end subroutine sda_size_allocation
  !---------------------------------------------------------
  !---distance: calculates the distance between two sites
  !---------------------------------------------------------   
  subroutine distance(jsite,ksite,jkdistance)
    use variables
    implicit none
    
    integer,intent(in) :: jsite, ksite
    integer :: deltax,deltay,deltaz
    real,intent(out) :: jkdistance
    
    deltax = rx(jsite) - rx(ksite)
    deltay = ry(jsite) - ry(ksite)
    deltaz = rz(jsite) - rz(ksite)

    deltax = deltax - 2*lx*nint(real(deltax)/real(2*lx))
    deltay = deltay - 2*ly*nint(real(deltay)/real(2*ly))
    deltaz = deltaz - 2*lz*nint(real(deltaz)/real(2*lz))
    
    jkdistance = sqrt(real(deltax**2 + deltay**2 + deltaz**2))
    
    return
  end subroutine distance
  !---------------------------------------------------------
  !--- LATTICE_INITIALIZATION: SETUP THE INITIAL CONDITION
  !--------------------------------------------------------- 
  subroutine random_configuration(comm_rank)
    use variables
    implicit none
    
    integer,intent(in) :: comm_rank
    integer :: isite, ivertex
    logical :: can_insert            
    integer :: molecule,bond,occfn         
    integer,dimension(4) :: fn,fnsite

    !set the occupancy and spin of sites to that od water: 
    !all sites are occupied by water
    !the sites of lattice initially have no pointer
    occupancy = occW
    spin = spinW
    head=0
    noccupy = 0
    resite = 0
    
    !-----FIRST: PLACE ALL THE SDA MOLECULES ON THE LATTICE
    molecule = 1
    do while(molecule.le.ni(3))
       !select a random site
       isite = int(1+ran3(seed)*real(nsites))
       
       !just for placing one molecule, make sure system size is 1,1,1
       !isite = tag_site((1-1)*2*ly*2*lz + (1-1)*2*lz + 1)
       !ivertex = 5
       
       !check if the SDA molecule can be inserted in isite
       call can_insert_SDA(isite,can_insert)
       
       if(can_insert)then!if the molecule can be inserted
          
          !change the occupancy of isite
          occupancy(isite) = occSDA
          !change the spin of isite
          spin(isite) = spinSDA
          !update the noccupy
          noccupy(ni(1)+ni(2)+molecule) = isite
          !update the resite arrays
          resite(isite) = molecule+ni(1)+ni(2)          
          
          molecule = molecule + 1
       end if!endif: the molecule can be inserted
    end do
    if(molecule-1.ne.ni(3))then
       print*,'BUG: molecule-1.ne.ni(3) :: subroutine initial_configuration'
       call freeandclose
       stop
    end if

    !-----SECOND: place all the SN on the sites
    molecule=1
    do while(molecule.le.ni(2))
       !select a random site among nsites
       isite = int(1 + ran3(seed)*real(nsites))
       !select a random vertex of the monomer
       ivertex = int(1 + ran3(seed)*real(2))
       
       !just for placing one molecule, make sure system size is 1,1,1
       !isite = tag_site((1-1)*2*ly*2*lz + (1-1)*2*lz + 1)
       !ivertex = 5

       !check for isite insertion and get the logical flag can_insert
       call can_insert_SN(isite,ivertex,can_insert)
       
       !if insertion of the monomer is possible
       if(can_insert) then

          !change the occupancy of the site to that of SN
          occupancy(isite) = occSN
          !update the occupancy of the vertices of the monomer
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          occupancy(fnsite(1)) = occupancy(fnsite(1)) + occSNO
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)
             occupancy(fnsite(bond)) = occupancy(fnsite(bond)) + occSNO
          end do

          !update the noccupy,ioxy,head and spin arrays
          noccupy(ni(1)+molecule) = isite
          head(isite) = ivertex
          spin(isite) = spinSN
          resite(isite) = ni(1) + molecule
          
          !proceed to inserting next molecule
          molecule = molecule+1          

       end if  !end checking can_insert if statement

    end do     !end loop over all the SN molecules
    if(molecule-1.ne.ni(2))then
       print*,'BUG: molecule-1.ne.ni(2) :: subroutine initial_configuration'
       call freeandclose
       stop
    end if
    !-----

    !-----THIRD:place all the SI molecules on the lattice
    molecule = 1
    do while(molecule.le.ni(1))
       ! select a random site
       isite = int(1 + ran3(seed)*real(nsites))
       !select a random orientation of the O-
       ivertex = int(1 + ran3(seed)*real(nc))

       !just for placing one molecule, make sure system size is 1,1,1
       !isite = tag_site((1-1)*2*ly*2*lz + (1-1)*2*lz + 1)
       !ivertex = 5
       
       !get the logical variable can_insert
       call can_insert_SI(isite,ivertex,can_insert)
       
       if(can_insert) then

          !change the occupancy of the site to that of SN
          occupancy(isite) = occSI
          !update the occupancy of the vertices of the monomer
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          occupancy(fnsite(1)) = occISIO
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)
             occupancy(fnsite(bond)) = occupancy(fnsite(bond)) + occSIO
          end do
          
          !update the noccupy,ioxy,head and spin  arrays
          noccupy(molecule) = isite
          head(isite) = ivertex
          spin(isite) = spinSI
          resite(isite) = molecule
          
          !proceed to inserting next molecule
          molecule = molecule + 1 
          
       end if!can_insert=yes if statement ends
       
    end do!end loop over all molecules of SI
        
    if(molecule-1.ne.ni(1))then
       print*,'BUG: molecule-1.ne.ni(1) :: subroutine initial_configuration'
       call freeandclose
       stop
    end if
    !-----
    
    !output the psf file
    call psf_out(comm_rank)
    
    !output initial configuration
    call config_out_atom
    
  end subroutine random_configuration
  !---------------------------------------------------------
  !--- INITIAL_CRITOBALITE CONFIGURATION
  !--- only apply this subroutine with SI/SDA = 0
  !--------------------------------------------------------- 
  subroutine initial_cristo_config
    use variables
    implicit none
    
    integer :: i,j,k,bond,opsttus,k1,k2
    integer,dimension(8) :: x,y,z,fn
    integer :: isite,m,ivertex,molecule
    integer,dimension(4) :: fn1,fnsite
    logical :: can_insert
    
    !setup the initial lattice state
    spin = spinW      !all water on site
    occupancy = occW  !all sites occupied by water
    head = 0          !no pointer for any site
    k1 = 1            !start with molecule 1
        
    !Only place neutral silica there
    do i = 1, lx / 4 
       !do i = 1, lx / 4 
       do j = 1, ly / 4
          do k = 1, lz / 4
         
             !unit cell
             !Point 1 (7, 7, 7)
             x(1) = i * 8 - 1
             y(1) = j * 8 - 1
             z(1) = k * 8 - 1
             fn(1) = 1
             
             !Point 2 (7, 3, 3)
             x(2) = i * 8 - 1
             y(2) = j * 8 - 5
             z(2) = k * 8 - 5
             fn(2) = 1
             
             !Point 3 (3, 3, 7)
             x(3) = i * 8 - 5
             y(3) = j * 8 - 5
             z(3) = k * 8 - 1
             fn(3) = 1
             
             !Point 4 (3, 7, 3)
             x(4) = i * 8 - 5
             y(4) = j * 8 - 1
             z(4) = k * 8 - 5
             fn(4) = 1
             
             !Point 5 (1, 1, 1)
             x(5) = i * 8 - 7
             y(5) = j * 8 - 7
             z(5) = k * 8 - 7
             fn(5) = 2
             
             !Point 6 (1, 5, 5)
             x(6) = i * 8 - 7
             y(6) = j * 8 - 3
             z(6) = k * 8 - 3
             fn(6) = 2
             
             !Point 7 (5, 1, 5)
             x(7) = i * 8 - 3
             y(7) = j * 8 - 7
             z(7) = k * 8 - 3
             fn(7) = 2
             
             !!	Point 8 (5, 5, 1)
             x(8) = i * 8 - 3
             y(8) = j * 8 - 3
             z(8) = k * 8 - 7
             fn(8) = 2
             
             do m = 1, 8
                
                isite = tag_site((x(m)-1)*(2*ly)*(2*lz) + (y(m)-1) * (2*lz)+z(m))
                noccupy(k1) = isite
                resite(isite) = k1
                
                call can_insert_SN(isite, fn(m), can_insert)
                if(can_insert)then
                   call insert_site(isite,fn(m),spinSN,k1)
                else
                   write(*,*) 'Initialize beta-cristobalite fail'
                   call freeandclose
                   stop 
                endif
                
                !Increase occupied site counters
            	
                k1 = k1 + 1
                
                if(k1.eq.(ni(2)+1)) then
                   !write the psf file
                   !call psf_out
                   !output initial configuration
                   call config_out_atom
                   return
                end if

             enddo
             !Coordinates for silica, on a bigger BCC lattice	
             
          enddo	!!	2*lz
          
       enddo	!!	2*ly
       
    enddo	!	2*lx
    
    return
  end subroutine initial_cristo_config
  !---------------------------------------------------------
  !--- PSF_OUT: GENERATES THE PSF FILE FOR VMD
  !--------------------------------------------------------- 
  subroutine psf_out(comm_rank)
    use variables
    implicit none
    
    integer,intent(in) :: comm_rank
    integer :: bond,i,j,k1,k2,sumoccupy,opstatus
    character(len=200) :: proc_num,proc_direc,ifile,pwd
        
    write(proc_num,'(I2.2)')comm_rank
    proc_direc = trim('temp_')//'_'//trim(proc_num)
    call getcwd(pwd)

    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('config.dat.psf')
    open(1005,file=ifile,status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'ERROR: opening config.dat.psf file'
       call freeandclose
       stop
    end if
        
    bond = ni(1)*4 + ni(2)*4
    sumoccupy = ni(1)*5 + ni(2)*5 + ni(3)
    
    write(1005,*) 'PSF   '
    write(1005,*) '      '
    write(1005,700) 1, '!NTITLE'
    write(1005,*) 'REMARKS original generated silica nanoparticle psf file'
    write(1005,*) '       '
    write(1005,700) sumoccupy,'!NATOM'
    
    !ionic silica
    j = 0
    do i=1,ni(1)
       j = j + 1
       write(1005,200) j
       j = j + 1
       write(1005,800) j
       do k1 = 1,3
          j = j + 1
          write(1005,200) j
       end do
    end do
    
    !neutral silica
    do i=1,ni(2)
       j = j + 1
       write(1005,100) j
       do k1 = 1,4
          j = j + 1
          write(1005,200) j
       end do
    end do
    
    !SDA
    do i=1,ni(3)
       j = j + 1
       write(1005,300) j
       !do k1=1,4
       !   j = j + 1
       !   write(1005,300)j
       !end do
    end do
    
    if(j.ne.natom)then
       print*,'BUG:j.ne.natom::subroutine psf_out'
       call freeandclose
       stop
    end if
    
    write(1005,*) '       '
    write(1005,700) bond, '!NBOND: bonds' 
    
    !for ionic silica
    k2 = 0
    do i=1,ni(1)
       write(1005,600) k2+1,k2+2,k2+1,k2+3,k2+1,k2+4,k2+1,k2+5
       k2 = k2 + 5
    end do
    
    !for neutral silica
    do i=1,ni(2)
       write(1005,600) k2+1,k2+2,k2+1,k2+3,k2+1,k2+4,k2+1,k2+5
       k2 = k2 + 5
    end do
    
    !for SDA
    do i=1,ni(3)
       write(1005,600) k2+1,k2+2,k2+1,k2+3,k2+1,k2+4,k2+1,k2+5
       k2 = k2 + 1
    end do
    
    if(k2.ne.natom)then
       print*,'BUG:k2.ne.natom::subroutine psf_out'
       call freeandclose
       stop
    end if
    
    write(1005,*) '       '
    write(1005,*) '       '
    write(1005,700) 0, '!NTHETA: angles'
    write(1005,*) '       '
    write(1005,700) 0, '!NPHI: dihedrals'
    write(1005,*) '       '
    write(1005,700) 0, '!NIMPHI: impropers'
    write(1005,*) '       '
    write(1005,700) 0, '!NCRTERM: cross-terms'
    
    close(1005)
    
100 format(I8,x,'CHAINS 1',x,'CHS',x,'CHS',x,'CHS',x,'0.000000 1.000000     0')
200 format(I8,x,'CHAINO 1',x,'CHO',x,'CHO',x,'CHO',x,'0.000000 1.000000     0')    
300 format(I8,x,'CHAINN 1',x,'CHN',x,'CHN',x,'CHN',x,'0.000000 1.000000     0')
400 format(I8,x,'CHAINC 1',x,'CHC',x,'CHC',x,'CHC',x,'0.000000 1.000000     0')
500 format(I8,x,'CHAINP 1',x,'CHP',x,'CHP',x,'CHP',x,'0.000000 1.000000     0')
800 format(I8,x,'CHAINH 1',x,'CHH',x,'CHH',x,'CHH',x,'0.000000 1.000000     0')
600 format(I8, I8, I8, I8, I8, I8, I8, I8)
700 format(I8,x,A) 
    
    return
  end subroutine psf_out
  !---------------------------------------------------------
  !--- CAN_INSERT_SN: CHECKS WHETHER SN CAN BE INSETED AT A SITE
  !---------------------------------------------------------   
  subroutine can_insert_SN(isite,ivertex,can_insert)
    use variables
    implicit none
    
    integer,intent(in) :: isite, ivertex
    logical,intent(out) :: can_insert
    integer :: occfn,neighs,overlap
    integer :: jsite,sn_neigh,jneigh,k
    integer :: bridge,i,j,bond,snbond         
    integer,dimension(4) :: fn,fnsite,sn,snsite
    
    !assume insertion is possible
    can_insert = .false.
    overlap = 0
    
    if(occupancy(isite).eq.occW)then!if isite has water on it
       
       !There's no need to uncomment the following lines,
       !since, these conditions would be taken care by the constraints
       !imposed on OH and their 1st,2nd and 3rd neighbors
       
       !check the sda_size neighbors of Si for any SDA
       do j = 1,sda_size
          jneigh = neighbors(j)
          do k = 1,jneigh
             jsite = nlist(j,k,isite)
             if(occupancy(jsite).eq.occSDA)then
                can_insert = .false.
                return
             end if
          end do
       end do
       
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       if((occupancy(fnsite(1)).eq.occSNO).or.((occupancy(fnsite(1)).eq.occSIO)))&
            overlap = overlap + 1
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          if((occupancy(fnsite(bond)).eq.occSNO).or.((occupancy(fnsite(bond)).eq.occSIO)))&
               overlap = overlap + 1
       end do
       
       !check if there are any water,SN,SNO or SIO molecules on the fnsite() sites
       do bond = 1,4
          occfn = occupancy(fnsite(bond))
          if((occfn.eq.occW).or.(occfn.eq.occSNO).or.(occfn.eq.occSIO))then
             can_insert = .true.
          else
             can_insert = .false.
             return
          end if
       end do
       
       do bond = 1,4!check sda_size neighbors of OH for any SDA occupancy
          
          !check the sda_size neighbors of Si for any SDA
          do j = 1,sda_size
             jneigh = neighbors(j)
             do k = 1,jneigh
                jsite = nlist(j,k,fnsite(bond))
                if(occupancy(jsite).eq.occSDA)then
                   can_insert = .false.
                   return
                end if
             end do
          end do

       end do!enddo: !check sda_size neighbors of OH for any SDA occupancy
             
       !-----CHECK FOR 2 MEMBERED RING FORMATION
       if(overlap.ge.2)then!if overlap>=2 chance of 2 membered ring

          do sn_neigh = 1,6!loop over 2 neighbors of isite
             
             jsite = n2list(sn_neigh,isite)
             !if there SN on neighboring site
             if((occupancy(jsite).eq.occSN).or.(occupancy(jsite).eq.occSI))then
                
                !find the bonded OH of jsite
                sn(1) = head(jsite)
                snsite(1) = n1list(sn(1),jsite)
                do snbond = 2,4
                   sn(snbond) = Si_O(snbond,sn(1))
                   snsite(snbond) = n1list(sn(snbond),jsite)
                end do
                
                !check for any bridging oxygens
                bridge=0
                do i=1,4
                   do j=1,4
                      if((fnsite(i).eq.snsite(j)).and.((occupancy(fnsite(i)).eq.occSNO).or.&
                           (occupancy(fnsite(i)).eq.occSIO)))then 
                         bridge = bridge + 1
                      end if
                   end do
                end do

                !if there are more than 1 bridging oxygens
                if(bridge.gt.1) then
                   can_insert = .false.
                   return
                end if
                
             end if!if there SN on neighboring site

          end do!enddo:loop over 2 neighbors of isite

       end if!if overlap.ge.2 chance of 2 membered ring
       !-----

    end if!if isite has water on it

    return
  end subroutine can_insert_SN  
  !----------------------------------------------------------
  !--- CAN_INSERT_SI: CHECKs WHETHER SI CAN BE INSERTED AT A SITE
  !---------------------------------------------------------- 
  subroutine can_insert_SI(isite,ivertex,can_insert)
    
    integer,intent(in) :: isite, ivertex
    logical,intent(out) :: can_insert
    integer :: jsite,jocc,neighs,overlap,ihead
    integer :: bridge,i,j,bond,snbond,sn_neigh
    integer,dimension(4) :: fn,fnsite,sn,snsite
    
    !assume insertion is possible
    can_insert = .false.
    overlap = 0
    
    !check if there's water on the site
    if(occupancy(isite).eq.occW)then!if isite has water on it
       
       !check first neighbors for SDA occupancy
       do i=1,8
          jocc = occupancy(n1list(i,isite))
          if(jocc.eq.occSDA)then
             can_insert = .false.
             return
          end if
       end do
       
      fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       ihead = fnsite(1)
       if(occupancy(fnsite(1)).eq.occW)then
       else
          can_insert = .false.
          return
       end if
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          if(occupancy(fnsite(bond)).eq.occSNO)&
               overlap = overlap + 1
       end do
       
       !check the first neighbors of OH for any SDA occupancy
       do bond = 1,4
          do j=1,8
             jocc = occupancy(n1list(j,fnsite(bond)))
             if(jocc.eq.occSDA)then
                can_insert = .false.
                return
             end if
          end do
       end do
       
       !search for the second and third neighbors of O- for any O-
       !if there are then don't insert SI
       !do j=1,6!loop over second neighbors
       !   jocc = occupancy(n2list(j,ihead))
       !   if(jocc.eq.occISIO)then
       !      can_insert = .false.
       !      return
       !   end if
       !end do!enddo:loop over second neighbors
       do j=1,12!loop over third neighbors
          jocc = occupancy(n3list(j,ihead))
          if(jocc.eq.occISIO)then
             can_insert = .false.
             return
          end if
       end do!enddo:loop over third neighbors
              
       !check if there are any water,SN,SI,SNO or SIO molecules on the 
       !fnsite() sites
       !check for the location of the pointer variable

       if(occupancy(fnsite(1)).eq.occW)then!if the pointer has water on its site
          do bond = 2,4
             jocc = occupancy(fnsite(bond))
             if((jocc.eq.occW).or.(jocc.eq.occSNO))then
                can_insert = .true.
             else
                can_insert = .false.
                return
             end if
          end do
          !-----CHECK FOR 2 MEMBERED RING FORMATION
          if(overlap.ge.2)then!if overlap.ge.2 chance of 2 membered ring
             
             do sn_neigh = 1,6!loop over 2 neighbors of isite
                
                jsite = n2list(sn_neigh,isite)
                !if there SN on neighboring site
                if(occupancy(jsite).eq.occSN.or.(occupancy(jsite).eq.occSI))then
                   
                   !find the bonded OH of jsite
                   sn(1) = head(jsite)
                   snsite(1) = n1list(sn(1),jsite)
                   do snbond = 2,4
                      sn(snbond) = Si_O(snbond,sn(1))
                      snsite(snbond) = n1list(sn(snbond),jsite)
                   end do
                   
                   !check for any bridging oxygens
                   bridge=0
                   do i=1,4
                      do j=1,4
                         if((fnsite(i).eq.snsite(j)).and.(occupancy(fnsite(i)).eq.occSNO))then 
                            bridge = bridge + 1
                         end if
                      end do
                   end do
                   
                   !if there are more than 1 bridging oxygens
                   if(bridge.gt.1) then
                      can_insert = .false.
                      return
                   end if
                   
                end if!if there SN on neighboring site
                
             end do!enddo:loop over 2 neighbors of isite
             
          end if!if overlap.ge.2 chance of 2 membered ring
          !-----
       else
          can_insert = .false.
          return
       end if!endif: the pointer has water on its site
       
    end if!endif: isite has water on it
    
    return
  end subroutine can_insert_SI
  !---------------------------------------------------------
  !--- CAN_INSERT_SDA(isite,can_insert): FIND THE Qn DISTRIBUTION
  !--------------------------------------------------------- 
  subroutine can_insert_SDA(isite,can_insert)
    use variables
    implicit none
    
    integer,intent(in) :: isite
    logical,intent(out) :: can_insert
    integer :: j,jsite,k,jneigh,jocc,occfn
    
    can_insert = .false.
    
    !check the occupancy of isite
    if(occupancy(isite).eq.occW)then!if there's water on isite
       !check the sda_size neighbors of Si for any SDA
       do j = 1,sda_size
          jneigh = neighbors(j)
          do k = 1,jneigh
             jsite = nlist(j,k,isite)
             if(occupancy(jsite).eq.occW)then
                can_insert = .true.
             else
                can_insert = .false.
                return
             end if
          end do
       end do
       
    else

       can_insert = .false.
       return

    end if!endif there's water on isite       
    
    return
  end subroutine can_insert_SDA
  !---------------------------------------------------------
  !--- Qn_DIST(Q0,Q1,Q2,Q3,Q4): FIND THE Qn DISTRIBUTION
  !--------------------------------------------------------- 
  subroutine Qn_dist
    use variables
    implicit none
    
    integer :: isite,ivertex,ispin,bond,occfn,centerocc
    integer :: overlap,monomer
    integer,dimension(4) :: fn,fnsite
    
    !set the accumulators to 0
    Qn = 0;Qni = 0;Qnn = 0
    
    do monomer=1,allmolecules!loop over all monomers
       
       isite = noccupy(monomer)
       ivertex = head(isite)
       ispin = spin(isite)
       centerocc = occupancy(isite)
       
       if((ispin.eq.spinSN).or.(ispin.eq.spinSI))then!if the molecule is SN/SI
          fn(1) = ivertex
          fn(2) = Si_O(2,ivertex)
          fn(3) = Si_O(3,ivertex)
          fn(4) = Si_O(4,ivertex)
          
          overlap = 0
          do bond=1,4
             fnsite(bond) = n1list(fn(bond),isite)
             occfn = occupancy(fnsite(bond))
             if((occfn.eq.occSNOSNO).or.(occfn.eq.occSIOSNO))then
                overlap = overlap + 1
             end if
          end do !end loop over hydroxide molecules
          
          !update the Qn variables
          if(overlap.eq.0) then
             Qn(0) = Qn(0) + 1
          elseif(overlap.eq.1) then
             Qn(1) = Qn(1) + 1
          elseif(overlap.eq.2) then
             Qn(2) = Qn(2) + 1
          elseif(overlap.eq.3) then
             Qn(3) = Qn(3) + 1
          elseif(overlap.eq.4) then
             Qn(4) = Qn(4) + 1
          else
             print*,'BUG: More than 4 bridging O to Si : subroutine Qn_distribution (Qn)'
             call freeandclose
             stop
          end if

       end if!endif the molecule is SN/SI
       
       !Qn distribution only for neutral silica
       if(ispin.eq.spinSN)then
          fn(1) = ivertex
          fn(2) = Si_O(2,ivertex)
          fn(3) = Si_O(3,ivertex)
          fn(4) = Si_O(4,ivertex)
          
          overlap = 0
          do bond=1,4
             fnsite(bond) = n1list(fn(bond),isite)
             occfn = occupancy(fnsite(bond))
             if((occfn.eq.occSNOSNO).or.(occfn.eq.occSIOSNO))then
                overlap = overlap + 1
             end if
          end do !end loop over hydroxide molecules
          
          !update the Qn variables
          if(overlap.eq.0) then
             Qnn(0) = Qnn(0) + 1
          elseif(overlap.eq.1) then
             Qnn(1) = Qnn(1) + 1
          elseif(overlap.eq.2) then
             Qnn(2) = Qnn(2) + 1
          elseif(overlap.eq.3) then
             Qnn(3) = Qnn(3) + 1
          elseif(overlap.eq.4) then
             Qnn(4) = Qnn(4) + 1
          else
             print*,'BUG: More than 4 bridging O to Si : subroutine Qn_distribution (Qnn)'
             call freeandclose
             stop
          end if
       end if

    end do!end loop over all the monomers
    
    return
  end subroutine Qn_dist
  !---------------------------------------------------------
  !---returnneigh: returns neighboring Si to a site
  !---------------------------------------------------------
  subroutine returnneigh(isite,isites,arethereSi)
    use variables
    implicit none
    
    integer,intent(in) :: isite
    !integer,dimension(:),allocatable,intent(out) :: neighSi
    logical,intent(out) :: arethereSi
    integer :: ispin,ivertex,bond,fnbond,OH,Obond,fnn,i,occfn
    integer :: neighboringSi,index
    integer,dimension(4) :: fn, fnsite, sn,snsite
    integer,dimension(4) :: isites
    logical :: connect
    
    ispin = spin(isite)
    ivertex = head(isite)
    
    !if(allocated(neighSi)) deallocate(neighSi)
    
    !set the accumulator to zero
    neighboringSi = 0
    isites=-1
    arethereSi = .false.
    
    if((ispin.eq.spinSI).or.(ispin.eq.spinSN))then
       !find the neighbors of isite
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       do bond = 2,4
          fn(bond) = Si_O(bond,ivertex)
          fnsite(bond) = n1list(fn(bond),isite)
       end do
    else
       print*,'BUG : molecule other than SI/SN :: subroutine returnneigh'
       print*,'identity:',ispin
       call freeandclose
       stop
    end if
    
    do bond = 1,4!find the db occ. O on nearest sites to isite

       OH = fnsite(bond)
       occfn = occupancy(OH)

       if((occfn.eq.occSNOSNO).or.(occfn.eq.occSIOSNO))then!if there are double occupies Os
          
          do Obond = 1,8!loop over the first neighbors of OH

             fnn = n1list(Obond,OH)!fnn:nearest neighboring Si to isite

             if(((occupancy(fnn).eq.occSN).or.(occupancy(fnn).eq.occSI)).and.&
                  (fnn.ne.isite))then!if there's Si neighboring OH

                !find the neighbors of fnn
                sn(1) = head(fnn)
                snsite(1) = n1list(sn(1),fnn)
                do fnbond = 2,4
                   sn(fnbond) = Si_O(fnbond,sn(1))
                   snsite(fnbond) = n1list(sn(fnbond),fnn)
                end do

                !check for a connection
                connect = .false.
                do fnbond = 1,4
                   if(snsite(fnbond).eq.OH)then
                      connect = .true.
                   end if
                end do

                if(connect)then
                   neighboringSi = neighboringSi + 1
                   isites(neighboringSi) = fnn
                end if

             end if!endif: there's SN/SI neighboring OH
             
          end do!loop over the first neighbors of OH

       end if!endif there are double occupies Os

    end do!enddo:find the db occ. O on nearest sites to isite
    
    index = 0
    if(neighboringSi.ge.1)then!isite is a part of a ring system
       arethereSi = .true.
       !allocate(neighSi(neighboringSi))
       !do i = 1,size(isites)
       !   if(isites(i).ne.0)then
       !      index = index + 1
       !      neighSi(index) = isites(i)
       !   end if
       !end do
    end if!endif:isite is a part of a ring system
    
    return
  end subroutine returnneigh
  !---------------------------------------------------------
  !---rings3and4: 3 and 4 membered ring calculation
  !--------------------------------------------------------- 
  subroutine rings3and4(rings3,rings4)
    use variables
    implicit none

    integer,intent(out) :: rings3,rings4
    integer :: rings3_temp,rings4_temp,molecule
    integer :: isite,ivertex,iocc,ispin,bond1,bond2
    integer :: jsite,jvertex,jocc
    integer :: ksite,kvertex,kocc
    integer :: lsite,lvertex,locc
    integer :: OH11,OH12,occOH11,occOH12,OH21,occOH21
    integer :: bondOH11,bondOH12,bondj,bondk,bondl,p,q,bondOH21
    integer,dimension(4) :: ifn,jfn,kfn,lfn,ifnsite,jfnsite,kfnsite,lfnsite
    logical :: connectedij,connectedik,connectedjl,connectedkl
    
    !set accumulators to zero
    rings3_temp = 0
    rings4_temp = 0
    rings3 = 0
    rings4 = 0
    
    ! if(pen3.eq.real(0))then
    !    rings3 = 0
    !    return
    ! end if
    ! if(pen4.eq.real(0))then
    !    rings4 = 0
    !    return
    ! end if
    
    do molecule = 1,allmolecules!loop over all molecules
       
       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if((ispin.eq.spinSN).or.(ispin.eq.spinSI))then!if the molecule is SN/SI
          
          do bond1 = 1,4
             ifn(bond1) = Si_O(bond1,ivertex)
             ifnsite(bond1) = n1list(ifn(bond1),isite)              
          end do
          
          do bond1 = 1,4!loop over the OH of isite
             OH11 = ifnsite(bond1);occOH11 = occupancy(OH11)
             do bond2 = 1,4!loop over the OH of isite
                OH12 = ifnsite(bond2);occOH12 = occupancy(OH12)
                
                !if there are more then 2 brigdes
                if(((occOH11.eq.occSNOSNO).or.(occOH11.eq.occSIOSNO)).and.&
                     ((occOH12.eq.occSNOSNO).or.(occOH12.eq.occSIOSNO)).and.(OH11.ne.OH12))then
                   
                   !find if there are any SN/SI connected to OH11
                   do bondOH11 = 1,8!loop over first neighbors of OH11
                      jsite = n1list(bondOH11,OH11);jocc = occupancy(jsite);jvertex = head(jsite)
                      
                      if((jvertex.gt.0).and.(jsite.ne.isite))then!if jsite has SN/SI
                         
                         !find the first neighbors of jsite
                         do bondj = 1,4
                            jfn(bondj) = Si_O(bondj,jvertex)
                            jfnsite(bondj) = n1list(jfn(bondj),jsite)
                         end do
                         
                         !check if there's a bridge between jsite and isite
                         connectedij = .false.
                         do bondj = 1,4
                            if(OH11.eq.jfnsite(bondj)) connectedij = .true.
                         end do
                         
                         if(connectedij)then!if there's a bridge between isite and jsite
                            
                            !find the neighboring SI/SN to OH12
                            do bondOH12 = 1,8!loop over first neighbors of OH12
                               ksite = n1list(bondOH12,OH12);kocc = occupancy(ksite);kvertex = head(ksite)
                               
                               if((kvertex.gt.0).and.(ksite.ne.isite))then!if ksite has SN/SI
                                  
                                  !find the first neighbors of ksite
                                  do bondk = 1,4
                                     kfn(bondk) = Si_O(bondk,kvertex)
                                     kfnsite(bondk) = n1list(kfn(bondk),ksite)
                                  end do
                                  
                                  !check if there's a brigde between ksite and isite
                                  connectedik = .false.
                                  do bondk = 1,4
                                     if(kfnsite(bondk).eq.OH12)connectedik = .true.
                                  end do
                                  
                                  if(connectedik)then!if there's a bridge between isite and ksite
                                     
                                     !---CALCULATE THE 3 MEMBERED RINGS---
                                     !older version
                                     !check if there are any common OH to ksite and jsite
                                     !do p=1,4
                                     !   do q=1,4
                                     !      if(kfnsite(q).eq.jfnsite(p))rings3_temp = rings3_temp + 1
                                     !   end do
                                     !end do
                                     !------------------------------------
                                     
                                     !---CALCULATE THE 4 MEMBERED RINGS
                                     do bondj = 1,4!loop over first neighbors of jsite
                                        OH21 = jfnsite(bondj);occOH21 = occupancy(OH21)
                                        
                                        if((OH21.ne.OH11).and.((occOH21.eq.occSNOSNO).or.&
                                             (occOH21.eq.occSIOSNO)))then!if there's another bridge of jsite
                                           
                                           do bondOH21 = 1,8!loop over first neighbors of OH21
                                              lsite = n1list(bondOH21,OH21);locc = occupancy(lsite);lvertex = head(lsite)
                                              
                                              if((lvertex.gt.0).and.(lsite.ne.jsite).and.(lsite.ne.ksite))then!if there's SN/SI on lsite
                                                 
                                                 !find the OH of lsite
                                                 do bondl = 1,4
                                                    lfn(bondl) = Si_O(bondl,lvertex)
                                                    lfnsite(bondl) = n1list(lfn(bondl),lsite)
                                                 end do
                                                 
                                                 !loop through the first neighbors of OH21 to find the bridge between lsite and jsite
                                                 connectedjl = .false.
                                                 do bondl = 1,4
                                                    if(lfnsite(bondl).eq.OH21)connectedjl = .true.
                                                 end do
                                                 
                                                 if(connectedjl)then
                                                    !check if there's any common OH between lsite and ksite
                                                    do p = 1,4
                                                       do q=1,4
                                                          if((kfnsite(p).eq.lfnsite(q)).and.(kfnsite(p).ne.OH12))&
                                                               rings4_temp = rings4_temp + 1
                                                       end do
                                                    end do
                                                 end if
                                                 
                                              elseif((lvertex.gt.0).and.(lsite.ne.jsite).and.(lsite.eq.ksite))then!if there's SN/SI on lsite
                                                 
                                                 !---CALCULATE 3 MEMBERED RINGS---
                                                 !newer version
                                                 !loop over OH of lsite and check if it is bonded to OH12
                                                 do bondk = 1,4
                                                    if(kfnsite(bondk).eq.OH12)rings3_temp = rings3_temp + 1
                                                 end do
                                                 !--------------------------------
                                                 
                                              end if!endif: there's SN/SI on lsite
                                              
                                           end do!enddo: loop over first neighbors of OH21
                                           
                                        end if!endif: there's another bridge of jsite site
                                        
                                     end do!enddo: loop over first neighbors of jsite
                                     !---------------------------------
                                     
                                  end if!endif: there's a bridge between isite and ksite
                                  
                               end if!endif: ksite has SN/SI
                               
                            end do!enddo: loop over first neighbors of OH12
                            
                         end if!endif: there's a bridge between isite and jsite
                         
                      end if!endif jsite has SN/SI
                      
                   end do!enddo: loop over first neighbors of OH11
                   
                end if!endif: there are more then 2 brigdes
                
             end do!enddo: loop over the OH of isite
          end do!enndo: loop over the OH of isite
          
       end if!endif: the molecule is SN/SI
       
    end do!enddo: loop over all molecules
    
    rings3 = rings3_temp/6
    rings4 = rings4_temp/8
    
    return
  end subroutine rings3and4
  !---------------------------------------------------------
  !---site_rings3and4
  !--------------------------------------------------------- 
  subroutine site_rings3and4(isite,ivertex,rings3,rings4)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex
    integer,intent(out) :: rings3,rings4
    integer :: rings3_temp,rings4_temp
    integer :: bond1,bond2,iocc
    integer :: jsite,jvertex,jocc
    integer :: ksite,kvertex,kocc
    integer :: lsite,lvertex,locc
    integer :: OH11,OH12,occOH11,occOH12,OH21,occOH21
    integer :: bondOH11,bondOH12,bondj,bondk,bondl,p,q,bondOH21
    integer,dimension(4) :: ifn,jfn,kfn,lfn,ifnsite,jfnsite,kfnsite,lfnsite
    logical :: connectedij,connectedik,connectedjl,connectedkl
    
    !set the accumulators to zero
    rings3_temp = 0
    rings4_temp = 0
    rings3 = 0
    rings4 = 0

    ! if(pen3.eq.real(0))then
    !    rings3 = 0
    !    return
    ! end if
    ! if(pen4.eq.real(0))then
    !    rings4 = 0
    !    return
    ! end if
    
    do bond1 = 1,4
       ifn(bond1) = Si_O(bond1,ivertex)
       ifnsite(bond1) = n1list(ifn(bond1),isite)              
    end do
    
    do bond1 = 1,4!loop over the OH of isite
       OH11 = ifnsite(bond1);occOH11 = occupancy(OH11)
       do bond2 = 1,4!loop over the OH of isite
          OH12 = ifnsite(bond2);occOH12 = occupancy(OH12)
          
          !if there are more then 2 brigdes
          if(((occOH11.eq.occSNOSNO).or.(occOH11.eq.occSIOSNO)).and.&
               ((occOH12.eq.occSNOSNO).or.(occOH12.eq.occSIOSNO)).and.(OH11.ne.OH12))then
             
             !find if there are any SN/SI connected to OH11
             do bondOH11 = 1,8!loop over first neighbors of OH11
                jsite = n1list(bondOH11,OH11);jocc = occupancy(jsite);jvertex = head(jsite)
                
                if((jvertex.gt.0).and.(jsite.ne.isite))then!if jsite has SN/SI
                   
                   !find the first neighbors of jsite
                   do bondj = 1,4
                      jfn(bondj) = Si_O(bondj,jvertex)
                      jfnsite(bondj) = n1list(jfn(bondj),jsite)
                   end do
                   
                   !check if there's a bridge between jsite and isite
                   connectedij = .false.
                   do bondj = 1,4
                      if(OH11.eq.jfnsite(bondj)) connectedij = .true.
                   end do
                   
                   if(connectedij)then!if there's a bridge between isite and jsite
                      
                      !find the neighboring SI/SN to OH12
                      do bondOH12 = 1,8!loop over first neighbors of OH12
                         ksite = n1list(bondOH12,OH12);kocc = occupancy(ksite);kvertex = head(ksite)
                         
                         if((kvertex.gt.0).and.(ksite.ne.isite))then!if ksite has SN/SI
                            
                            !find the first neighbors of ksite
                            do bondk = 1,4
                               kfn(bondk) = Si_O(bondk,kvertex)
                               kfnsite(bondk) = n1list(kfn(bondk),ksite)
                            end do
                            
                            !check if there's a brigde between ksite and isite
                            connectedik = .false.
                            do bondk = 1,4
                               if(kfnsite(bondk).eq.OH12)connectedik = .true.
                            end do
                            
                            if(connectedik)then!if there's a bridge between isite and ksite
                               
                               !---CALCULATE THE 3 MEMBERED RINGS---
                               !older version
                               !check if there are any common OH to ksite and jsite
                               !do p=1,4
                               !   do q=1,4
                               !      if(kfnsite(q).eq.jfnsite(p))rings3_temp = rings3_temp + 1
                               !   end do
                               !end do
                               !------------------------------------
                               
                               !---CALCULATE THE 4 MEMBERED RINGS
                               do bondj = 1,4!loop over first neighbors of jsite
                                  OH21 = jfnsite(bondj);occOH21 = occupancy(OH21)
                                  
                                  if((OH21.ne.OH11).and.((occOH21.eq.occSNOSNO).or.&
                                       (occOH21.eq.occSIOSNO)))then!if there's another bridge of jsite
                                     
                                     do bondOH21 = 1,8!loop over first neighbors of OH21
                                        lsite = n1list(bondOH21,OH21);locc = occupancy(lsite);lvertex = head(lsite)
                                        
                                        if((lvertex.gt.0).and.(lsite.ne.jsite).and.(lsite.ne.ksite))then!if there's SN/SI on lsite
                                           
                                           !find the OH of lsite
                                           do bondl = 1,4
                                              lfn(bondl) = Si_O(bondl,lvertex)
                                              lfnsite(bondl) = n1list(lfn(bondl),lsite)
                                           end do
                                           
                                           !loop through the first neighbors of OH21 to find the bridge between lsite and jsite
                                           connectedjl = .false.
                                           do bondl = 1,4
                                              if(lfnsite(bondl).eq.OH21)connectedjl = .true.
                                           end do
                                           
                                           if(connectedjl)then
                                              !check if there's any common OH between lsite and ksite
                                              do p = 1,4
                                                 do q=1,4
                                                    if((kfnsite(p).eq.lfnsite(q)).and.(kfnsite(p).ne.OH12))&
                                                         rings4_temp = rings4_temp + 1
                                                 end do
                                              end do
                                           end if
                                           
                                        elseif((lvertex.gt.0).and.(lsite.ne.jsite).and.(lsite.eq.ksite))then!if there's SN/SI on lsite
                                              
                                           !---CALCULATE 3 MEMBERED RINGS---
                                           !newer version
                                           !loop over OH of lsite and check if it is bonded to OH12
                                           do bondk = 1,4
                                              if(kfnsite(bondk).eq.OH12)rings3_temp = rings3_temp + 1
                                           end do
                                           !--------------------------------
                                           
                                        end if!endif: there's SN/SI on lsite
                                        
                                     end do!enddo: loop over first neighbors of OH21
                                                                          
                                  end if!endif: there's another bridge of jsite site
                                  
                               end do!enddo: loop over first neighbors of jsite
                               !---------------------------------
                               
                            end if!endif: there's a bridge between isite and ksite
                            
                         end if!endif: ksite has SN/SI
                                                  
                      end do!enddo: loop over first neighbors of OH12
                      
                   end if!endif: there's a bridge between isite and jsite
                   
                end if!endif jsite has SN/SI

             end do!enddo: loop over first neighbors of OH11
             
          end if!endif: there are more then 2 brigdes
          
       end do!enddo: loop over the OH of isite
    end do!enndo: loop over the OH of isite
    
    rings3 = rings3_temp/2
    rings4 = rings4_temp/2
    
    return
  end subroutine site_rings3and4
  !---------------------------------------------------------
  !---site_energy(isite,ivertex,overlep)
  !--------------------------------------------------------- 
  subroutine site_energy(isite,ivertex,ispin,rings3,rings4,en)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex,ispin
    real,intent(out) :: en
    integer,intent(out) :: rings3,rings4
    integer :: i,j,bond,bond2,jocc,ihead,jsite,ksite
    integer :: interSNSN,interSISN,interBOSDA,jneigh
    integer :: lrx,lry,lrz,jrx,jry,jrz,mrx,mry,mrz,mtag
    integer,dimension(4) :: fn, fnsite
    integer :: lsite, msite
    logical :: linked = .false.
    
    !set the accumulators to zero
    en = real(0)
    interSNSN = 0;interSISN = 0;interBOSDA = 0
    rings3 = 0;rings4 = 0
    
    !if isite has SN or SI, then calculate its neighbors
    if((ispin.eq.spinSI).or.(ispin.eq.spinSN))then
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       ihead = fnsite(1)
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
       end do
    end if
   
    if(ispin.eq.spinSN)then!if isite has SN
       
       !loop over OHs and search for a bridge with SN or SI
       do bond = 1,4
          jocc = occupancy(fnsite(bond))
          if(jocc.eq.occSNOSNO)then
             interSNSN = interSNSN + 1!SN SN condensation
          elseif(jocc.eq.occSIOSNO)then
             interSISN = interSISN + 1!SI SN condensation
          end if
       end do

       !now setup (sda_size+1)th neighbor interactions of -O- with SDA
       do bond = 1,4!loop OH group of Si
          if(occupancy(fnsite(bond)).eq.occSNOSNO)then!if there's -O-
             jneigh = neighbors(sda_size+1)
             do j = 1,jneigh
                jsite = nlist(sda_size+1,j,fnsite(bond))
                if(occupancy(jsite).eq.occSDA)interBOSDA = interBOSDA + 1
             end do
          end if!endif: there's -O-
       end do!enddo: loop OH group of Si
       
    elseif(ispin.eq.spinSI)then!if isite has SI
       
       !loop over OHs and search for a bridge with SN
       do bond = 2,4
          jocc = occupancy(fnsite(bond))
          if(jocc.eq.occSIOSNO)interSISN = interSISN + 1
       end do

       !search for SI-SDA linkage
       !lsite = SIlink(ihead)
       !if(lsite.gt.0)then
       !   msite = SDAlink(lsite)
       !   !if(msite.eq.ihead)interISISDA = interISISDA + 1
       !   if(msite.eq.ihead)interISISDA = 1
       !end if
                    
    elseif(ispin.eq.spinSDA)then!if isite has SDA
       
       !setup interactions of of SDA with -O- at (sda_size+1)th neighbor
       jneigh = neighbors(sda_size+1)
       do j = 1,jneigh
          jsite = nlist(sda_size+1,j,isite)
          if(occupancy(jsite).eq.occSNOSNO)interBOSDA = interBOSDA + 1
       end do
       

       !search for SI-SDA linkage
       !msite = SDAlink(isite)
       !if(msite.gt.0)then
          !lsite = SIlink(msite)
          !if(lsite.eq.isite)interISISDA = interISISDA + 1
          !if(lsite.eq.isite)interISISDA = 1
          !interISISDA = 1
       !end if
       
    else!otherwise
       print*,'BUG :: isite has neither SI/SN/SDA :: subroutine site_energy'
       print*,'jsite:',isite
       print*,'jvertex:',ivertex
       print*,'jspin:',ispin
       call freeandclose
       stop
    end if!endif: isite has SN
    
    !calculate the 3 rings and 4 rings the site is associated with
    if((ispin.eq.spinSN).or.(ispin.eq.spinSI))then
       call site_rings3and4(isite,ivertex,rings3,rings4)
    else
       rings3=0;rings4=0
    end if

    en = eSNSN*real(interSNSN) + eSISN*real(interSISN) + eOSDA*real(interBOSDA) + &
         pen3*real(rings3) + pen4*real(rings4)
    
    return
  end subroutine site_energy
  !---------------------------------------------------------
  !--- total_energy: CALCULATE TOTAL ENERGY OF THE LATTICE
  !---------------------------------------------------------
  subroutine total_energy(tot_en)
    use variables
    implicit none
        
    real,intent(out) :: tot_en
    integer :: isite,ivertex,ispin,ihead
    integer :: i,j,molecule,bond,jocc,jsite
    integer :: rings3,rings4,jneigh
    integer :: interSNSN,interSISN,interBOSDA
    integer :: lrx,lry,lrz,jrx,jry,jrz,mrx,mry,mrz,mtag
    integer,dimension(4) :: fn,fnsite
    integer :: lsite, msite    
    
    call rings3and4(rings3,rings4)
    
    !set the accumulators to zero
    tot_en = real(0)
    interSNSN = 0;interSISN = 0;interBOSDA = 0
        
    do molecule = 1,allmolecules!loop over all molecules
       
       isite = noccupy(molecule)
       ispin = spin(isite)
       ivertex = head(isite)
       msite = -1;lsite = -1
       
       if((ispin.eq.spinSI).or.(ispin.eq.spinSN))then
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          ihead = fnsite(1)
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)
          end do
       end if
       
       if(ispin.eq.spinSN)then!if isite has SN
          
          !loop over OHs and search for a bridge with SN or SI
          do bond = 1,4
             jocc = occupancy(fnsite(bond))
             if(jocc.eq.occSNOSNO)then
                interSNSN = interSNSN + 1!SN SN condensation
             elseif(jocc.eq.occSIOSNO)then
                interSISN = interSISN + 1!SI SN condensation
             end if
          end do
          
          !include the -O-=SDA interaction with SDA molecule not with SN to 
          !eliminate double counting
          
       elseif(ispin.eq.spinSI)then!if isite has SI
          
          !interactions taken care by SN and SDA
          !loop over OHs and search for a bridge with SN
          !do bond = 2,4
          !   jocc = occupancy(fnsite(bond))
          !   if(jocc.eq.occSIOSNO)interSISN = interSISN + 1
          !end do

          !search for SI-SDA linkage
          !lsite = SIlink(ihead)
          !if(lsite.gt.0)then
          !   msite = SDAlink(lsite)
          !   if(msite.eq.ihead)interISISDA = interISISDA + 1
          !end if
                    
       elseif(ispin.eq.spinSDA)then!if isite has SDA

          !setup interactions of of SDA with -O- at (sda_size+1)th neighbor
          jneigh = neighbors(sda_size+1)
          do j = 1,jneigh
             jsite = nlist(sda_size+1,j,isite)
             if(occupancy(jsite).eq.occSNOSNO)interBOSDA = interBOSDA + 1
          end do

          !setup interactions of of SDA with -O-
          ! do bond = 1,6!loop over second neighbors of Si
          !    jsite = n2list(bond,isite)
          !    if(occupancy(jsite).eq.occSNOSNO)interBOSDA = interBOSDA + 1
          ! end do!enddo: loop over second neighbors of Si
          
          !search for SI-SDA linkage
          !msite = SDAlink(isite)
          !if(msite.gt.0)then
          !   lsite = SIlink(msite)
          !   if(lsite.eq.isite)interISISDA = interISISDA + 1
          !   interISISDA = interISISDA + 1
          !end if
                    
       else!otherwise
          print*,'BUG :: isite has neither SI/SN/SDA :: subroutine site_energy'
          print*,'jsite:',isite
          print*,'jvertex:',ivertex
          print*,'jspin:',ispin
          call freeandclose
          stop
       end if!endif: isite has SN
       
    end do!enddo:loop over all molecules
    
    !calculate the energy enforced by the 3/4 membered ring penalties
    !tot_en = real(0.5)*(eSNSN*real(interSNSN) + eSISN*real(interSISN) + &
    !     eISISDA*real(interISISDA)) + pen3*real(rings3) + pen4*real(rings4)
    
    tot_en = real(0.5)*eSNSN*real(interSNSN) + eSISN*real(interSISN) + &
         eOSDA*real(interBOSDA) + pen3*real(rings3) + pen4*real(rings4)
    
    return
  end subroutine total_energy
  !---------------------------------------------------------
  !---insert_SN(isite,ivertex,ispin): INSERT SN at isite
  !---------------------------------------------------------   
  subroutine insert_site(isite,ivertex,ispin,molecule)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex,ispin,molecule
    integer :: bond,occfn
    integer,dimension(4) :: fn,fnsite
    
    !update the occupancy of isite
    if(ispin.eq.spinSN)then

       occupancy(isite) = occSN
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       occupancy(fnsite(1)) = occupancy(fnsite(1)) + occSNO
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          occupancy(fnsite(bond)) = occupancy(fnsite(bond)) + occSNO
       end do

    elseif(ispin.eq.spinSI)then

       occupancy(isite) = occSI
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       occupancy(fnsite(1)) = occISIO
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          occupancy(fnsite(bond)) = occupancy(fnsite(bond)) + occSIO
       end do
       
    elseif(ispin.eq.spinSDA)then

       occupancy(isite) = occupancy(isite) + occSDA

    else
       print*,'BUG:molecule is not SN/SI/SDA :: subroutine insert_site'
       print*,'ispin:',ispin
       call freeandclose
       stop
    end if
    
    !update head,spin,noccupy,resite arrays
    if(ispin.eq.spinSI.or.ispin.eq.spinSN)then
       head(isite) = ivertex
    elseif(ispin.eq.spinSDA)then
       head(isite) = 0
    end if
    spin(isite) = ispin
    noccupy(molecule) = isite
    resite(isite) = molecule

    return
  end subroutine insert_site
  !---------------------------------------------------------
  !---remove_SN(isite,ivertex): REMOVE SN from isite
  !---------------------------------------------------------   
  subroutine remove_site(isite,ivertex,ispin,molecule)
    use variables
    implicit none
    
    integer,intent(in) :: isite, ivertex,ispin,molecule
    integer :: bond,occfn
    integer,dimension(4) :: fn,fnsite
    
    !find the 1st neighbors of isite and update their occupancy
    if(ispin.eq.spinSN)then

       occupancy(isite) = occupancy(isite) - occSN
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       occupancy(fnsite(1)) = occupancy(fnsite(1)) - occSNO
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          occupancy(fnsite(bond)) = occupancy(fnsite(bond)) - occSNO
       end do

    elseif(ispin.eq.spinSI)then
       
       occupancy(isite) = occupancy(isite) - occSI
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       occupancy(fnsite(1)) = occupancy(fnsite(1)) - occISIO
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          occupancy(fnsite(bond)) = occupancy(fnsite(bond)) - occSIO
       end do

    elseif(ispin.eq.spinSDA)then
       
       occupancy(isite) = occupancy(isite) - occSDA
       
    else
       print*,'BUG:molecule is not SN/SI/SDA :: subroutine remove_site'
       print*,'ispin:',ispin
       call freeandclose
       stop
    end if

    !update the spin and head arrays
    !note: resite and noccupy arrays are not updated
    spin(isite) = 0
    head(isite) = 0
    
    !note:noccupy/resite arrays are not updated
    !they are updated in the subroutine insert_site
    return
  end subroutine remove_site
  !---------------------------------------------------------
  !---mc_algo: Parallel Tempering MONTE CARLO ALGORITHM
  !--------------------------------------------------------- 
  subroutine mc_algo(comm_rank,comm_size)
    use variables
    use MPI
    implicit none
    
    integer,intent(in) :: comm_rank,comm_size
    integer :: molecule,bond, lsite
    integer :: isite,ispin,ivertex,imolecule
    integer :: jsite,jspin,jvertex,jmolecule,jneigh
    integer :: ksite,kspin,kvertex,kmolecule
    !integer :: occfn,overlap_old,overlap_new
    !integer :: rings3_old,rings3_new,rings3_temp,rings4_old,rings4_new,rings4_temp
    integer :: rings3,rings4,ierror
    !real :: d0
    real :: old_en,new_en,delta_en,tot_en,temp_en,norm_en
    real :: AR1,AR2,AR3
    real,dimension(0:4) :: norm_Qn
    
    !calculate total energy, Qn distribution and 3 and 4 ring 
    !distribution before simulation
    call Qn_dist
    call total_energy(tot_en)
    call rings3and4(rings3,rings4)
    !call cluster_stat
    !call average_cluster_dimension(d0)
    norm_Qn = real(Qn)/real(size(noccupy))
    norm_en = tot_en/real(nsites)
    c = (norm_Qn(1) + norm_Qn(2)*2.0d0 + norm_Qn(3)*3.0d0 + norm_Qn(4)*4.0d0)/4.0d0
    if(initial_step.eq.1)then
       write(1001,*) step,norm_Qn,c;flush(1001)
       write(1002,*) step,rings3,rings4;flush(1002)
       write(1003,*) step,tot_en,norm_en;flush(1003)
       !write(1007,*) step,max_size,avg_size,cluster_num,d0;flush(1007)
    end if
 
    do step = initial_step,nsweeps!do all MC steps
          
       do molecule = 1,allmolecules!loop over all molecules
          
          !*********************************START SWAP MOVES*********************************          
          imolecule = int(1 + ran3(seed)*real(allmolecules))
          isite = noccupy(imolecule)
          ivertex = head(isite)
          ispin = spin(isite)

          AR1 = real(0); AR2 = real(0); AR3 = real(0)
          
          !select a random site
          !SN: select a random site anywhere on the lattice
          !SDA: select a random site out of the first 8 nearest neighbors
          jsite = int(1 + ran3(seed)*real(nsites))
          do while(jsite.eq.isite)
             jsite = int(1 + ran3(seed)*real(nsites))
          end do
          jvertex = head(jsite)
          jspin = spin(jsite)
          jmolecule = resite(jsite)
          
          ! nearest neighbor of translation for SDA
          lsite = int(1 + ran3(seed)*real(8))
          lsite = n1list(lsite,isite)
          
          
          if(occupancy(isite).eq.occSN)then
             
             if(occupancy(jsite).eq.occW)then
                call translate(isite,ivertex,ispin,imolecule,jsite,rings3,rings4,tot_en,AR1)
                SN_jumps_attempt = SN_jumps_attempt + 1
             elseif(occupancy(jsite).eq.occSI)then
                call swap(jsite,jvertex,jspin,jmolecule,isite,ivertex,ispin,imolecule,&
                     rings3,rings4,tot_en,AR2)
             elseif(occupancy(jsite).eq.occSDA)then
                call swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,&
                     rings3,rings4,tot_en,AR2)
                SN_SDA_swaps_attempt = SN_SDA_swaps_attempt + 1
             end if
             
          ! elseif(occupancy(isite).eq.occSI)then
             
          !    if(occupancy(jsite).eq.occW)then
          !       call translate(isite,ivertex,ispin,imolecule,jsite,rings3,rings4,tot_en,AR1)
          !    elseif(occupancy(jsite).eq.occSDA)then
          !       call swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,&
          !            rings3,rings4,tot_en,AR2)
          !    elseif(occupancy(jsite).eq.occSN)then
          !       !call swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,&
          !       !     rings3,rings4,tot_en,AR2)
          !    end if
             
          elseif(occupancy(isite).eq.occSDA)then
             
             if(occupancy(lsite).eq.occW)then
                call translate(isite,ivertex,ispin,imolecule,lsite,rings3,rings4,tot_en,AR1)
                SDA_translations_attempt = SDA_translations_attempt + 1
             elseif(occupancy(jsite).eq.occSI)then
                call swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,&
                     rings3,rings4,tot_en,AR2)
             elseif(occupancy(jsite).eq.occSN)then
                call swap(jsite,jvertex,jspin,jmolecule,isite,ivertex,ispin,imolecule,&
                     rings3,rings4,tot_en,AR2)                
                SN_SDA_swaps_attempt = SN_SDA_swaps_attempt + 1
             end if
             
          else
             print*,'BUG :: molecule neither SI/SN/SDA :: subroutine mc_algo'
             print*,'isite:',isite
             print*,'ivertex:',ivertex
             print*,'imolecule:',imolecule
             print*,'ispin:',ispin
             call freeandclose
             stop
          end if
          !*********************************END SWAP MOVES***********************************
          
          !*********************************START ROTATION MOVES*****************************
          !select a monomer amongst all
          kmolecule = int(1 + ran3(seed)*real(allmolecules))
          ksite = noccupy(kmolecule)
          kspin = spin(ksite)
          kvertex = head(ksite)
          if((occupancy(ksite).eq.occSI).or.(occupancy(ksite).eq.occSN))then!if ksite has SI/SN on it
             call rotate(ksite,kvertex,kspin,kmolecule,rings3,rings4,tot_en,AR3)
             SN_rotations_attempt = SN_rotations_attempt + 1
          elseif(kspin.eq.spinSDA)then
             !reject the move
          else
             print*,'BUG :: molecule for rotation neither SI/SN/SDA :: subroutine mc_algo'
             print*,'isite:',ksite
             print*,'ivertex:',kvertex
             print*,'imolecule:',kmolecule
             print*,'ispin:',kspin     
             call freeandclose
             stop
          end if!endif: ksite has SI/SN on it
          !*********************************END ROTATION MOVES*******************************
          
       end do!endloop over all molecules
       
       !print the output of the program
       call print_output(AR1,AR2,AR3,rings3,rings4,tot_en,comm_rank)
              
       if((step.le.int(0.90*nsweeps)).and.(mod(step,nexchange).eq.0))then
          call pt_move(comm_rank,comm_size,tot_en)       
       end if
       
    end do!enddo all MC steps
    
    return
  end subroutine mc_algo
  !---------------------------------------------------------
  !---PARALLEL TEMPERING MOVE
  !--------------------------------------------------------- 
  subroutine pt_move(comm_rank,comm_size,tot_en)
    use variables
    use MPI
    implicit none
    
    integer,intent(in) :: comm_rank,comm_size
    real,intent(in) :: tot_en
    integer :: ierror,i,j,krank,error,sort_sweeps = 100
    integer :: jump_probability,max_jump,jump,ijump_rank
    real :: rand,PT_AR,temp1,temp2
    real :: itstar, jtstar, ratio,epsilon = real(0.00001)
    integer,dimension(1:MPI_status_size) :: status
    real,dimension(2,0:comm_size-1) :: dyn_temp
    integer,dimension(0:comm_size-1) :: sorted_pid
    
    !make the process synchronous by making all the processors wait till they reach here
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    
    !*********************************************************************************************
    !to ensure that only the adjacent temperatures are swapped, we first create an array
    !of processors IDs with increasing temperature

    !send all the temperatures from the processors to the master processor(PID = 0)
    if(comm_rank.gt.0)call MPI_SSEND(tstar,1,MPI_REAL,0,comm_rank,MPI_COMM_WORLD,ierror)
    
    if(comm_rank.eq.0)then
       dyn_temp(1,0) = real(0)
       dyn_temp(2,0) = tstar
       do i = 1,comm_size - 1
          call MPI_RECV(dyn_temp(2,i),1,MPI_REAL,i,i,MPI_COMM_WORLD,status,ierror)    
          dyn_temp(1,i) = i
       end do
    end if
    
    ! if(comm_rank.eq.0)then
    !    write(1000000,*)
    !    write(1000000,*)'UNSORTED ARRAY'
    !    do i = 0,comm_size - 1
    !       write(1000000,*)dyn_temp(1,i),dyn_temp(2,i),i,proc_temp(i)
    !    end do
    ! end if
    
    !now sort the dyn array ar master id
    if(comm_rank.eq.0)then
       do i = 1,sort_sweeps
          do j = 0,comm_size - 2
             
             if(dyn_temp(2,j).gt.(dyn_temp(2,j+1)))then
                temp1 = dyn_temp(2,j)
                dyn_temp(2,j) = dyn_temp(2,j+1)
                dyn_temp(2,j+1) = temp1
                temp1 = dyn_temp(1,j)
                dyn_temp(1,j) = dyn_temp(1,j+1)
                dyn_temp(1,j+1) = temp1
             end if
             
          end do
       end do
    end if
    
    ! if(comm_rank.eq.0)then
    !    write(1000000,*)
    !    write(1000000,*)'SORTED ARRAY'
    !    do i = 0,comm_size - 1
    !       write(1000000,*)dyn_temp(1,i),dyn_temp(2,i),i,proc_temp(i)
    !    end do
    ! end if
    
    !now make the integral sorted_pid array
    if(comm_rank.eq.0)then
       do i = 0,comm_size - 1
          sorted_pid(i) = int(dyn_temp(1,i))
       end do
    end if
    
    ! if(comm_rank.eq.0)then
    !    write(2000000,*)'SORTED PID ARRAY'
    !    do i = 0,comm_size - 1
    !       write(2000000,*)sorted_pid(i)
    !    end do
    ! end if
    
    !*********************************************************************************************    
    
    max_jump = int(int(comm_size)/int(4))
    
    if(mod(step,nexchange).eq.0)then!if step is a multiple of nexchange

       do krank = 0,comm_size-2!loop over all processors,
          
          !first the master selects the 2 processors
          if(comm_rank.eq.0)then!if the processor is a master processor
             iproc = sorted_pid(krank)
             jproc = sorted_pid(krank+1)
             
             ! ----------------------
             ! ---HIGH TEMP EXCHANGE-
             ! ----------------------
             ! here we associate a higher temperature exchange (other than neighboring)
             ! 10% of the total moves go to exchanging high 
             ! The Swaps are made between temperatures selected from 0 to 25% of comm_size
             jump_probability = int(1 + ran3(seed)*real(10))
             if(jump_probability.eq.1)then!if the jump probability is 10%
                
                !select a random temperature distance to jump 1->max_jump
                jump = int(1 + ran3(seed)*real(max_jump))
                
                !select a random processor to make the jump
                ijump_rank = int(ran3(seed)*real(comm_size-1-jump))
                iproc = sorted_pid(ijump_rank)
                jproc = sorted_pid(ijump_rank+jump)
                
             end if!endif: the jump probability is 10%
             ! ----------------------
             
          end if!endif: the processor is a master processor
          !call MPI_BARRIER(MPI_COMM_WORLD,ierror)       
          !broadcast iproc,jproc variables to all processors from master
          call MPI_BCAST(iproc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(jproc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
          
          !since the master processor is doing the calculations,
          !three cases arise
          !iproc = master; jproc = slave
          !iproc = slave; jproc = master
          !iproc = slave; jproc = slave
          
          if((iproc.eq.0).and.(jproc.gt.0))then!if iproc is master and jproc is slave

             !---------------------CASE I-----------------------------------------------
             !--------------------------------------------------------------------------
             !send temperature and energy from jproc to iproc => slave to master
             !---send signal from jproc to iproc
             if(comm_rank.eq.jproc)then
                call MPI_SSEND(tstar,1,MPI_REAL,0,1000,MPI_COMM_WORLD,ierror)!temperature
                call MPI_SSEND(tot_en,1,MPI_REAL,0,2000,MPI_COMM_WORLD,ierror)!energy
             end if
             !---receive data to master from jproc
             if(comm_rank.eq.0)then
                call MPI_RECV(procj_temp,1,MPI_REAL,jproc,1000,MPI_COMM_WORLD,status,ierror)!temperature
                call MPI_RECV(procj_en,1,MPI_REAL,jproc,2000,MPI_COMM_WORLD,status,ierror)!energy
                proci_temp = tstar
                proci_en = tot_en
             end if
             !call MPI_BARRIER(MPI_COMM_WORLD,ierror)          
             
             !-----CONDITION FOR PARALLEL TEMPERINGO ON MASTER-----
             if(comm_rank.eq.0)then!if the processor is a master processor
                PT_AR = exp( (proci_en - procj_en) * (real(1/proci_temp) - real(1/procj_temp)) )
                rand = ran3(seed)
                if(PT_AR.ge.rand)then!IF THE MOVE IS ACCEPTED
                   !exchange the temperatures and update the temperature array
                   temp1 = proci_temp
                   proci_temp = procj_temp
                   procj_temp = temp1
                   proc_temp(iproc) = proci_temp
                   proc_temp(jproc) = procj_temp
                else!IF THE MOVE IS REJECTED
                   !do nothing
                end if!ENDIF: THE MOVE IS ACCEPTED
             end if!endif: processor is a master processor
             
             !send the temperatures back to the jproc from master
             !---send signal from master
             if(comm_rank.eq.0)then
                tstar = proci_temp
                call MPI_SSEND(procj_temp,1,MPI_REAL,jproc,1300,MPI_COMM_WORLD,ierror)
             end if
             !receive signal to slave
             if(comm_rank.eq.jproc)call MPI_RECV(tstar,1,MPI_REAL,0,1300,MPI_COMM_WORLD,status,ierror)
             !--------------------------------------------------------------------------
             !--------------------------------------------------------------------------
          
          elseif((iproc.gt.0).and.(jproc.eq.0))then!if jproc is master and iproc is slave

             !----------------------CASE II---------------------------------------------
             !--------------------------------------------------------------------------
             !send temperature and energy from iproc to jproc => slave to master
             !---send signal from jproc to iproc
             if(comm_rank.eq.iproc)then
                call MPI_SSEND(tstar,1,MPI_REAL,0,1000,MPI_COMM_WORLD,ierror)!temperature
                call MPI_SSEND(tot_en,1,MPI_REAL,0,2000,MPI_COMM_WORLD,ierror)!energy
             end if
             !---receive data to master from iproc
             if(comm_rank.eq.0)then
                call MPI_RECV(proci_temp,1,MPI_REAL,iproc,1000,MPI_COMM_WORLD,status,ierror)!temperature
                call MPI_RECV(proci_en,1,MPI_REAL,iproc,2000,MPI_COMM_WORLD,status,ierror)!energy
                procj_temp = tstar
                procj_en = tot_en
             end if
             !call MPI_BARRIER(MPI_COMM_WORLD,ierror)     
             
             !-----CONDITION FOR PARALLEL TEMPERINGO ON MASTER-----
             if(comm_rank.eq.0)then!if the processor is a master processor
                PT_AR = exp( (proci_en - procj_en) * (real(1/proci_temp) - real(1/procj_temp)) )
                rand = ran3(seed)
                if(PT_AR.ge.rand)then!IF THE MOVE IS ACCEPTED
                   !exchange the temperatures and update the temperature array
                   temp1 = proci_temp
                   proci_temp = procj_temp
                   procj_temp = temp1
                   proc_temp(iproc) = proci_temp
                   proc_temp(jproc) = procj_temp
                else!IF THE MOVE IS REJECTED
                   !do nothing
                end if!ENDIF: THE MOVE IS ACCEPTED
             end if!endif: processor is a master processor
             
             !send the temperatures back to the iproc from master
             !---send signal from master
             if(comm_rank.eq.0)then
                tstar = procj_temp
                call MPI_SSEND(proci_temp,1,MPI_REAL,iproc,1300,MPI_COMM_WORLD,ierror)
             end if
             !receive signal to slave
             if(comm_rank.eq.iproc)call MPI_RECV(tstar,1,MPI_REAL,0,1300,MPI_COMM_WORLD,status,ierror)        
             !--------------------------------------------------------------------------
             !--------------------------------------------------------------------------

          else!if iproc and jproc both are slaves

             !----------------------CASE III--------------------------------------------
             !--------------------------------------------------------------------------
             !iproc and jproc first send their energy and temperature to master
             !---send signal from jproc to iproc
             if(comm_rank.eq.iproc)then
                call MPI_SSEND(tstar,1,MPI_REAL,0,1000,MPI_COMM_WORLD,ierror)!temperature
                call MPI_SSEND(tot_en,1,MPI_REAL,0,2000,MPI_COMM_WORLD,ierror)!energy
             end if
             if(comm_rank.eq.jproc)then
                call MPI_SSEND(tstar,1,MPI_REAL,0,3000,MPI_COMM_WORLD,ierror)!temperature
                call MPI_SSEND(tot_en,1,MPI_REAL,0,4000,MPI_COMM_WORLD,ierror)!energy
             end if
             !---receive data to master from iproc and jproc
             if(comm_rank.eq.0)then
                !from iproc
                call MPI_RECV(proci_temp,1,MPI_REAL,iproc,1000,MPI_COMM_WORLD,status,ierror)!temperature
                call MPI_RECV(proci_en,1,MPI_REAL,iproc,2000,MPI_COMM_WORLD,status,ierror)!energy
                !from jproc
                call MPI_RECV(procj_temp,1,MPI_REAL,jproc,3000,MPI_COMM_WORLD,status,ierror)!temperature
                call MPI_RECV(procj_en,1,MPI_REAL,jproc,4000,MPI_COMM_WORLD,status,ierror)!energy
             end if
             !call MPI_BARRIER(MPI_COMM_WORLD,ierror)    
             
             !-----CONDITION FOR PARALLEL TEMPERINGO ON MASTER-----
             if(comm_rank.eq.0)then!if the processor is a master processor
                PT_AR = exp( (proci_en - procj_en) * (real(1/proci_temp) - real(1/procj_temp)) )
                rand = ran3(seed)
                if(PT_AR.ge.rand)then!IF THE MOVE IS ACCEPTED
                   !exchange the temperatures and update the temperature array
                   temp1 = proci_temp
                   proci_temp = procj_temp
                   procj_temp = temp1
                   proc_temp(iproc) = proci_temp
                   proc_temp(jproc) = procj_temp
                else!IF THE MOVE IS REJECTED
                   !do nothing
                end if!ENDIF: THE MOVE IS ACCEPTED
             end if!endif: processor is a master processor
             
             !send the temperatures back from master to iproc and jproc
             !---send signal from master
             if(comm_rank.eq.0)then
                call MPI_SSEND(proci_temp,1,MPI_REAL,iproc,5000,MPI_COMM_WORLD,ierror)
                call MPI_SSEND(procj_temp,1,MPI_REAL,jproc,6000,MPI_COMM_WORLD,ierror)
             end if
             !receive signal to slave
             if(comm_rank.eq.iproc)call MPI_RECV(tstar,1,MPI_REAL,0,5000,MPI_COMM_WORLD,status,ierror)        
             if(comm_rank.eq.jproc)call MPI_RECV(tstar,1,MPI_REAL,0,6000,MPI_COMM_WORLD,status,ierror)        
             !--------------------------------------------------------------------------
             !--------------------------------------------------------------------------

          end if
          
       end do!enddo: loop over all processors
       
       !-----WRITE OUT THE EXCHANGE DETAIL FILE-----
       if(comm_rank.eq.0)then
          
          if(mod(step,nprint).eq.0)then                
             write(5000,*)'########''STEP:',step,'###################';flush(5000)
             if(PT_AR.ge.rand)then
                write(5000,*) 'Exchange ACCEPTED';flush(5000)
             else
                write(5000,*) 'Exchange REJECTED';flush(5000)
             end if
             write(5000,*) 'Acceptance ratio,random number',PT_AR,rand;flush(5000)
             write(5000,*) 'iproc,iproc_temp',iproc,proci_temp,proci_en;flush(5000)
             write(5000,*) 'jproc,jproc_temp',jproc,procj_temp,procj_en;flush(5000)
             write(5000,*)'###########################';flush(5000)
             !write the temperature out on the temperature file
             write(5000,*)'---------------------------';flush(5000)
             write(5000,*)'---------------------------';flush(5000)
             do i=0,comm_size-1
                write(5000,*) i,proc_temp(i);flush(5000)
             end do
          end if
       end if
       
    else!step is not a multiple of nexchange
       return
    end if!endif: step is a multiple of nexchange  
          
    return
  end subroutine pt_move
  !---------------------------------------------------------
  !---translate: translates molecule from isite to jsite(water site)
  !--------------------------------------------------------- 
  subroutine translate(isite,ivertex,ispin,imolecule,jsite,rings3,rings4,tot_en,AR)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex,ispin,imolecule,jsite
    integer,intent(inout) :: rings3,rings4
    real,intent(inout) :: tot_en
    real,intent(out) :: AR
    integer :: i,rings3_old,rings3_new,rings4_old,rings4_new
    real :: old_en,new_en,delta_en,rand
    logical :: can_insert, neighbor

    !first calculate the molecule energy and them remove the molecule
    !calculate the energy associated with the molecule at isite
    call site_energy(isite,ivertex,ispin,rings3_old,rings4_old,old_en)
    !remove molecule from isite and insert it at jsite
    call remove_site(isite,ivertex,ispin,imolecule)
    
    !first check if the molecule can be inserted at jsite
    if(ispin.eq.spinSN)then
       call can_insert_SN(jsite,ivertex,can_insert)
    elseif(ispin.eq.spinSI)then
       call can_insert_SI(jsite,ivertex,can_insert)
    elseif(ispin.eq.spinSDA)then
       call can_insert_SDA(jsite,can_insert)
    else
       print*,'BUG :: molecule neither SI/SN/SDA :: subroutine translate'
       print*,'isite:',isite
       print*,'ivertex:',ivertex
       print*,'imolecule:',imolecule
       print*,'ispin:',ispin
    end if
    
    if(can_insert)then!if molecule at isite can be inserted at jsite
       
       !insert the molecule at jsite and calculate the new energy associated with it
       call insert_site(jsite,ivertex,ispin,imolecule)
       call site_energy(jsite,ivertex,ispin,rings3_new,rings4_new,new_en)
       
       !calculate the change in energy and apply the metropolis criteria
       delta_en = new_en - old_en
       AR = exp(-1.000d0*delta_en/tstar)
       rand = ran3(seed)
       if(AR.gt.rand)then!if metropolis criteria is satisfied
          tot_en = tot_en + delta_en
          rings3 = rings3 + (rings3_new - rings3_old)
          rings4 = rings4 + (rings4_new - rings4_old)
          
          !update the move variables
          if(ispin.eq.spinSN)then
             SN_jumps_success = SN_jumps_success + 1
          elseif(ispin.eq.spinSI)then
             !do nothing
          elseif(ispin.eq.spinSDA)then
             SDA_translations_success = SDA_translations_success + 1
          end if
          
       else!if metropolis criteria is not satisfied
          !remove the molecule from jsite and put it back on isite
          call remove_site(jsite,ivertex,ispin,imolecule)
          call insert_site(isite,ivertex,ispin,imolecule)
       end if!endif:metropolis criteria is satisfied
       
    else
       !remove the molecule from jsite and put it back on isite
       call insert_site(isite,ivertex,ispin,imolecule)
    end if!if molecule at isite can be inserted at jsite
    
    !safety check
    !call safetycheck_translate(tot_en,ispin,AR,rand,old_en,new_en,rings3,rings4)
    !call linkage
    
    return
  end subroutine translate
  !---------------------------------------------------------
  !---swap: swaps two molecules
  !---note: molecule at jsite would be removed first 
  !---------------------------------------------------------   
  subroutine swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,rings3,rings4,tot_en,AR)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex,ispin,imolecule
    integer,intent(in) :: jsite,jvertex,jspin,jmolecule
    integer,intent(inout) :: rings3,rings4
    real,intent(inout) :: tot_en
    real,intent(out) :: AR
    integer :: rings3_old,rings3_new,rings4_old,rings4_new,rings3_temp,rings4_temp
    real :: old_en,new_en,delta_en,new_temp_en,old_temp_en,rand
    logical :: can_inserti, can_insertj
    !variables for checking double couting in energy
    integer :: ihead,jhead,hrx,hry,hrz,irx,iry,irz,jrx,jry,jrz
    real :: d_ij,d_ih,d_jh,d_ij_sq,d_ih_sq,d_jh_sq
    
        
    !calculate the site energy associated with the molecule located at isite
    call site_energy(jsite,jvertex,jspin,rings3_temp,rings4_temp,old_temp_en)
    !remove the molecule from isite
    call remove_site(jsite,jvertex,jspin,jmolecule)
    
    !check if the molecule at isite can be inserted at jsite
    if(ispin.eq.spinSN)then
       call can_insert_SN(jsite,ivertex,can_inserti)
    elseif(ispin.eq.spinSI)then
       call can_insert_SI(jsite,ivertex,can_inserti)
    elseif(ispin.eq.spinSDA)then
       call can_insert_SDA(jsite,can_inserti)
    else
       print*,'BUG :: imolecule neither SI/SN/SDA :: subroutine swap'
       print*,'isite:',isite
       print*,'ivertex:',ivertex
       print*,'imolecule:',imolecule
       print*,'ispin:',ispin
    end if
    
    if(can_inserti)then!if the molecule at isite can be inserted at jsite
       
       !calculate the energy associated with the molecule located at isite
       call site_energy(isite,ivertex,ispin,rings3_old,rings4_old,old_en)
       old_en = old_en + old_temp_en
       rings3_old = rings3_old + rings3_temp
       rings4_old = rings4_old + rings4_temp
       !now remove imolecule from isite
       call remove_site(isite,ivertex,ispin,imolecule)
       
       !now check if jmolecule can be inserted at isite
       if(jspin.eq.spinSN)then
          call can_insert_SN(isite,jvertex,can_insertj)
       elseif(jspin.eq.spinSI)then
          call can_insert_SI(isite,jvertex,can_insertj)
       elseif(jspin.eq.spinSDA)then
          call can_insert_SDA(isite,can_insertj)
       else
          print*,'BUG :: jmolecule neither SI/SN/SDA :: subroutine swap'
          print*,'jsite:',jsite
          print*,'jvertex:',jvertex
          print*,'jmolecule:',jmolecule
          print*,'jspin:',jspin
       end if

       if(can_insertj)then!if jmolecule can be inserted at isite

          !put imolecule on jsite and calculate its energy
          call insert_site(jsite,ivertex,ispin,imolecule)
          call site_energy(jsite,ivertex,ispin,rings3_temp,rings4_temp,new_temp_en)
          !put imolecule on isite and calculate its energy
          call insert_site(isite,jvertex,jspin,jmolecule)
          call site_energy(isite,jvertex,jspin,rings3_new,rings4_new,new_en)

          ! calculate the new energy
          new_en = new_en + new_temp_en
          rings3_new = rings3_new + rings3_temp
          rings4_new = rings4_new + rings4_temp
          SN_SDA_swaps_success = SN_SDA_swaps_success + 1
          
          !calculate the change in energy after the move and apply the metropolis criteria
          delta_en = new_en - old_en
          AR = exp(-1.000000d0*delta_en/tstar)
          rand = ran3(seed)
          if(AR.gt.rand)then!if the netropolis criteria is satisfied
             tot_en = tot_en + delta_en
             rings3 = rings3 + (rings3_new - rings3_old)
             rings4 = rings4 + (rings4_new - rings4_old)
          else!if the netropolis criteria is not satisfied
             !remove imolecule from jsite and jmolecule from isite
             call remove_site(jsite,ivertex,ispin,imolecule)
             call remove_site(isite,jvertex,jspin,jmolecule)
             !put the molecules back in their original places
             call insert_site(jsite,jvertex,jspin,jmolecule)
             call insert_site(isite,ivertex,ispin,imolecule)
          end if!endif: the netropolis criteria is satisfied
                    
       else!if jmolecule cannot be inserted at isite

          !put both the molecules on their original location
          call insert_site(jsite,jvertex,jspin,jmolecule)
          call insert_site(isite,ivertex,ispin,imolecule)

       end if!endif: jmolecule can be inserted at isite
       
    else!if the molecule at isite cannot be inserted at jsite
       !put the jmolecule back on jsite
       call insert_site(jsite,jvertex,jspin,jmolecule)
    end if!endif: the molecule at isite can be inserted at jsite
    
    !safety check
    !call safetycheck_swap(tot_en,isite,jsite,AR,rand,old_en,new_en,rings3,rings4)
        
    return
  end subroutine swap
  !---------------------------------------------------------
  !---rotate: rotates the molecule in its position
  !---------------------------------------------------------  
  subroutine rotate(ksite,kvertex,kspin,kmolecule,rings3,rings4,tot_en,AR)
    use variables
    implicit none
    
    integer,intent(in) :: ksite,kvertex,kspin,kmolecule
    integer,intent(inout) :: rings3,rings4
    real,intent(inout) :: tot_en
    real,intent(out) :: AR
    integer :: kvertex_new
    integer :: rings3_old,rings3_new,rings4_old,rings4_new
    real :: old_en,new_en,delta_en,rand
    logical :: can_rotate
    
    if(kspin.eq.spinSN)then
       kvertex_new = int(1 + ran3(seed) * real(2))
    elseif(kspin.eq.spinSI)then
       kvertex_new = int(1 + ran3(seed) * real(nc))
    else
       print*,'BUG :: molecule not SI/SN :: subroutine rotate'
       print*,'kspin:',kspin
       print*,'kvertex:',kvertex
       print*,'kmolecule:',kmolecule
       print*,'ksite:',ksite
       call freeandclose
       stop
    end if
    
    !calculate the initial energy associated with the molecule ksite and then remove it
    call site_energy(ksite,kvertex,kspin,rings3_old,rings4_old,old_en)
    call remove_site(ksite,kvertex,kspin,kmolecule)
    
    !check if the molecule can be inserted again at ksite after rotation
    !rotation is defined as the assignment of the new pointer variable
    if(kspin.eq.spinSN)then
       call can_insert_SN(ksite,kvertex_new,can_rotate)
    elseif(kspin.eq.spinSI)then
       call can_insert_SI(ksite,kvertex_new,can_rotate)
    end if
    
    if(can_rotate)then
       
       !rotate the molecule at ksite and calculate the new energy associated with the molecule
       call insert_site(ksite,kvertex_new,kspin,kmolecule)
       call site_energy(ksite,kvertex_new,kspin,rings3_new,rings4_new,new_en)
       
       !calculate the change in energy and apply the metropolis criteria
       delta_en = new_en - old_en
       AR = exp(-1.000d0*delta_en/tstar)
       rand = ran3(seed)
       if(AR.gt.rand)then!if the metropolis criteria is satisfied
          tot_en = tot_en + delta_en
          rings3 = rings3 + (rings3_new - rings3_old)
          rings4 = rings4 + (rings4_new - rings4_old)
          SN_rotations_success = SN_rotations_success + 1
          
       else!if the metropolis criteria is not satisfied
          !put the molecule back in its original orientation
          call remove_site(ksite,kvertex_new,kspin,kmolecule)
          call insert_site(ksite,kvertex,kspin,kmolecule)          
          return
       end if!endif: the metropolis criteria is satisfied
       
    else
       !put the molecule back in its original orientation
       call insert_site(ksite,kvertex,kspin,kmolecule)
       return
    end if
    
    !safety check
    !call safetycheck_rotate(tot_en,kspin,kvertex,kvertex_new,AR,rand,old_en,new_en,rings3,rings4)
    
    return
  end subroutine rotate
  !---------------------------------------------------------
  !---print_output: prints the output into a file
  !--------------------------------------------------------- 
  subroutine print_output(AR1,AR2,AR3,rings3,rings4,en,comm_rank)
    use MPI
    use variables
    implicit none
    
    integer :: istep,points,comm_rank
    integer :: rings3,rings4,ierror
    real :: deltapoint,ratio,avg_particle_size
    real :: d0
    real :: AR1,AR2,AR3,en,norm_en,temp1,temp2,temp3,temp4
    integer,dimension(0:4) :: old_Qn,old_Qnn
    real,dimension(0:4) :: norm_Qn,new_Qn,norm_Qnn,new_Qnn
    
    !safety checking
    if(mod(step,2*nprint).eq.0)then
       call safetycheck(rings3,rings4,en)
    end if
    
    !write the output configuration
    if(mod(step,nprint).eq.0)then
       call config_out_atom
       !call config_out_molecule
       call output_state(comm_rank)
    end if

    !output Qn distribution
    !points for log plot of Qn distribution b/w 100 and 1000
    points = 8;ratio = real(1)/real(points)
    if(step.le.10)then
       call Qn_dist
       norm_Qn = real(Qn)/real(size(noccupy))
       norm_Qnn = real(Qnn)/real(ni(2))
       c = (norm_Qn(1) + norm_Qn(2)*2.0d0 + norm_Qn(3)*3.0d0 + norm_Qn(4)*4.0d0)/4.0d0
       cn = (norm_Qnn(1) + norm_Qnn(2)*2.0d0 + norm_Qnn(3)*3.0d0 + norm_Qnn(4)*4.0d0)/4.0d0
       write(1001,*) step,norm_Qn,c;flush(1001)
       write(1015,*) step,norm_Qnn,cn;flush(1015)
    else
       deltapoint = 10**(ratio)
       istep = int(deltapoint**qn_point)
       if((istep.eq.step).or.(mod(step,nprint).eq.0))then
          old_Qn = Qn
          old_Qnn = Qnn
          call Qn_dist
          new_Qn = real(Qn + old_Qn)/real(2)!average the Qn with the old value to get a smooth curve
          new_Qnn = real(Qnn + old_Qnn)/real(2)!average the Qnn with the old value to get a smooth curve
          norm_Qn = real(new_Qn)/real(size(noccupy))
          norm_Qnn = real(new_Qnn)/real(ni(2))
          c = (norm_Qn(1) + norm_Qn(2)*2.0d0 + norm_Qn(3)*3.0d0 + norm_Qn(4)*4.0d0)/4.0d0
          cn = (norm_Qnn(1) + norm_Qnn(2)*2.0d0 + norm_Qnn(3)*3.0d0 + norm_Qnn(4)*4.0d0)/4.0d0
          write(1001,*) step,norm_Qn,c;flush(1001)
          write(1015,*) step,norm_Qnn,cn;flush(1015)
          qn_point = qn_point + 1
       end if
    end if

    if(mod(step,nprint).eq.0)then

       !output cluster statistics only when all SN are present
       !call cluster_stat
       
       call take_snapshot(comm_rank)
       
       !find the negative charge per SI and sublattice ordering
       !call negchargeperSI       
       !call sublattice_ordering
       !call shell_order
              
       norm_en = en/real(nsites)
       write(*,8000) 'MC step:',step
       write(*,9000) 'Acceptance Ratios:',AR1,AR2,AR3
       write(*,10000) 'Energy per lattice site:',norm_en
       write(*,*)
    
       write(1008,8000) 'MC step:',step;flush(1008)
       write(1008,9000) 'Acceptance Ratios:',AR1,AR2,AR3;flush(1008)
       write(1008,10000) 'Energy per lattice site:',norm_en;flush(1008)
       write(1008,*);flush(1008)
       
       write(1002,*) step,tstar,rings3,rings4;flush(1002)
       write(1003,*) step,tstar,en,norm_en;flush(1003)
       write(1007,*) step,tstar,max_size,avg_size,cluster_num,d0;flush(1007)
       !write(1012,*) step,frac_SISDA,frac_SISNSDA;flush(1012)
    end if
    
    !calculate the average values of various quantities
    ! if((step.gt.neqsweeps).and.(mod(step,nprint).eq.0))then
    !    final_avg_energy = final_avg_energy + norm_en
    !    final_avg_max_size = final_avg_max_size + max_size
    !    final_avg_cluster_num = final_avg_cluster_num + cluster_num
    !    final_avg_avg_size = final_avg_avg_size + avg_size
       
    !    temp1 = final_avg_energy/real((step-neqsweeps)/nprint)
    !    temp2 = final_avg_max_size/real((step-neqsweeps)/nprint)
    !    temp3 = final_avg_cluster_num/real((step-neqsweeps)/nprint)
    !    temp4 = final_avg_avg_size/real((step-neqsweeps)/nprint)
       
    !    write(1014,12000) step,temp1,temp2,temp4,temp3;flush(1014)

    ! end if
    
    call cpu_time(end_time)
    runtime = (end_time-start_time)/3600.0d0
    !if((runtime.gt.time_limit).or.(step.ge.nsweeps))then
    if(step.ge.nsweeps)then

       call MPI_BARRIER(MPI_COMM_WORLD,ierror)       
       
       !output the final state of the syste
       call output_state(comm_rank)
       
       write(*,10000),'Program runtime(hrs):',runtime
       write(1008,10000)'Program runtime(hrs):',runtime
       write(*,*);write(*,*)
       
       !call freeandclose
       !stop

    end if
           
1000 format(x,A)
2000 format(x,A,x,i8,x,i8)
8000 format(x,A,x,i8)
9000 format(x,A,5x,es11.4,5x,es11.4,5x,es11.4)
10000 format(x,A,es12.4)
11000 format(x,A,5x,es10.4,5x,es10.4,5x,es10.4,5x,es10.4)
12000 format(x,i10,x,f10.4,x,f10.4,x,f10.4,x,f10.4)
13000 format(x,i10,x,f10.4)
    
    return
  end subroutine print_output
  !---------------------------------------------------------
  !---UPDATE_TEMPERATURE: updates a temperature according to a profile
  !--------------------------------------------------------- 
  subroutine update_temperature(step)
    implicit none
    
    integer,intent(in) :: step
    real :: frac

    frac = real(step)/real(nsweeps)
    
    if((frac.gt.m0).and.(frac.lt.m1))then
       tstar = tstar_initial + (frac - m0)*(up_rate)
    elseif((frac.gt.m1).and.(frac.lt.m2))then
       !do nothing
    elseif(frac.gt.m2)then
       tstar = tstar_initial + (tstar_final-tstar_initial)/exp(down_rate*real(frac-m2))
    end if
   
    write(1009,*)step,tstar
    
    return
  end subroutine update_temperature
  !---------------------------------------------------------
  !---cluster_stat(max_size):cluster counting subroutine
  !---------------------------------------------------------   
  subroutine cluster_stat
    use variables
    implicit none
    
    integer :: isite,ivertex,ispin,ihead
    integer :: p,q,monomer,neighs,scanned_neighs,cluster_size,same,molecule
    integer :: jsite,jlabel,pjlabel,jocc
    integer :: ksite,klabel,pklabel
    integer :: min_label,largest_label
    integer :: cluster_units
    logical :: arethereSi
    integer,dimension(4) :: mark
    !integer,dimension(:),allocatable :: neighSi
    integer,dimension(4) :: neighSi
    integer :: check_size,check_label
    
    csize = -1
    clabel = -1
    largest_label = 0
        
    do monomer = 1,allmolecules!loop over all monomers
       
       isite = noccupy(monomer)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if((occupancy(isite).eq.occSN).or.(occupancy(isite).eq.occSI))then!if isite has SN/SI

          neighSi=0
          call returnneigh(isite,neighSi,arethereSi)
          
          if(arethereSi)then!if there are neighboring SI/SN to isite

             neighs = size(neighSi)
             !check for scanned neighbors
             mark=-1;scanned_neighs = 0
             do p=1,neighs
                if(neighSi(p).gt.0)then
                   jsite = neighSi(p);jlabel = clabel(jsite)
                   if(jlabel.gt.0)then
                      scanned_neighs = scanned_neighs + 1
                      mark(p) = jsite
                   end if
                end if
             end do
             
             if(scanned_neighs.eq.1)then!if there's one scanned neighbor
                
                !if there only one scanned neighbor
                !assign isite the same label as the neighbor
                !and increment the size of the proper label of the neighbor by 1
                do p=1,neighs
                   if(mark(p).gt.0) jsite = mark(p)
                end do
                jlabel = clabel(jsite)
                call classify(jlabel,pjlabel)
                clabel(isite) = jlabel
                csize(pjlabel) = csize(pjlabel) + 1             
                
             elseif(scanned_neighs.ge.2)then!if there's more than one scanned neighbors
                !----------------------------------------------------------------------
                !find out the minimum proper label of the scanned neighbors
                min_label = allmolecules
                do p=1,neighs!loop over neighbors
                   if(mark(p).gt.0)then!if the site is scanned

                      jsite = mark(p);jlabel = clabel(jsite)
                      call classify(jlabel,pjlabel)
                      if(min_label.gt.pjlabel)then
                         min_label = pjlabel
                      end if

                   end if!endif: the site is scanned
                end do!enddo: loop over neighbors
                
                !coalesce the cluster into a larger one
                cluster_size = 0
                do p =1,neighs !loop over all neighbors
                   if(mark(p).gt.0)then!if the site is scanned
                      
                      jsite = mark(p);jlabel = clabel(jsite)
                      call classify(jlabel,pjlabel)
                      !check if pjlabel is a repetitive label amongst other labels
                      same = 0!variable to check if the label has already been operated upon
                      do q = 1,p-1
                         if(mark(q).gt.0)then
                            ksite = mark(q);klabel = clabel(ksite)
                            call classify(klabel,pklabel)
                            if(pklabel.eq.pjlabel) same = 1
                         end if
                      end do
                      
                      if((same.eq.0).and.(pjlabel.ne.min_label))then
                         cluster_size = cluster_size + csize(pjlabel)
                      end if
                   end if!endif: the site is scanned
                end do!enddo: loop over all neighbors
                
                !update the size of min_label
                csize(min_label) = csize(min_label) + cluster_size + 1
                
                !loop over the neighbors and set their sizes to zero except the min_label
                do p=1,neighs
                   if(mark(p).gt.0)then
                      jsite = mark(p);jlabel = clabel(jsite)
                      call classify(jlabel,pjlabel)
                      if(pjlabel.ne.min_label) csize(pjlabel) = -min_label
                   end if
                end do
                clabel(isite) = min_label!set the label of isite to the min_label
                !----------------------------------------------------------------------                
             else!if there's no scanned neighbor
                
                largest_label = largest_label + 1
                clabel(isite) = largest_label
                csize(largest_label) = 1
                
             end if!endif there's one scanned neighbor
             
          else!if there are no neighboring SI/SN to isite
             
             largest_label = largest_label + 1
             clabel(isite) = largest_label
             csize(largest_label) = 1
             
          end if!endif: there are neighboring SI/SN to isite
                    
       end if!endif: isite has SN/SI
       !if(allocated(neighSi)) deallocate(neighSi)

    end do!enddo: loop over all monomers
    
    !safety check
    check_size = 0
    do p=1,size(csize)
       if(csize(p).gt.0)check_size = check_size + csize(p)
    end do
    if(check_size.ne.(allmolecules-ni(3)))then
       print*,'BUG: cluster calculation error :: subroutine cluster_stat'
       print*,'check_size:',check_size
       print*,'monomer:',allmolecules
       call freeandclose
       stop       
    end if

    !loop over all the templates and label them according to the connected SI
    do molecule = 1,allmolecules
       
       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if((ispin.eq.spinSI).and.(occupancy(isite).eq.occSI))then
          ihead = n1list(ivertex,isite)
          do p = 1,6!loop through the second neigbors of O-
             jsite = n2list(p,ihead);jocc = occupancy(jsite)
             if(jocc.eq.occSDA)then
                clabel(jsite) = clabel(isite)
             end if
          end do!enddo: loop through the second neigbors of O-
          
       end if
    end do

    !calculate the various cluster statistics
    avg_size = 0;cluster_num = 0;cluster_units=0
    do p=1,size(csize)
       if(csize(p).ge.cluster_threshold)then
          cluster_num = cluster_num + 1
          cluster_units = cluster_units + csize(p)
       end if
    end do
    
    if(cluster_num.gt.0)then
       avg_size = real(cluster_units)/real(cluster_num)
       max_size = maxval(csize)
    else
       avg_size = real(0)
       max_size = int(0)
    end if
    
    return
  end subroutine cluster_stat
  !---------------------------------------------------------
  !---average_cluster_dimension: calculates average particle sizes
  !--------------------------------------------------------- 
  subroutine average_cluster_dimension(d0)
    use variables
    implicit none
    
    real,intent(out) :: d0
    real :: di,tot_di,tot_clus_units
    integer :: ilabel

    !call the cluster statistics
    call cluster_stat
    
    d0 = real(0)
    tot_di = real(0)
    tot_clus_units = real(0)
    
    if(max_size.ge.cluster_threshold)then!if max_size is greater than zero
       
       do ilabel = 1,nsites!loop through all the cluster labels
          
          if(csize(ilabel).ge.cluster_threshold)then
             call particle_size(di,ilabel)
             !tot_di = tot_di + di*csize(ilabel)!normalized by cluster units
             tot_di = tot_di + di
             tot_clus_units = tot_clus_units + csize(ilabel)!normalized by equal weights
          end if
          
       end do!loop through all the cluster labels
       
    end if!endif: max_size is greater than zero
    
    if(max_size.ge.cluster_threshold)then
       d0 = tot_di/(real(cluster_num))!normalized by equal weights
       !d0 = tot_di/(real(tot_clus_units))!normalized by cluster units
    else
       d0 = real(0)
    end if
    
    return
  end subroutine average_cluster_dimension
  !---------------------------------------------------------
  !---particle_size(d0) :: calculates the spatial dimensions of the particles
  !--------------------------------------------------------- 
  subroutine particle_size(d0,label)
    use variables
    implicit none
    
    real,intent(out) :: d0
    integer,intent(in) :: label
    real, parameter :: d=real(0.16),a0 = d/sqrt(real(3))
    integer :: molecule,jsite,jlabel,label_size,pjlabel,jocc
    integer :: xcm,ycm,zcm,cmsite,cmtag
    real :: tempd0,jcmdistance
    
    call cluster_stat
        
    label_size = csize(label)
    xcm = 0;ycm = 0;zcm = 0
    d0 = real(0);tempd0 = real(0)
    !first find the center of mass of the cluster with label j
    do molecule = 1,allmolecules
       jsite = noccupy(molecule);jocc = occupancy(jsite);jlabel = clabel(jsite)
       if((jlabel.gt.0).and.((jocc.eq.occSN).or.(jocc.eq.occSI)))then
       
          call classify(jlabel,pjlabel)
          if((pjlabel.eq.label).and.(csize(pjlabel).eq.label_size))then
             xcm = xcm + rx(jsite)
             ycm = ycm + ry(jsite)
             zcm = zcm + rz(jsite)
          end if

       end if
    end do

    xcm = int(real(xcm)/real(label_size))
    ycm = int(real(ycm)/real(label_size))
    zcm = int(real(zcm)/real(label_size))
    
    !calculate the closest site to CM on the lattice
    !cmcheck=mod(xcm,2)+mod(ycm,2)+mod(zcm,2)
    if(mod(xcm,2).eq.0)then
       if(mod(ycm,2).eq.0)then
          !do nothing
       elseif(mod(ycm,2).eq.1)then
          ycm = ycm + 1
       end if
       if(mod(zcm,2).eq.0)then
          !do nothing
       elseif(mod(zcm,2).eq.1)then
          zcm = zcm + 1
       end if
    elseif(mod(xcm,2).eq.1)then
       if(mod(ycm,2).eq.0)then
          ycm = ycm + 1
       elseif(mod(ycm,2).eq.1)then
          !do nothing
       end if
       if(mod(zcm,2).eq.0)then
          zcm = zcm + 1
       elseif(mod(zcm,2).eq.1)then
          !do nothing
       end if
    else
       print*,'BUG: lattice defect :: subroutine particle_size'
       call freeandclose
       stop
    end if
    
    cmtag = (xcm-1)*2*ly*2*lz + (ycm - 1)*2*lz + zcm
    cmsite = tag_site(cmtag)
    
    !now calculate the distance of every site with cm and calculate the
    !radius of gyration
    do molecule = 1,allmolecules
       jsite = noccupy(molecule);jocc = occupancy(jsite);jlabel = clabel(jsite)
       if((jlabel.gt.0).and.((jocc.eq.occSN).or.(jocc.eq.occSI)))then
          
          call classify(jlabel,pjlabel)
          if((pjlabel.eq.label).and.(csize(pjlabel).eq.label_size))then
             call distance(jsite,cmsite,jcmdistance)
             tempd0 = tempd0 + jcmdistance**2
          end if
          
       end if
    end do
    
    d0 = sqrt(tempd0/real(label_size))
    d0 = real(2)*d0*a0
    
    return
  end subroutine particle_size
  !---------------------------------------------------------
  !---CLASSIFY(JLABEL): CALCULATES THE PROPER LABEL OF A SITE
  !--------------------------------------------------------- 
  subroutine classify(jlabel,pjlabel)
    use variables
    implicit none
    
    integer,intent(in) :: jlabel
    integer,intent(out) :: pjlabel
    integer :: t,r
    
    r = jlabel
    t = r
    t = -csize(t)
    
    if(t.lt.0)then
       pjlabel = r
       return
    else
       r=t
       t = -csize(t)
       if(t.lt.0)then
          pjlabel = r
          return
       else
          
          do while(t.gt.0)
             r = t
             t = -csize(t)
          end do

          pjlabel = r
          csize(jlabel) = -r
          return

       end if
    end if
    
    return
  end subroutine classify
  !---------------------------------------------------------
  !---negchargeperSI: calculates the charges SI in solution
  !--------------------------------------------------------- 
  subroutine negchargeperSI
    use variables
    implicit none
    
    integer :: isite,ivertex,ispin,molecule
    integer :: i,occfn,bond
    integer,dimension(4) :: fn,fnsite
    logical :: dissolved_SI
    
    chargeperSI = 0
    
    do molecule = 1,allmolecules!loop over all molecules
       
       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if(ispin.eq.spinSI)then!if the molecule is SDA
          
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)             
          end do
          
          !check if SI is free or bonded by SN or SDA in the solid phase
          dissolved_SI = .true.
          do bond = 1,4
             occfn = occupancy(fnsite(bond))
             if((occfn.eq.occISIO).or.(occfn.eq.occSIO))then
                !do nothing
             elseif((occfn.eq.occSIOSNO).or.(occfn.eq.occSIOSIO))then
                dissolved_SI = .false.
             end if
          end do
          
          if(dissolved_SI)then
             !do nothing
          else
             chargeperSI = chargeperSI + 1.00
          end if
          
       end if!endif: the molecule is SDA
              
    end do!enddo:loop over all molecules
    
    chargeperSI = chargeperSI/real(ni(1))
    
    return
  end subroutine negchargeperSI
  !---------------------------------------------------------
  !---sublattice_ordering: calculates the sublattice ordering
  !                        in the system
  !--------------------------------------------------------- 
  subroutine sublattice_ordering
    use variables
    implicit none
    
    integer :: molecule,isite,ivertex,ispin
    integer :: i,occfn,bond
    integer,dimension(4) :: fn,fnsite
    integer :: different_connection

    sublat_ordering = 0.00
    
    do molecule = 1,allmolecules!loop over all molecules
       
       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if(ispin.eq.spinSI)then!if the molecule is SI
          
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)             
          end do
          
          !check the no. of SI connected to SN
          different_connection = 0
          do bond = 1,4
             occfn = occupancy(fnsite(bond))
             if(occfn.eq.occSIOSNO) different_connection = different_connection + 1
          end do
          
          if(different_connection.eq.3)then
             sublat_ordering = sublat_ordering + 1.00
          elseif(different_connection.eq.4)then
             print*,'BUG: occupancy of OH and O- :: subroutine sublattice_ordering'
             call freeandclose
             stop
          end if
          
       end if!endif: the molecule is SI
       
    end do!enddo: loop over all molecules
    
    sublat_ordering = sublat_ordering/real(ni(1))
    
    return
  end subroutine sublattice_ordering
  !---------------------------------------------------------
  !---shell_ordering: checks the fraction of SI connected to SI&SDA and only to SDA
  !--------------------------------------------------------- 
  subroutine shell_order
    use variables
    implicit none
    
    integer :: bond,molecule,isite,iocc,ivertex
    integer :: occfn
    integer,dimension(4) :: fn,fnsite
    logical :: connectedtoSN,connectedtoSDA
    
    frac_SISDA = real(0)
    frac_SISNSDA = real(0)
    
    do molecule = 1,allmolecules!loop over all the molecules
       isite = noccupy(molecule)
       iocc = occupancy(isite)
       ivertex = head(isite)
       
       if(iocc.eq.occSI)then!if isite has SI
          
          !find the neighbors
          do bond = 1,4
             fn(bond) = Si_O(bond,ivertex)
             fnsite(bond) = n1list(fn(bond),isite)
          end do
          
          !check for any connection with SN          
          connectedtoSN = .false.
          do bond = 2,4
             occfn = occupancy(fnsite(bond))
             if(occfn.eq.occSIOSNO) connectedtoSN = .true.
          end do

          !check for any connection with SDA          
          connectedtoSDA = .false.
          do bond = 1,6
             occfn = occupancy(n2list(bond,fnsite(1)))
             if(occfn.eq.occSDA) connectedtoSDA = .true.
          end do
          
          if(connectedtoSN.and.connectedtoSDA)then
             frac_SISNSDA = frac_SISNSDA + real(1)
          elseif(connectedtoSDA.and.(connectedtoSN.eqv..false.))then
             frac_SISDA = frac_SISDA + real(1)
          end if
          
       end if!endif: isite has SI       
       
    end do!enddo: loop over all the molecules

    frac_SISDA = frac_SISDA/real(ni(1))
    frac_SISNSDA = frac_SISNSDA/real(ni(1))
    
    return
  end subroutine shell_order
  !---------------------------------------------------------
  !---FAULTY: DO NOT USE
  !---config_out_molecule: output the configuration of molecule
  !--------------------------------------------------------- 
  subroutine config_out_molecule
    use variables
    implicit none
    
    integer :: isite, ivertex, ispin, ilabel, isize, pilabel
    integer :: monomer, bond, tot_atoms
    integer,dimension(4) :: fn
    integer :: i, j
    integer :: x, y, z, nx, ny, nz
    
    tot_atoms = ni(1)+ni(2) + ni(3)
    write(1021,*) tot_atoms;flush(1021)
    !write(1021,*) 'index, x, y, z';flush(1021)
    write(1021,*) 'STEP:',step;flush(1021)
    
    j=0
    do monomer=1,allmolecules
       
       isite = noccupy(monomer)
       ivertex = head(isite)
       ispin = spin(isite)
       
       ilabel = clabel(isite)
       if(ilabel.gt.0)then
          call classify(ilabel,pilabel)
          isize = csize(pilabel)

          if(isize.gt.cluster_threshold)then

             if(ispin.eq.spinSN)then
                !output SN configuration
                fn(1) = ivertex
                do bond = 2,4
                   fn(bond) = Si_O(bond,fn(1))
                end do
                x = rx(isite)
                y = ry(isite)
                z = rz(isite)
                j = j + 1
                write(1021,1000) 'Si',x,y,z;flush(1021)
                !do bond = 1,4
                !   j = j+1
                !   call output(fn(bond),x,y,z,nx,ny,nz)
                !   write(1021,1000) 'Si',nx,ny,nz;flush(1021)
                !end do
                
             elseif(ispin.eq.spinSI)then
                !output SI configuration
                fn(1) = ivertex
                do bond = 2,4
                   fn(bond) = Si_O(bond,fn(1))
                end do
                x = rx(isite)
                y = ry(isite)
                z = rz(isite)
                j = j + 1
                write(1021,1000) 'O',x,y,z;flush(1021)
                !do bond = 1,4
                !   j = j+1
                !   call output(fn(bond),x,y,z,nx,ny,nz)
                !   write(1021,1000) 'O',nx,ny,nz;flush(1021)
                !end do
             elseif(ispin.eq.spinSDA)then
                !output SDA configuration
                x = rx(isite)
                y = ry(isite)
                z = rz(isite)
                j = j + 1
                write(1021,1000) 'N',x,y,z;flush(1021)
                !do bond = 1,4
                !   j = j+1
                !   call output(fn(bond),x,y,z,nx,ny,nz)
                !   write(1021,1000) 'N',nx,ny,nz;flush(1021)
                !end do
             else
                print*,'Error: molecule other than SN/SI/SDA :: subroutine config_out_molecule'
                print*,'ispin',ispin
                call freeandclose
                stop          
             end if

          else

             if(ispin.eq.spinSN)then
                !output SN configuration
                fn(1) = ivertex
                do bond = 2,4
                   fn(bond) = Si_O(bond,fn(1))
                end do
                x = rx(isite)
                y = ry(isite)
                z = rz(isite)
                j = j + 1
                write(1021,1000) 'P',x,y,z;flush(1021)
                !do bond = 1,4
                !   j = j+1
                !   call output(fn(bond),x,y,z,nx,ny,nz)
                !   write(1021,1000) 'P',nx,ny,nz;flush(1021)
                !end do
                
             elseif(ispin.eq.spinSI)then
                !output SI configuration
                fn(1) = ivertex
                do bond = 2,4
                   fn(bond) = Si_O(bond,fn(1))
                end do
                x = rx(isite)
                y = ry(isite)
                z = rz(isite)
                j = j + 1
                write(1021,1000) 'P',x,y,z;flush(1021)
                !do bond = 1,4
                !   j = j+1
                !   call output(fn(bond),x,y,z,nx,ny,nz)
                !   write(1021,1000) 'P',nx,ny,nz;flush(1021)
                !end do
             elseif(ispin.eq.spinSDA)then
                !output SDA configuration
                x = rx(isite)
                y = ry(isite)
                z = rz(isite)
                j = j + 1
                write(1021,1000) 'P',x,y,z;flush(1021)
                !do bond = 1,4
                !   j = j+1
                !   call output(fn(bond),x,y,z,nx,ny,nz)
                !   write(1021,1000) 'P',nx,ny,nz;flush(1021)
                !end do
             else
                print*,'Error: molecule other than SN/SI/SDA :: subroutine config_out_molecule'
                print*,'ispin',ispin
                call freeandclose
                stop          
             end if
             
          end if
       else
          
          if(ispin.eq.spinSN)then
             !output SN configuration
             fn(1) = ivertex
             do bond = 2,4
                fn(bond) = Si_O(bond,fn(1))
             end do
             x = rx(isite)
             y = ry(isite)
             z = rz(isite)
             j = j + 1
             write(1021,1000) 'P',x,y,z;flush(1021)
             do bond = 1,4
                j = j+1
                call output(fn(bond),x,y,z,nx,ny,nz)
                write(1021,1000) 'P',nx,ny,nz;flush(1021)
             end do
             
          elseif(ispin.eq.spinSI)then
             !output SI configuration
             fn(1) = ivertex
             do bond = 2,4
                fn(bond) = Si_O(bond,fn(1))
             end do
             x = rx(isite)
             y = ry(isite)
             z = rz(isite)
             j = j + 1
             write(1021,1000) 'P',x,y,z;flush(1021)
             do bond = 1,4
                j = j+1
                call output(fn(bond),x,y,z,nx,ny,nz)
                write(1021,1000) 'P',nx,ny,nz;flush(1021)
             end do
          elseif(ispin.eq.spinSDA)then
             !output SDA configuration
             x = rx(isite)
             y = ry(isite)
             z = rz(isite)
             j = j + 1
             write(1021,1000) 'P',x,y,z;flush(1021)
             !do bond = 1,4
             !   j = j+1
             !   call output(fn(bond),x,y,z,nx,ny,nz)
             !   write(1021,1000) 'P',nx,ny,nz;flush(1021)
             !end do
          else
             print*,'Error: molecule other than SN/SI/SDA :: subroutine config_out_molecule'
             print*,'ispin',ispin
             call freeandclose
             stop          
          end if
          
       end if
       
    end do
    
    if(j.ne.tot_atoms)then
       print*,'Error: j.ne.allatoms :: subroutine config_out_molecule'
       print*,'j:',j
       print*,'allatoms:',tot_atoms
       call freeandclose
       stop
    end if
    
900 format(I8, I8, I8, I8) 
1000 format(A, I8, I8, I8)
    
    return
  end subroutine config_out_molecule
  !---------------------------------------------------------
  !---config_out_atom: output the configuration of atoms
  !--------------------------------------------------------- 
  subroutine config_out_atom
    use variables
    implicit none
    
    integer :: isite, ivertex, ispin
    integer :: monomer, bond
    integer,dimension(4) :: fn
    integer :: i, j
    integer :: x, y, z, nx, ny, nz
    
    write(1006,*) natom;flush(1006)
    write(1006,*) 'index, x, y, z';flush(1006)
    
    j=0
    do monomer=1,allmolecules
       
       isite = noccupy(monomer)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if(ispin.eq.spinSN)then
          !output SN configuration
          fn(1) = ivertex
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
          end do
          x = rx(isite)
          y = ry(isite)
          z = rz(isite)
          j = j + 1
          write(1006,900) j,x,y,z;flush(1006)
          do bond = 1,4
             j = j+1
             call output(fn(bond),x,y,z,nx,ny,nz)
             write(1006,900) j,nx,ny,nz;flush(1006)
          end do

       elseif(ispin.eq.spinSI)then
          !output SI configuration
          fn(1) = ivertex
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
          end do
          x = rx(isite)
          y = ry(isite)
          z = rz(isite)
          j = j + 1
          write(1006,900) j,x,y,z;flush(1006)
          do bond = 1,4
             j = j+1
             call output(fn(bond),x,y,z,nx,ny,nz)
             write(1006,900) j,nx,ny,nz;flush(1006)
          end do
       elseif(ispin.eq.spinSDA)then
          !output SDA configuration
          x = rx(isite)
          y = ry(isite)
          z = rz(isite)
          j = j + 1
          write(1006,900) j,x,y,z;flush(1006)
          !do bond = 1,4
          !   j = j+1
          !   call output(fn(bond),x,y,z,nx,ny,nz)
          !   write(1006,900) j,nx,ny,nz;flush(1006)
          !end do
       else
          print*,'Error: molecule other than SN/SI/SDA :: subroutine config_out_atom'
          print*,'ispin',ispin
          call freeandclose
          stop          
       end if
       
    end do
    
    if(j.ne.natom)then
       print*,'Error: j.ne.allatoms :: subroutine config_out_atom'
       print*,'j:',j
       print*,'allatoms:',natom
       call freeandclose
       stop
    end if
    
900 format(I8, I8, I8, I8)     
    
    return
  end subroutine config_out_atom
  !---------------------------------------------------------
  !---take_snapshot(snapshot)
  !---------------------------------------------------------   
  subroutine take_snapshot(comm_rank)
    use variables
    implicit none
    
    integer,intent(in) :: comm_rank
    integer :: isite, ivertex, ispin
    integer :: monomer, bond, molecule
    integer,dimension(4) :: fn
    character(len=16) :: filenum,prefix,extension,temp1,filename
    character(len=200) :: pwd,proc_num,ifile

    snapshot = snapshot + 1
    
    write(proc_num,'(I2.2)')comm_rank
    write(filenum,'(I4.4)')snapshot    
    call getcwd(pwd)
    
    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('snap')//'_'//trim(filenum)//'.'//trim('cfg')
    open(1016,file=ifile,status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening snapshot file'
       call freeandclose
       stop
    end if
            
    write(1016,*)step,snapshot,qn_point,tstar
    
    !output the configuration of the molecules
    do molecule = 1,allmolecules
       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       write(1016,*)molecule,isite,ivertex,ispin
    end do

200 format(i8,x,i8,x,i8)
300 format(f6.2)
400 format(f10.2,x,f10.2,x,f10.2)
500 format(i8,x,i8,x,i8,x,i8)
600 format(i8,x,i8,x,i8,x,i8)
700 format(f6.2,x,f6.2)
800 format(f6.4,x,f6.4,x,f6.4,x,f6.4)
900 format(i15)    
    
    close(1016)
    return
  end subroutine take_snapshot
  !---------------------------------------------------------
  !---output_state: outputs the state of the system in the file system_state.cfg
  !--------------------------------------------------------- 
  subroutine output_state(comm_rank)
    use variables
    implicit none
    
    integer,intent(in) :: comm_rank
    integer :: opstatus, molecule, isite, ispin ,ivertex
    character(len=200) :: pwd,proc_num,ifile
    
    write(proc_num,'(I2.2)')comm_rank
    call getcwd(pwd)
    
    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('system_state.cfg')
    open(1017,file=ifile,status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening system_state.cfg file'
       call freeandclose
       stop
    end if
    
    write(1017,*)step,snapshot,qn_point,tstar
    
    !output the configuration of the molecules
    do molecule = 1,allmolecules
       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       write(1017,*)molecule,isite,ivertex,ispin
    end do

200 format(i8,x,i8,x,i8)
300 format(f6.2)
400 format(f10.2,x,f10.2,x,f10.2)
500 format(i8,x,i8,x,i8,x,i8)
600 format(i8,x,i8,x,i8,x,i8)
700 format(f6.2,x,f6.2)
800 format(f6.4,x,f6.4,x,f6.4,x,f6.4)
900 format(i15)
    
    close(1017)
    return
  end subroutine output_state
  !---------------------------------------------------------
  !---read_state: outputs the state of the system in the file system_state.cfg
  !--------------------------------------------------------- 
  subroutine read_state(comm_rank)
    use variables
    use input_SiO4_SDA
    implicit none
    
    integer,intent(in) :: comm_rank
    integer :: opstatus, data
    integer :: molecule, isite, ispin ,ivertex
    integer,dimension(4) :: fn,fnsite
    logical :: can_insert
    character(len=200) :: pwd,proc_num,ifile
    
    occupancy = occW
    spin = spinW
    head=0
    noccupy = 0
    resite = 0
    
    write(proc_num,'(I2.2)')comm_rank
    call getcwd(pwd)
    ifile = trim(pwd)//'/'//trim('temp')//'_'//trim(proc_num)//'/'//trim('system_state.cfg')
    open(1017,file=ifile,status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening system_state.cfg file :: subroutine read_state'
       call freeandclose
       stop
    end if
        
    read(1017,*)initial_step,snapshot,qn_point,tstar
    
    !output the configuration of the molecules
    do data = 1,allmolecules
       
       read(1017,*)molecule,isite,ivertex,ispin

       if(data.ne.molecule)then
          print*,'BUG: data.ne.molecule :: subroutine read_output'
          print*,'data/molecule',data/molecule
          call freeandclose
          stop
       end if
              
       if(ispin.eq.spinSN)then!if the molecule is SN
          call can_insert_SN(isite,ivertex,can_insert)
          if(can_insert)then
             call insert_site(isite,ivertex,ispin,molecule)
          else
             print*,'BUG: Reading SN from system_state.cfg :: subroutine read_state'
             call freeandclose
             stop
          end if
       elseif(ispin.eq.spinSI)then!if the molecule is SI
          call can_insert_SI(isite,ivertex,can_insert)
          if(can_insert)then
             call insert_site(isite,ivertex,ispin,molecule)
          else
             print*,'BUG: Reading SI from system_state.cfg :: subroutine read_state'
             call freeandclose
             stop
          end if
       elseif(ispin.eq.spinSDA)then!if the molecule is SDA
          call can_insert_SDA(isite,can_insert)
          if(can_insert)then
             call insert_site(isite,ivertex,ispin,molecule)
          else
             print*,'BUG: Reading SDA from system_state.cfg :: subroutine read_state'
             call freeandclose
             stop
          end if
       end if!endif: the molecule is SN
       
    end do
    
    !output the psf file
    call psf_out(comm_rank)

    !output initial configuration
    call config_out_atom

200 format(i8,x,i8,x,i8)
300 format(f6.2)
400 format(f10.2,x,f10.2,x,f10.2)
500 format(i8,x,i8,x,i8,x,i8)
600 format(i8,x,i8,x,i8,x,i8)
700 format(f6.2,x,f6.2)
800 format(f6.4,x,f6.4,x,f6.4,x,f6.4)
900 format(i15)

    close(1017)    
    return
  end subroutine read_state
  !---------------------------------------------------------
  !---output(ivertex,x,y,z,nx,ny,nz)
  !--------------------------------------------------------- 
  subroutine output(ivertex, x, y, z, nx, ny, nz)
    implicit none
    
    integer :: ivertex,x,y,z,nx,ny,nz
    
    if(ivertex.eq.1)then
       nx = x - 1
       ny = y - 1
       nz = z + 1
    elseif(ivertex.eq.2)then
       nx = x + 1
       ny = y - 1
       nz = z + 1
    elseif(ivertex.eq.3)then
       nx = x + 1
       ny = y + 1
       nz = z + 1
    elseif(ivertex.eq.4)then
       nx = x - 1
       ny = y + 1
       nz = z + 1
    elseif(ivertex.eq.5)then
       nx = x - 1
       ny = y - 1
       nz = z - 1
    elseif(ivertex.eq.6)then
       nx = x + 1
       ny = y - 1
       nz = z - 1
    elseif(ivertex.eq.7)then
       nx = x + 1
       ny = y + 1
       nz = z - 1
    elseif(ivertex.eq.8)then
       nx = x - 1
       ny = y + 1
       nz = z - 1
    else
       print*,'Error:: ivertex must be between 1 and 8:: subroutine output'
       call freeandclose
       stop
    end if
    
    return
  end subroutine output
  !---------------------------------------------------------
  !---LCG RANDOM NUMBER GENERATOR
  !--------------------------------------------------------- 
  function ran3(idum)
    implicit none
    integer idum
    integer mbig,mseed,mz
    !real*8 mbig,mseed,mz
    real*8 ran3,fac
    parameter(mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)  
    !parameter (mbig=4000000.,mseed=1618033.,mz=0.0,fac=1./mbig)
    integer i,iff,ii,inext,inextp,k
    integer mj,mk,ma(55)
    !real*8 mj,mk,ma(55)
    save iff,inext,inextp,ma
    data iff  /0/
    if(idum.lt.0.or.iff.eq.0) then
       iff=1
       mj=abs(mseed-abs(idum))
       mj=mod(mj,mbig)
       ma(55)=mj
       mk=1
       do i = 1,54
          ii=mod(21*i,55)
          ma(ii) = mk
          mk = mj-mk
          if(mk.lt.mz) mk = mk + mbig
          mj = ma(ii)
       enddo
       do  k=1,4
          do  i = 1,55
             ma(i) = ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.mz)ma(i) = ma(i) + mbig
          enddo
       enddo
       inext = 0
       inextp = 31
       idum = 1
    endif
    inext = inext + 1
    if(inext.eq.56) inext = 1
    inextp = inextp + 1
    if(inextp.eq.56) inextp = 1
    mj = ma(inext) - ma(inextp)
    if(mj.lt.mz)mj = mj + mbig
    ma(inext) = mj
    ran3 = mj*fac
    return
  end function ran3
  !---------------------------------------------------------
  !--- FREE_ALL: DEALLOCATE ALL MEMORY
  !---------------------------------------------------------
  subroutine freeandclose
    use variables
    use MPI
    implicit none
    
    !deallocate all arrays
    deallocate(occupancy)
    deallocate(tag_site)
    
    deallocate(noccupy)
    deallocate(head)
    deallocate(spin)
    deallocate(resite)
    deallocate(clabel)
    deallocate(csize)
    
    deallocate(n1list)
    deallocate(n2list)
    deallocate(n3list)
    deallocate(n4list)
    deallocate(n5list)

    deallocate(rx)
    deallocate(ry)
    deallocate(rz)
    
    deallocate(proc_temp)
    
    !close all opened files
    close(1001);close(1002);close(1003);close(1004)
    close(1006);close(1007);close(1008)
    !close(1009);close(1010);close(1011)
    close(1012);close(1013);close(1014);close(1015)
    close(1018);close(1019);close(1021);close(5000)
    
    !call MPI_FINALIZE

    return
  end subroutine freeandclose
  !***************************************************************************
  !***SAEFTY CHECKING
  !***************************************************************************
  subroutine safetycheck(rings3,rings4,tot_en)
    use variables
    implicit none
    
    integer,intent(in) :: rings3,rings4
    real,intent(in) :: tot_en
    integer :: rings3_temp=0,rings4_temp=0,i,occen,occfn
    integer :: j,molecule,bond,isite,ispin,ivertex,ihead,iocc
    real :: temp_en,ratio,old_en
    real*8 :: epsilon=1E-3
    integer,dimension(4) :: fn,fnsite
    
    do i=1,nsites
       occen = occupancy(i)
       if((occen.eq.occW).or.(occen.eq.occSN).or.(occen.eq.occSI).or.(occen.eq.occSNOSNO)&
            .or.(occen.eq.occISIO).or.(occen.eq.occSIOSNO).or.(occen.eq.occSNO).or.&
            (occen.eq.occSIO).or.(occen.eq.occSDA))then
       else
          print*,'BUG : site has invalid occupancy :: subroutine safetycheck'
          print*,i,rx(i),ry(i),rz(i),occen
          call freeandclose
          stop
       end if
    end do
    
    call total_energy(temp_en)
    call rings3and4(rings3_temp,rings4_temp)
    ratio = abs(1-(temp_en/tot_en))
    !ratio = abs(tot_en-temp_en)
    if((ratio.gt.epsilon).or.(rings3_temp.ne.rings3).or.(rings4_temp.ne.rings4))then
       print*,'BUG :: energy/rings3/rings4 calculation :: subroutine safetycheck'
       print*,'Energy from subroutine:',temp_en
       print*,'Energy from simulation:',tot_en
       print*,'ratio/epsilon:',ratio,epsilon
       print*,'rings3/rings4 from subroutine:',rings3_temp,rings4_temp
       print*,'rings3/rings4 from simulation:',rings3,rings4
       print*,
       call freeandclose
       stop
    end if
    
    do molecule=1,allmolecules!loop over all molecules

       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       iocc = occupancy(isite)
       
       if(iocc.eq.occSN)then!if molecule is SN
          
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)
          end do
          
          do bond = 1,8
             occfn = occupancy(n1list(bond,isite))
             if((occfn.eq.occW).or.(occfn.eq.occISIO).or.(occfn.eq.occSIO).or.&
                  (occfn.eq.occSNO).or.(occfn.eq.occSNOSNO).or.&
                  (occfn.eq.occSI).or.(occfn.eq.occSN).or.(occfn.eq.occSIOSNO))then
             else
                print*,'BUG:invalid occupancy in neighborhood of SI/SN :: subroutine safeycheck' 
                print*,'identity:',ispin
                print*,'occfn',bond,occfn
                call freeandclose
                stop
             end if
          end do
          
          do bond = 1,4
             ihead = fnsite(bond)
             do j=1,8
                occfn = occupancy(n1list(j,ihead))
                if(occfn.eq.occSDA)then
                   print*,'BUG:incorrect SDA occupancy :: subroutine safetycheck'
                   print*,'occfn',occfn
                   call freeandclose
                   stop
                end if
             end do
          end do
          
       elseif(iocc.eq.occSI)then!if molecule is SI

          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)
          end do
          
          do bond = 1,8
             occfn = occupancy(n1list(bond,isite))
             if((occfn.eq.occW).or.(occfn.eq.occISIO).or.(occfn.eq.occSIO).or.&
                  (occfn.eq.occSNO).or.(occfn.eq.occSNOSNO).or.&
                  (occfn.eq.occSI).or.(occfn.eq.occSN).or.(occfn.eq.occSIOSNO))then
             else
                print*,'BUG:invalid occupancy in neighborhood of SI/SN :: subroutine safeycheck' 
                print*,'identity:',ispin
                print*,'occfn',bond,occfn
                call freeandclose
                stop
             end if
          end do
          
          do bond = 1,4
             ihead = fnsite(bond)
             do j=1,8
                occfn = occupancy(n1list(j,ihead))
                if(occfn.eq.occSDA)then
                   print*,'BUG:incorrect SDA occupancy :: subroutine safetycheck'
                   print*,'occfn',occfn
                   call freeandclose
                   stop
                end if
             end do
          end do
          
          ihead = fnsite(1)
          !check second neighbors of O-
          !do j=1,6
          !   occfn = occupancy(n2list(j,ihead))
          !   if(occfn.eq.occISIO)then
          !      print*,'BUG: second neighbor of O- is O-:: subroutine safetycheck'
          !      print*,'occfn/occISIO:',occfn,occISIO
          !      call freeandclose
          !      stop
          !   end if
          !end do          
          !check third neighbors of O-
          do j=1,12
             occfn = occupancy(n3list(j,ihead))
             if(occfn.eq.occISIO)then
                print*,'BUG: third neighbor of O- is O-:: subroutine safetycheck'
                print*,'occfn/occISIO:',occfn,occISIO
                call freeandclose
                stop
             end if
          end do
          
       elseif(iocc.eq.occSDA)then!if molecule is SDA
          !check first neighbors of template
          do j=1,8
             occfn = occupancy(n1list(j,isite))
             if(occfn.eq.occW)then
             else
                print*,'BUG: first neighbor of SDA not water:: subroutine safetycheck'
                call freeandclose
                stop
             end if
          end do
          !check second neighbors of template
          !do j=1,6
          !   occfn = occupancy(n2list(j,isite))
          !   if(occfn.eq.occSDA)then
          !      print*,'BUG: second neighbor of SDA is SDA:: subroutine safetycheck'
          !      print*,'occfn/occSDA:',occfn,occSDA
          !      call freeandclose
          !      stop
          !   end if
          !end do          
          !check third neighbors of template
          !do j=1,12
          !   occfn = occupancy(n3list(j,isite))
          !   if(occfn.eq.occSDA)then
          !      print*,'BUG: third neighbor of SDA is SDA:: subroutine safetycheck'
          !      print*,'occfn/occISIO:',occfn,occISIO
          !      call freeandclose
          !      stop
          !   end if
          !end do
       else
          print*,'BUG: molecule neither SI/SN/SDA: subroutine safetycheck'
          print*,'iocc:',iocc
          call freeandclose
          stop
       end if
       
    end do
    
  end subroutine safetycheck
  !---------------------------------------------------------------------------
  subroutine safetycheck_translate(tot_en,ispin,AR1,rand,old_en,new_en,rings3,rings4)
    use variables
    implicit none
    
    integer :: ispin,rings3,rings4,rings3_temp,rings4_temp
    real :: tot_en,old_en,new_en,temp_en,delta_en,AR1,rand,ratio
    real*8 :: epsilon = 0.00000001d0
    
    call total_energy(temp_en)
    call rings3and4(rings3_temp,rings4_temp)
    ratio = abs(1-(temp_en/tot_en))
    if((ratio.gt.epsilon).or.(rings3_temp.ne.rings3).or.(rings4_temp.ne.rings4))then
       print*,'BUG :: energy calculation :: subroutine safetycheck_translate'
       print*,'ispin:',ispin
       print*,'AR1:',AR1
       print*,'ran3(seed)',rand
       print*,'Energy from subroutine:',temp_en
       print*,'Energy from simulation:',tot_en
       print*,'old_en',old_en
       print*,'new_en:',new_en
       print*,'delta_en',new_en-old_en
       print*,'rings3/rings4 from subroutine:',rings3_temp,rings4_temp
       print*,'rings3/rings4 from simulation:',rings3,rings4
       print*,
       call freeandclose
       stop
    end if
    
    return
  end subroutine safetycheck_translate
  !---------------------------------------------------------------------------
  subroutine safetycheck_swap(tot_en,isite,jsite,AR1,rand,old_en,new_en,rings3,rings4)
    use variables
    implicit none
    
    integer :: isite,rings3,rings4,rings3_temp,rings4_temp,jsite
    real :: tot_en,old_en,new_en,temp_en,delta_en,AR1,rand,ratio
    real*8 :: epsilon = 0.00000001d0
    
    integer :: bond,ihead
    
    call total_energy(temp_en)
    call rings3and4(rings3_temp,rings4_temp)
    ratio = abs(1-(temp_en/tot_en))
    if((ratio.gt.epsilon).or.(rings3_temp.ne.rings3).or.(rings4_temp.ne.rings4))then
       print*,'BUG :: energy calculation :: subroutine safetycheck_swap'
       print*,'ispin/jspin:',spin(isite),spin(jsite)
       print*,'AR1:',AR1
       print*,'ran3(seed):',rand
       print*,'Energy from subroutine:',temp_en
       print*,'Energy from simulation:',tot_en
       print*,'old_en',old_en
       print*,'new_en',new_en
       print*,'delta_en',new_en-old_en
       print*,'rings3/rings4 from subroutine:',rings3_temp,rings4_temp
       print*,'rings3/rings4 from simulation:',rings3,rings4
       print*,
       
       !neighbor information
       print*,'---imolecule information---'
       print*,'ivertex/jocc/spin',head(isite),occupancy(isite),spin(isite)
       print*,'rx(isite),ry(isite),rz(isite)',rx(isite),ry(isite),rz(isite)
       print*,'first neighbors'
       do bond = 1,8
          print*,bond,occupancy(n1list(bond,isite)),rx(n1list(bond,isite)),&
               ry(n1list(bond,isite)),rz(n1list(bond,isite))
       end do
       print*,'second neighbors'
       do bond = 1,6
          print*,bond,occupancy(n2list(bond,isite)),rx(n2list(bond,isite)),&
               ry(n2list(bond,isite)),rz(n2list(bond,isite))
       end do
       
       !neighbor information
       print*,'---jmolecule information---'
       print*,'jvertex/jocc/jspin',head(jsite),occupancy(jsite),spin(jsite)
       print*,'rx(jsite),ry(jsite),rz(jsite)',rx(jsite),ry(jsite),rz(jsite)
       print*,'first neighbors'
       do bond = 1,8
          print*,bond,occupancy(n1list(bond,jsite)),rx(n1list(bond,jsite)),&
               ry(n1list(bond,jsite)),rz(n1list(bond,jsite))
       end do
       print*,'second neighbors'
       do bond = 1,6
          print*,bond,occupancy(n2list(bond,jsite)),rx(n2list(bond,jsite)),&
               ry(n2list(bond,jsite)),rz(n2list(bond,jsite))
       end do
       
       call freeandclose
       stop
    end if
    
    return
  end subroutine safetycheck_swap
  !---------------------------------------------------------------------------
  subroutine safetycheck_rotate(tot_en,kspin,kvertex,kvertex_new,AR2,rand,old_en,new_en,rings3,rings4)
    use variables
    implicit none
    
    integer :: kspin,rings3,rings4,rings3_temp,rings4_temp,kvertex,kvertex_new
    real :: tot_en,old_en,new_en,temp_en,delta_en,AR2,rand,ratio
    real*8 :: epsilon = 0.00000001d0
    
    call total_energy(temp_en)
    call rings3and4(rings3_temp,rings4_temp)
    ratio = abs(1-(temp_en/tot_en))
    if((ratio.gt.epsilon).or.(rings3_temp.ne.rings3).or.(rings4_temp.ne.rings4))then
       print*,'BUG :: energy calculation :: subroutine safetycheck_rotate'
       print*,'kspin:',kspin
       print*,'AR2:',AR2
       print*,'ran3(seed)',rand
       print*,'Energy from subroutine:',temp_en
       print*,'Energy from simulation:',tot_en
       print*,'old_en',old_en
       print*,'new_en:',new_en
       print*,'delta_en',new_en-old_en
       print*,'rings3/rings4 from subroutine:',rings3_temp,rings4_temp
       print*,'rings3/rings4 from simulation:',rings3,rings4
       print*,
       call freeandclose
       stop
    end if
    
    return
  end subroutine safetycheck_rotate
  !---------------------------------------------------------
  !---manual_check: CHECK FOR CORRECTNESS MANUALLY
  !--------------------------------------------------------- 
  subroutine manual_check
    use variables
    implicit none
    integer,dimension(13) :: s
    integer,dimension(18) :: o
    integer,dimension(17) :: d
    integer :: ivertex,ispin,monomer,tag
    logical :: arethereSi,can_insert
    integer :: isite,jsite,rings3,rings4,i,j,overlap
    integer,dimension(32) :: neighSi1
    integer,dimension(:),allocatable :: neighSi
    real :: en,en_temp
    !
    s(1) = tag_site(73326)
    s(2) = tag_site(73806)
    s(3) = tag_site(101886)
    s(4) = tag_site(102366)
    s(5) = tag_site(102846)
    s(6) = tag_site(130446)
    s(7) = tag_site(130926)
    s(8) = tag_site(131406)
    s(9) = tag_site(159006)
    s(10) = tag_site(159966)
    s(11) = tag_site(101648)
    s(12) = tag_site(158768)
    s(13) = tag_site(187566)
    !
    o(1) = tag_site(58807)
    o(2) = tag_site(59045)
    o(3) = tag_site(59287)
    o(4) = tag_site(59525)
    o(5) = tag_site(88565)
    o(6) = tag_site(117367)
    o(7) = tag_site(144967)
    o(8) = tag_site(145205)
    o(9) = tag_site(145927)
    o(10) = tag_site(173527)
    o(11) = tag_site(174487)
    o(12) = tag_site(174245)
    o(13) = tag_site(116169)
    o(14) = tag_site(87129)
    o(15) = tag_site(144249)
    o(16) = tag_site(173289)
    o(17) = tag_site(202087)
    o(18) = tag_site(201845)
    !
    d(1) = tag_site(88327)
    d(2) = tag_site(87847)
    d(3) = tag_site(87367)
    d(4) = tag_site(87605)
    d(5) = tag_site(88085)
    d(6) = tag_site(116887)
    d(7) = tag_site(116407)
    d(8) = tag_site(115927)
    d(9) = tag_site(116165)
    d(10) = tag_site(116645)
    d(11) = tag_site(117125)
    d(12) = tag_site(145447)
    d(13) = tag_site(144487)
    d(14) = tag_site(144725)
    d(15) = tag_site(145685)
    d(16) = tag_site(173047)
    d(17) = tag_site(173285)    
    !
    !set the occupancy and spin of sites to that of water: 
    !all sites are occupied by water
    !the sites of lattice initially have no pointer
    occupancy = occW
    spin = spinW
    head = 0
    !
    do i=1,size(s)
       occupancy(s(i)) = occSN
    end do
    do i=1,size(d)
       occupancy(d(i)) = occSNOSNO
    end do
    !
    do i=1,size(o)
       occupancy(o(i)) = occSNO
    end do
    !
    do i=1,size(s)
       noccupy(i) = s(i)
    end do
    !
    do i = 1,size(s)
       head(s(i)) = 1
    end do
    !
    do i=1,size(s)
       spin(s(i)) = spinSN
    end do
    
    !call psf_out(comm_rank)
    call config_out_atom
    
    do monomer=1,allmolecules
       isite = noccupy(monomer);ivertex=head(isite)
       call site_rings3and4(isite,head(isite),rings3,rings4)
       print*,monomer,rings3,rings4
    end do
    !    
    return    
  end subroutine manual_check
  !********************************************************************
end program PT_SiO4_TAA
