module variables
  implicit none
  
  !dimensions of the lattice
  integer :: lx,ly,lz
  real,parameter :: a = real(2)*real(0.16)/sqrt(real(3)) !in nanometers

  !tag_site: array for site number <--> site label transformation
  integer,dimension(:),allocatable :: tag_site

  !arrays storing the neighbor list:
  !n1list: first neighbor list
  !n2list: second neighbor list
  !n3list: third neighbor listweq
  integer,dimension(:,:),allocatable :: n1list,n2list,n3list,n4list,n5list
  integer,dimension(:),allocatable :: neighbors
  integer,dimension(:,:,:),allocatable :: nlist

  !nsites:total sites
  !natom: total atoms 
  !allmonomers: all the monomers on the lattice
  integer :: nsites,natom,allmolecules,unitcells

  ! variables for SDA size
  integer :: sda_size
  real :: sda_dia

  !initial configuration flag,restart code flag, MC step, initial MC step
  integer :: initial_config_flag,step,initial_step
  !forerunnung string for configuration files
  character*10 :: snapshot_
  
  !coordinates of the particles 
  integer,dimension(:),allocatable :: rx,ry,rz
  
  !arrays storing the occupancy of sites
  integer,dimension(:),allocatable :: occupancy

  !occupancies parameters
  integer,parameter :: occSI = 7    !occupancy of Si in SI
  integer,parameter :: occSN = 2    !occupancy of Si in SN
  integer,parameter :: occW = 0     !occupancy of water(W)
  integer,parameter :: occSDA = 13   !occupancy of SDA site
  integer,parameter :: occISIO = 113 !occupancy of oxygen in SI
  integer,parameter :: occSNO = 17   !occupancy of oxygen in SN
  integer,parameter :: occSIO = 51   !occupancy of oxygen in SI
  !bridging O- between SI and SN
  integer,parameter :: occSIOSNO = occSIO + occSNO
  !bridging O- between SI and SI
  integer,parameter :: occSIOSIO = occSIO + occSIO
  !bridging O- between SN and SN
  integer,parameter :: occSNOSNO = occSNO + occSNO
  
  !array storing the identity of molecules on sites
  !1:SI,2:SN,3:SDA,4:Water(W)
  integer,dimension(:),allocatable :: spin
  integer,parameter :: spinSN = 2,spinSI = 1,spinW = 4,spinSDA = 3

  !noccupy:array storing the site of Si molecules
  !head:target of all the sites
  !resite: transforms isite to monomer no.
  integer,dimension(:),allocatable :: noccupy, head, resite

  !clabel: stores the labels of the clusters
  !csize: stores the size of the corresponding cluster
  integer,dimension(:),allocatable :: clabel,csize
  integer :: max_size,cluster_num,cluster_threshold
  real :: avg_size
  
  !array storing the orientation detail of the monomer unit
  integer,dimension(4,8) :: Si_O,Si_notO
  
  !variables for the no of sweeps,equilibrium sweeps
  !and snapshot taking variables
  integer :: nsweeps,nprint,qn_point,snapshot
  
  !temperature,total energy variable
  real :: tstar
  
  !variable to store the charge per SI
  real :: chargeperSI,sublat_ordering,frac_SISDA,frac_SISNSDA
  
  !Qn distribution variables
  integer,dimension(0:4) :: Qn,Qni,Qnn
  
  !rings size distribution variables
  !integer :: rings3,rings4
  
  integer :: seed!seed for the random number generator
  integer,parameter :: nc=8!coordination number
  real :: c,cn!degree of condensation
  
  !arrays for storing concentrations and the number of molecules
  !indices same as the spin number
  real,dimension(4) :: xi
  integer,dimension(4) :: ni
  real :: nTEOS,nSDA,nH2O,pH
  
  !penalties on 3 and 4 membered rings
  real :: pen3, pen4
  
  !interaction energy strengths
  !eISISDA: interaction energy b/w ISI and SDA
  !eSNSDA: interaction energy b/w SN and SDA
  !eSNSN:     "         "    "   SN and SN
  !eSISN      "         "    "   SI and SN
  !eSISI      "         "    "   SI and S
  !eSDASDA        "         "    "   SDA and SDA
  !real :: eISISDA,eSNSDA,eSNSN,eSISN,eSISI,eSDASDA
  real :: eOSDA, eSNSN, eSISN

  !variable for the temperature profile change
  real :: tstar_initial,tstar_final
  real :: up_rate,down_rate
  real :: m0,m1,m2,m3
  
  !final averages of various variables
  real :: final_avg_energy=0.00
  real :: final_avg_max_size=0.00,final_avg_cluster_num=0.00,final_avg_avg_size=0.00
  
  !time variables
  real :: start_time,end_time,runtime,time_limit
  
  !PARALLEL TEMPERING variables
  integer :: ntemp, iproc, jproc, nexchange, jump_flag
  real :: ftemp,dtemp
  real :: proci_temp,procj_temp,proci_en,procj_en
  real,dimension(:),allocatable :: proc_temp
  
  !move couters
  integer :: SN_jumps_attempt,SN_jumps_success
  integer :: SN_rotations_attempt, SN_rotations_success
  integer :: SDA_jumps_attempt, SDA_jumps_success
  integer :: SDA_translations_attempt, SDA_translations_success
  integer :: SN_SDA_swaps_attempt, SN_SDA_swaps_success
 
end module variables
  
  
