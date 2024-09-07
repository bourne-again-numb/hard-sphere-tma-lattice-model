module input_SiO4_SDA  
  use variables     
  implicit none
contains
  subroutine input
    implicit none
    integer :: ioerr,status
    integer :: i,j,manual_concen_flag
    real :: a,b,c,alpha
    
    qn_point = 11
    cluster_threshold = 10
    
    open(101,file='input.dat',status='unknown',iostat=ioerr)
    if(ioerr.ne.0) then
       print*,"Error: Reading input.dat file"
       stop
    end if
    
    !read lattice dimensions
    read(101,*,iostat=status) lx, ly, lz
    if(status.ne.0)then
       print*,'error reading lattice dimensions'
       stop
    end if
    
    !total number of sites and unit cells
    nsites = 2*lx*ly*lz
    unitcells = lx*ly*lz

    !read the temperature 
    read(101,*,iostat=status) tstar,sda_size
    if(tstar.lt.0.0d0.or.status.ne.0)then
       print*,'error reading tempersture'
       stop
    end if

    !read the no. of various sweeps
    read(101,*,iostat=status) nsweeps,nprint,nexchange
    if(status.ne.0)then
       print*,'error reading sweeps'
       stop
    end if
    !neqsweeps = int(0.7 * nsweeps)

    !read the concentrations
    read(101,*,iostat=status)nTEOS,nSDA,nH2O,manual_concen_flag
    !if manual_concen_flag is zero then
    !read the numbers in the next line as nTEOS,nSDAOH,nH2O
    !if the flag is non zero then
    !nTEOS => nSN
    !nSDAOH => nSI
    !nH2O => nSDA
    !nSN,nSI,nSDA are the molecule number of the respective species
        
    if(manual_concen_flag.eq.0)then!if manual concentration flag is zero

       !do the usual calculation of pH, ni, xi
       call molar_to_pH(alpha)
       call molar_fractions(alpha,xi(1),xi(2),xi(3),xi(4))
       call fractions_number(xi(1),xi(2),xi(3),xi(4),ni(1),ni(2),ni(3),ni(4))
       
    elseif(manual_concen_flag.eq.1)then!if the manual concentration flag is non zero

       !negelct the nH2O read from the file. only the first two values of nTEOS and nSDAOH
       !are important. Refer to job_queue_molecules.sh for more information.
       ni(1) = int(0)
       ni(2) = nTEOS
       ni(3) = nSDA
       ni(4) = nsites - (ni(1)+ni(2)+ni(3))
       xi = 1.0d0*ni/real(unitcells)
       xi(4) = 1.00000d0 - (xi(1)+xi(2)+xi(3))
       
    else
       print*,'ERROR: incorrect manual_config_flag :: subroutine input'
       stop
    end if!endif manual concentration flag is zero
        
    !particle numbers for manual_check
    !make sure that the system size is 60,60,60
    !ni(1) = 200
    !ni(2) = 200
    !ni(3) = 0
    !ni(4) = nsites - ni(2) - ni(1) - ni(3)
        
    read(101,*,iostat=status) pen3, pen4
    if(status.ne.0) then
       print*, 'error: reading ring penalties'
       stop
    end if
    
    read(101,*,iostat=status) eOSDA,eSNSN
    if(status.ne.0) then
       print*, 'error: reading interactions'
       stop
    end if
    
    read(101,*,iostat=status) tstar_final,m0,m1,m2
    if(status.ne.0)then
       print*,'error: reading temperature profile data'
       stop
    end if
    tstar_initial=tstar
    up_rate = real((tstar_final-tstar_initial)/(m1-m0))
    down_rate = 8.000

    read(101,*,iostat=status) initial_config_flag,time_limit
    if(status.ne.0)then
       print*,'error: reading initial configuration flag'
       stop
    end if

    read(101,*,iostat=status) ftemp, nexchange, jump_flag
    if(status.ne.0)then
       print*,'error: reading ftemp value'
       stop
    end if
    
    close(101)
    
    !total atoms on the lattice
    natom = ni(1)*5 + ni(2)*5 + ni(3)
    !all the monomers/molecules on the lattice
    allmolecules = ni(1) + ni(2) + ni(3)
    !allmonomers = allmolecules - ni(3)
        
    allocate(occupancy(nsites))
    allocate(tag_site(nsites * 4))
    
    allocate(noccupy(allmolecules))
    allocate(resite(nsites))
    allocate(head(nsites))
    allocate(spin(nsites))
    allocate(clabel(nsites))
    allocate(csize(nsites))
    
    allocate(n1list(8,nsites))
    allocate(n2list(6,nsites))
    allocate(n3list(12,nsites))
    allocate(n4list(8,nsites))
    allocate(n5list(6,nsites))    

    allocate(rx(nsites))
    allocate(ry(nsites))
    allocate(rz(nsites))

    !custom size of SDA and its neighbors/neighbor list
    ! here neighbors is sda_size+1 because the BO interactions 
    ! with SDA are at sda_size+1 
    allocate(neighbors(sda_size+1))
    allocate(nlist(sda_size+1,nsites,nsites))

    SN_jumps_attempt = 0; SN_jumps_success = 0
    SN_rotations_attempt = 0; SN_rotations_success = 0
    SDA_jumps_attempt = 0; SDA_jumps_success = 0
    SDA_translations_attempt = 0; SDA_translations_success = 0
    SN_SDA_swaps_attempt = 0; SN_SDA_swaps_success = 0
    
  end subroutine input
  !------------------------------------------------------------                       
  !--- molar_pH(x,y,z,pH): converts molar ratios to pH values                         
  !------------------------------------------------------------                       
  subroutine molar_to_pH(alpha)
    use variables
    implicit none

    real,intent(out) :: alpha
    real :: y, x, z!mole proportions                                                  
    real :: a1,b1,c1!quadratic equation coefficients                                  
    real*8,parameter :: c0=real(55.55),K=real(1.75*10**6),pKw = real(14)
    real :: alpha1,alpha2,beta,pOH

    y = nTEOS;x=nSDA;z=nH2O

    if(x.eq.real(0))then
       pH = real(2.5)
       alpha = real(0)
       return
    end if


    !concentrated system                                                              
    !quadratic equation is: (K-1)alpha^2 - (x*K + y*K + 4*y + z)alpha + x*y*k = 0     
    a1 = K-real(1)
    b1 = real(-1)*(x*K + y*K + real(4)*y + z)
    c1 = x*y*K

    !dilute system                                                                    
    !quadratic equation is: Kalpha^2 - (x*K+y*K+1) + K*x*y = 0                        
    !a1 = K                                                                           
    !b1 = real(-1)*(x*K + y*K + 1)                                                    
    !c1 = x*y*K                                                                       

    alpha1 = (-b1 + sqrt(b1**2 - real(4)*a1*c1))/real(2*a1)
    alpha2 = (-b1 - sqrt(b1**2 - real(4)*a1*c1))/real(2*a1)

    alpha = min(alpha1,alpha2)

    pH = pKw + log10((x - alpha) * c0)

    if((alpha.lt.0).or.(alpha.gt.x))then
       print*,'alpha1',alpha1
       print*,'alpha2',alpha2
       print*,'x',x
       print*,'y',y
       print*,'z',z
       stop
    end if

    return
  end subroutine molar_to_pH
  !------------------------------------------------------------
  !--- pH_to_molar(y,x,z,pH): converts pH to molar ratios
  !------------------------------------------------------------
  subroutine ph_to_molar(y,x,z,alpha,pH)
    implicit none
    
    real,intent(out) :: x, alpha!mole proportions
    real,intent(in) :: pH, y, z
    real :: a,b,c!quadratic equation coefficients
    real,parameter :: c0=real(55.55),K=real(1.75*10**6),pKw = real(14)
    real :: alpha1,alpha2,beta,pOH
    
    if(pH.eq.real(2.5))then
       x = real(0)
       alpha = real(0)
       return
    end if
    
    pOH = pKw - pH
    beta = (real(10)**(-pOH))/c0
    
    a = real(1)
    b = beta*K + real(4)*y + z
    c = real(-1)*beta*y*K
    
    alpha1 = (-b + sqrt(b**2 - real(4)*a*c))/real(2*a)
    alpha2 = (-b - sqrt(b**2 - real(4)*a*c))/real(2*a)
    
    if((alpha1.gt.0).and.(alpha2.lt.0))then
       alpha = alpha1
    elseif((alpha1.lt.0).and.(alpha2.gt.0))then
       alpha = alpha2
    else
       print*,'Error'
       stop
    end if
    
    x = alpha + beta
    
    return
  end subroutine ph_to_molar
  !------------------------------------------------------------
  !--- molar_fractions(y,x,z,alpha,xSI,xSN,xSDA,xH2O): converts molar ratios to fractions
  !------------------------------------------------------------
  subroutine molar_fractions(alpha,xSI,xSN,xSDA,xH2O)
    implicit none
    
    real,intent(in) :: alpha
    real :: y,x,z
    real,intent(out) :: xSI,xSN,xSDA,xH2O
    
    y = nTEOS;x=nSDA;z=nH2O
    
    if(x.eq.real(0))then
       
       xSN = (y)/(real(5)*y + z)
       xSI = real(0)
       xSDA = real(0)
       xH2O = real(1) - xSN - xSI - xSDA
       
    else
       
       xSN = (y-alpha)/(x + real(5)*y + z)
       xSI = (alpha)/(x + real(5)*y + z)
       xSDA = (x)/(x + real(5)*y + z)
       xH2O = real(1) - xSN - xSI - xSDA
       
    end if
    
    return
  end subroutine molar_fractions
  !------------------------------------------------------------
  !--- fractions_number(xSI,xSN,xSDA,xH2O,nSI,nSN,nSDA,nH2O): converts fractions to numbers
  !------------------------------------------------------------
  subroutine fractions_number(xSI,xSN,xSDA,xH2O,nSI,nSN,nSDA,nH2O)
    implicit none 

    integer,intent(out) :: nSI,nSN,nSDA,nH2O
    real,intent(in) :: xSI,xSN,xSDA,xH2O
    integer :: nsites
    real :: temp_SI,temp_SN,temp_SDA
    
    nsites = real(lx) * real(ly) * real(lz)
    temp_SI = xSI * nsites
    nSI = nint(temp_SI)
    temp_SN = xSN * nsites
    nSN = nint(temp_SN)
    temp_SDA = xSDA * nsites
    nSDA = nint(temp_SDA)
    nH2O = nsites - nSI - nSN - nSDA
        
    return
  end subroutine fractions_number
  !*************************************************************
end module input_SiO4_SDA
