module grid_m
    use dataset_m
    use json_m
    use atmosphere_m
    use panel_m
    implicit none
    
    
    type grid_t
        character(100) :: master_filename
        character(100) :: tagName
        character(100) :: tagUUID
        character(100) :: tagDate
        character(1000) :: afpath
        character(20) :: aftype
        character(20) :: afcode

        type(json_file) :: json     ! the JSON structure read from the file:
        type(panel_t) :: panel      ! panel representation for potential flow solution
    
        real :: afLength
        real :: alpha               ! angle of attack used to generate the grid
        real :: altitude
        real :: velocity
        real :: ReynoldsNumber
        real :: MachNumber
        real :: temperature
        real :: density
        real :: pressure
        real :: minBeta

        integer :: ni       !wraps around surface
        integer :: nj       !radial

        integer :: nSurfaceNodes
        integer :: nWakeNodes
        integer :: nRadialNodes
        
        character(10) :: teType

        character(10) :: flapType
        real :: flapFraction
        real :: flapDeflection
        real :: flapHingeVertPos

        real,allocatable,dimension(:) :: mergeFactor
        real,allocatable,dimension(:,:) :: surfacePoints
        real,allocatable,dimension(:,:,:) :: organicGrid,algebraicGrid,masterGrid,mediumGrid,coarseGrid
    end type grid_t

contains

!-----------------------------------------------------------------------------------------------------------
subroutine grid_allocate(t)
    type(grid_t) :: t

    t%ni = 2*t%nWakeNodes + t%nSurfaceNodes
    t%nj = t%nRadialNodes

    allocate(t%mergeFactor(t%nj))
    allocate(t%surfacePoints(t%nSurfaceNodes,2))
    allocate(t%organicGrid(  t%ni, t%nj, 2))
    allocate(t%algebraicGrid(t%ni, t%nj, 2))
    allocate(t%masterGrid(   t%ni, t%nj, 2))
    allocate(t%mediumGrid((t%ni-1)/2+1, t%nj, 2))
    allocate(t%coarseGrid((t%ni-1)/4+1, t%nj, 2))
end subroutine grid_allocate

!-----------------------------------------------------------------------------------------------------------
subroutine grid_deallocate(t)
    type(grid_t) :: t
    deallocate(t%surfacePoints)
    deallocate(t%organicGrid)
    deallocate(t%algebraicGrid)
    deallocate(t%masterGrid)
    deallocate(t%mediumGrid)
    deallocate(t%coarseGrid)
end subroutine grid_deallocate

!-----------------------------------------------------------------------------------------------------------
subroutine grid_set_defaults(t)
    type(grid_t) :: t    

    t%minBeta = 1.000000001

end subroutine grid_set_defaults

!-----------------------------------------------------------------------------------------------------------
subroutine grid_load_json(t)
    type(grid_t) :: t
    integer :: loc
    character(len=:),allocatable :: cval

    write(*,*) 'reading input file: ',t%master_filename
    call t%json%load_file(filename = t%master_filename); call json_check()
    loc = index(t%master_filename,'.json')
    t%master_filename = t%master_filename(1:loc-1) !deletes the .json file extension

    call t%json%get('tag.name',          cval);               call json_check(); t%tagName = trim(cval)
    call t%json%get('tag.UUID',          cval);               call json_check(); t%tagUUID = trim(cval)
    call t%json%get('tag.date',          cval);               call json_check(); t%tagDate = trim(cval)
    call t%json%get('mesh.surface.nCells',t%nSurfaceNodes);   call json_check()
    call t%json%get('mesh.wake.nCells',      t%nWakeNodes);   call json_check()
    call t%json%get('mesh.radial.nCells',  t%nRadialNodes);   call json_check()
    call t%json%get('airfoil.length',          t%afLength);   call json_check()

    t%nSurfaceNodes = t%nSurfaceNodes + 1
    t%nRadialNodes = t%nRadialNodes + 1
    call grid_allocate(t)
end subroutine grid_load_json

!-----------------------------------------------------------------------------------------------------------
subroutine grid_set_conditions(t)
    type(grid_t) :: t
    type(atmosphere_t) :: atm
    character(len=:),allocatable :: cval
    real :: value,rho,mu,a,ans(7)

    call t%json%get('condition.alpha',     t%alpha);  call json_check()
    t%alpha = t%alpha*pi/180.0
    call t%json%get('condition.altitude',     t%altitude);   call json_check()
    call t%json%get('condition.velocityType',       cval);   call json_check()
    call t%json%get('condition.velocityValue',     value);   call json_check()
    
    call atm_create(atm)
    call atm_get_properties(atm,t%altitude,ans)
    rho = ans(5)
    mu = ans(6)
    a = ans(7)

    select case (trim(cval))
        case('Velocity')
            t%velocity = value;
            t%ReynoldsNumber = rho*t%velocity*t%afLength/mu
            t%MachNumber = t%velocity/a
        case('ReynoldsNumber')
            t%ReynoldsNumber = value;
            t%velocity = t%ReynoldsNumber*mu/rho/t%afLength
            t%MachNumber = t%velocity/a
        case('MachNumber')
            t%MachNumber = value;
            t%velocity = t%MachNumber*a
            t%ReynoldsNumber = rho*t%velocity*t%afLength/mu
        case default
            write(*,*) 'Invalid velocity specification! Aborting program.'
            stop
    end select

    t%temperature = ans(3);
    t%pressure = ans(4);
    t%density = ans(5);
    
    write(*,*) '-------------------- Condition --------------------'
    write(*,*) '     Characteristic Length (m) = ',t%afLength
    write(*,*) '     Freestream Velocity (m/s) = ',t%velocity
    write(*,*) '        Freestream Mach Number = ',t%MachNumber
    write(*,*) '               Reynolds Number = ',t%ReynoldsNumber
    write(*,*) '                  Altitude (m) = ',t%altitude
    write(*,*) '---------------------------------------------------'

end subroutine grid_set_conditions

!-----------------------------------------------------------------------------------------------------------
subroutine grid_save_surface_geometry(t)
    type(grid_t) :: t
    character(100) :: filename
    integer :: ios,i
    real :: CL,CmLE,Cmc4

    call panel_coefficients(t%panel,t%MachNumber,CL,CmLE,Cmc4)
    write(*,*)
    write(*,*) '--------- Inviscid Results -----------'
    write(*,*) '     alpha (deg) = ',t%alpha*180.0/pi
    write(*,*) '              CL = ',CL
    write(*,*) '         Cm(c/4) = ',Cmc4
    write(*,*) '--------------------------------------'

    !Write results to file
    filename = trim(t%master_filename)//'_out.json'

    write(*,*) 'Saving ',trim(t%tagName),' airfoil to file: ',trim(filename)
    open(unit = 100, File = trim(filename), action = "write", iostat = ios)
    write(100,'(a)')           '{'
    write(100,'(a)')           '  "airfoil" : {'
    write(100,'(a)')           '      "name" : "'//trim(t%tagName)//'",'
    write(100,'(a)')           '      "inviscid" : {'
     write(100,'(a,ES25.16,a)')'          "alpha" : ',t%alpha*180.0/pi,','
     write(100,'(a,ES25.16,a)')'          "CL" : ',CL,','
     write(100,'(a,ES25.16,a)')'          "Cmc4" : ',Cmc4,''
    write(100,'(a)')           '      },'
    write(100,'(a)')           '      "geometry" : {'
    write(100,'(a,I10,a)')     '          "npoints" : ',t%nSurfaceNodes,','
    write(100,'(a)')           '          "xpoints" : ['
    do i=1,t%nSurfaceNodes-1
     write(100,'(a,ES25.16,a)')'                       ',t%surfacePoints(i,1),','
    end do
    write(100,'(a,ES25.16,a)') '                       ',t%surfacePoints(t%nSurfaceNodes,1),'],'
    write(100,'(a)')           '          "ypoints" : ['
    do i=1,t%nSurfaceNodes-1
     write(100,'(a,ES25.16,a)')'                       ',t%surfacePoints(i,2),','
    end do
    write(100,'(a,ES25.16,a)') '                       ',t%surfacePoints(t%nSurfaceNodes,2),']'
    write(100,'(a)')           '      }'
    write(100,'(a)')           '   }'
    write(100,'(a)')           '}'

    close(100)
    
end subroutine grid_save_surface_geometry

!-----------------------------------------------------------------------------------------------------------
subroutine grid_characterize_airfoil(t)
    type(grid_t) :: t
    character(100) :: filename
    integer :: ios,i
    real :: alphaL0, CL0,Cm0,CD0,CLa,Cma
    real :: alpha,y,yprime,dx,dummy,Cm1,Cm2,CL1,CL2

    !Characterize Airfoil
    dx =   0.0001
    alpha = 0.0 !initial guess
    y = grid_get_airfoil_CL(t,alpha)
    do while(abs(y) > 1.0e-13)
        yprime = 0.5*(grid_get_airfoil_CL(t,alpha+dx) - grid_get_airfoil_CL(t,alpha-dx))/dx
        alpha = alpha - y/yprime
        y = grid_get_airfoil_CL(t,alpha)
    end do
    alphaL0 = alpha
    
    t%panel%alpha = alphaL0
    call panel_solve(t%panel)
    call panel_coefficients(t%panel,t%MachNumber,CL0,dummy,Cm0)
    CD0 = 0.0
    
    t%panel%alpha = alphaL0 - dx
    call panel_solve(t%panel)
    call panel_coefficients(t%panel,t%MachNumber,CL1,dummy,Cm1)
    t%panel%alpha = alphaL0 + dx
    call panel_solve(t%panel)
    call panel_coefficients(t%panel,t%MachNumber,CL2,dummy,Cm2)
    CLa = 0.5*(CL2 - CL1)/dx
    Cma = 0.5*(Cm2 - Cm1)/dx

    write(*,*)
    write(*,*) '-------- Computed Airfoil Characteristics ---------'
    write(*,*) '   Zero-Lift Alpha (rad) = ',alphaL0
    write(*,*) '   Zero-Lift Alpha (deg) = ',alphaL0*180.0/pi
    write(*,*) '            Zero-Lift CL = ',CL0
    write(*,*) '            Zero-Lift Cm = ',Cm0
    write(*,*) '            Zero-Lift CD = ',CD0
    write(*,*) '                CL,alpha = ',CLa
    write(*,*) '                Cm,alpha = ',Cma

    !Write results to file
    filename = trim(t%master_filename)//'_out.json'

!    !Set Surface Points so airfoil has length = 1.0
    t%surfacePoints(:,:) = t%surfacePoints(:,:) / t%afLength

    write(*,*) 'Saving ',trim(t%tagName),' airfoil to file: ',trim(filename)
    open(unit = 100, File = trim(filename), action = "write", iostat = ios)
    write(100,'(a)')           '{'
    write(100,'(a)')           '  "airfoil" : {'
    write(100,'(a)')           '      "name" : "'//trim(t%tagName)//'",'
    write(100,'(a)')           '      "aeroData" : {'
    write(100,'(a)')           '          "default" : {'
    write(100,'(a)')           '              "type" : "coeffs",'
    write(100,'(a)')           '              "Mach" : 0.0,'
    write(100,'(a)')           '              "ReynoldsNumber" : 0.0,'
    write(100,'(a,ES25.16,a)') '              "zeroLiftAlpha" : ',alphaL0,','
    write(100,'(a,ES25.16,a)') '              "zeroLiftCm" : ',Cm0,','
    write(100,'(a,ES25.16,a)') '              "zeroLiftCD" : ',CD0,','
    write(100,'(a,ES25.16,a)') '              "CLalpha" : ',CLa,','
    write(100,'(a,ES25.16,a)') '              "Cmalpha" : ',Cma,','
    write(100,'(a)')           '              "units" : "radians"'
    write(100,'(a)')           '          }'
    write(100,'(a)')           '      },'
    write(100,'(a)')           '      "geometry" : {'
    write(100,'(a,I10,a)')     '          "npoints" : ',t%nSurfaceNodes,','
    write(100,'(a)')           '          "xpoints" : ['
    do i=1,t%nSurfaceNodes-1
     write(100,'(a,ES25.16,a)')'                       ',t%surfacePoints(i,1),','
    end do
    write(100,'(a,ES25.16,a)') '                       ',t%surfacePoints(t%nSurfaceNodes,1),'],'
    write(100,'(a)')           '          "ypoints" : ['
    do i=1,t%nSurfaceNodes-1
     write(100,'(a,ES25.16,a)')'                       ',t%surfacePoints(i,2),','
    end do
    write(100,'(a,ES25.16,a)') '                       ',t%surfacePoints(t%nSurfaceNodes,2),']'
    write(100,'(a)')           '      }'
    write(100,'(a)')           '   }'
    write(100,'(a)')           '}'

    close(100)
    
    !Set Surface Points back to correct length of airfoil
    t%surfacePoints(:,:) = t%surfacePoints(:,:) * t%afLength


end subroutine grid_characterize_airfoil
 
!-----------------------------------------------------------------------------------------------------------
real function grid_get_airfoil_CL(t,alpha)
    type(grid_t) :: t
    real :: alpha,CL,CmLE,Cmc4

    t%panel%alpha = alpha
    call panel_solve(t%panel)
    call panel_coefficients(t%panel,t%MachNumber,CL,CmLE,Cmc4)
    grid_get_airfoil_CL = CL
        
end function grid_get_airfoil_CL

!-----------------------------------------------------------------------------------------------------------
subroutine grid_create_airfoil_surface(t)
    type(grid_t) :: t
    type(dataset_t) :: afdata
    type(dataset_t) :: afdataFlap
    type(panel_t) :: temp
    character(len=:),allocatable :: cval
    real :: mindist,interpolate,lepercent,xilength,mypercent,theta
    character(6) :: clusterType
    real :: AA,BB,CC,beta,tau,A,ybar,scale
    real :: Amat(3,3), Bvec(3), coeffs(3)
    integer :: i, j, mincoord, nNodes
    real :: CL,CmLE,Cmc4
    integer :: local_npoints
    real,allocatable,dimension(:,:) :: local_points
    character(len=10) :: str

    real :: foil(2), myxi, XiSeg(4)
    integer :: Seg1nodes, Seg2nodes, Seg3nodes, Seg4nodes

    !Read Airfoil Specification
    call t%json%get('airfoil.type',                      cval);     call json_check(); t%aftype = trim(cval)
    call t%json%get('airfoil.teType',                    cval);     call json_check(); t%teType = trim(cval)

    !Create Basic airfoil
    select case (trim(t%aftype))
        case ('Local')
            call t%json%get('geomPoints.npoints',   local_npoints);     call json_check()
            allocate(local_points(local_npoints,2))
            do i=1,local_npoints
                call integer_to_string(i,str)
                call t%json%get('geomPoints.xpoints('//trim(str)//')',   local_points(i,1));     call json_check()
                call t%json%get('geomPoints.ypoints('//trim(str)//')',   local_points(i,2));     call json_check()
            end do
            call ds_create_from_data(afdata,local_npoints,2,local_points(:,:))

        case ('File')
            call t%json%get('airfoil.code',                  cval);     call json_check(); t%afpath = trim(cval)
            call ds_create_from_file(afdata,t%afpath,2)
            
        case ('NACA 4-Digit')
            call t%json%get('airfoil.code',                  cval);     call json_check(); t%afcode = trim(cval)
            temp%naca = trim(t%afcode) !Either NACA 4-digit or uniform-load NACA
            write(*,*) 'Creating NACA ',temp%naca(1:4),' airfoil.'
            read(temp%naca(1:2),*) temp%load
            if(temp%load.eq.'UL') then
                call t%json%get('airfoil.designCL',              temp%CLd);     call json_check()
                write(*,*) '    Design Lift Coefficient = ',temp%CLd
            end if
            temp%npts = t%nSurfaceNodes - 1 !must be an even number
            if(trim(t%teType).eq.'closed') then
                temp%OpenTE = 0 !closed trailing edge
            else
                temp%OpenTE = 1
            end if
            temp%PointsPlacement = 2
            call panel_create_from_naca(temp) !used only to get NACA 4-digit geometry
            call ds_create_from_data(afdata,temp%npts,2,temp%Points(:,:))
    end select

    if(trim(t%teType).eq.'closed') then !Close trailing edge
        afdata%RawData(1,:) = 0.5*(afdata%RawData(1,:) + afdata%RawData(afdata%datasize,:))
        afdata%RawData(afdata%datasize,:) = afdata%RawData(1,:)
    end if
    !Scale Airfoil to length = 1.0
    scale = max(afdata%RawData(1,1),afdata%RawData(afdata%datasize,1))
    afdata%RawData(:,:) = afdata%RawData(:,:)/scale
    call ds_calc_Xi(afdata)

    !At this point, the airfoil TE is at x=1.0
    call ds_cubic_setup(afdata,0,2,0.0,2,0.0)
    call ds_print_data(afdata)

    !Find Leading-Edge
    mindist = math_length(2,[0.0,0.0],afdata%RawData(1,:))
    mincoord = 1
    do i=1,afdata%datasize
        if(math_length(2,[0.0,0.0],afdata%RawData(i,:)) < mindist) then
            mindist = math_length(2,[0.0,0.0],afdata%RawData(i,:))
            mincoord = i
        end if
    end do
    xilength = afdata%Xi(afdata%datasize)
    lepercent = afdata%Xi(mincoord)/xilength
    write(*,*) "Leading Edge at i = ",mincoord, ", Distance from origin = ", mindist
    write(*,*) "Distance along airfoil spline = ",afdata%Xi(mincoord), " = ",lepercent*100.0,"%"
    
    !Read Cluster and Flap Information
    call t%json%get('mesh.surface.interpolate',     interpolate);  call json_check()
    call t%json%get('mesh.surface.clusterBeta',     beta);  call json_check()
    call t%json%get('airfoil.flapType',              cval);     call json_check(); t%flapType = trim(cval)
    beta = 1.0+ beta;
    beta = max(t%minBeta,beta);

    !Create Surface Points
    select case (trim(t%flapType))

        case ('None')
            tau = 1.0 !equally spaces about leading edge and trailing edge
            nNodes = (t%nSurfaceNodes+1)/2
            do i=1,nNodes
                mypercent = real(i-1)/real(nNodes-1)
                !Alley, p. 38
                A = ((beta+1.0)/(beta-1.0))**((2.0*mypercent-tau)/(2.0-tau))
                ybar = (A*(tau+beta) + (tau-beta))/(A+1.0)/(tau+1.0)
                call ds_weighted_interpolate(afdata,ybar*lepercent*xilength,interpolate,t%surfacePoints(i,:))
                call ds_weighted_interpolate(afdata,(1.0-ybar*(1.0-lepercent))*xilength,interpolate,&
                                            &t%surfacePoints(t%nSurfaceNodes+1-i,:))
!                write(*,*) i,mypercent*lepercent, (1.0-mypercent*(1.0-lepercent))
            end do
        case('Sealed')
            call t%json%get('airfoil.flapX',     t%flapFraction);    call json_check()
            call t%json%get('airfoil.flapY',     t%flapHingeVertPos);  call json_check()
            call t%json%get('airfoil.flapDeflection',   t%flapDeflection);  call json_check()

            t%flapFraction = 1.0-t%flapFraction;
            call grid_create_deflected_airfoil(t,afdata,lepercent,afdataFlap,XiSeg)
            
!            write(*,*) 'Seg4',XiSeg

            ! 4 Sections
            ! Allocate points by percentage of length along spline
            Seg1nodes = nint(real(t%nSurfaceNodes) * XiSeg(1)/XiSeg(4)) 
            Seg2nodes = nint(real(t%nSurfaceNodes) * (XiSeg(2) - XiSeg(1))/XiSeg(4)) + 1
            Seg3nodes = nint(real(t%nSurfaceNodes) * (XiSeg(3) - XiSeg(2))/XiSeg(4)) + 1 
            Seg4nodes = t%nSurfaceNodes - Seg3nodes - Seg2nodes - Seg1nodes + 3       

            write(*,*) Seg1nodes, Seg2nodes,Seg3nodes,Seg4nodes,t%nSurfaceNodes
            ! Each node : Percentage, A, ybar

            tau = 1.0
            j = 1
            do i = 1, Seg1nodes
                mypercent = real(i-1)/real(Seg1nodes-1)
                A = ((beta+1.0)/(beta-1.0))**((2.0*mypercent-tau)/(2.0-tau))
                ybar = (A*(tau+beta) + (tau-beta))/(A+1.0)/(tau+1.0)
                call ds_weighted_interpolate(afdataFlap,ybar*XiSeg(1),interpolate,t%surfacePoints(j,:))
                write(*,*) i, j, t%surfacePoints(j,:)
                j = j + 1
            end do
            j = j - 1
            write(*,*)
            do i = 1, Seg2nodes
                mypercent = real(i-1)/real(Seg2nodes-1)
                A = ((beta+1.0)/(beta-1.0))**((2.0*mypercent-tau)/(2.0-tau))
                ybar = (A*(tau+beta) + (tau-beta))/(A+1.0)/(tau+1.0)
                call ds_weighted_interpolate(afdataFlap,ybar*(XiSeg(2)-XiSeg(1))+XiSeg(1),interpolate,t%surfacePoints(j,:))
                write(*,*) i, j, t%surfacePoints(j,:)
                j = j + 1
            end do
            j = j - 1
            write(*,*)
            do i = 1, Seg3nodes
                mypercent = real(i-1)/real(Seg3nodes-1)
                A = ((beta+1.0)/(beta-1.0))**((2.0*mypercent-tau)/(2.0-tau))
                ybar = (A*(tau+beta) + (tau-beta))/(A+1.0)/(tau+1.0)
                call ds_weighted_interpolate(afdataFlap,ybar*(XiSeg(3)-XiSeg(2))+XiSeg(2),interpolate,t%surfacePoints(j,:))
                write(*,*) i, j, t%surfacePoints(j,:)
                j = j + 1
            end do
            j = j - 1
            write(*,*)
            do i = 1, Seg4nodes
                mypercent = real(i-1)/real(Seg4nodes-1)
                A = ((beta+1.0)/(beta-1.0))**((2.0*mypercent-tau)/(2.0-tau))
                ybar = (A*(tau+beta) + (tau-beta))/(A+1.0)/(tau+1.0)
                call ds_weighted_interpolate(afdataFlap,ybar*(XiSeg(4)-XiSeg(3))+XiSeg(3),interpolate,t%surfacePoints(j,:))
                write(*,*) i, j, t%surfacePoints(j,:)
                j = j + 1
            end do

    end select
    

    !Create Inviscid Panel Case for final geometry (includes flap deflection)
    t%panel%PointsPlacement = 1
    call panel_create_from_data(t%panel,t%nSurfaceNodes,t%surfacePoints)
    t%panel%alpha = t%alpha
    call panel_solve(t%panel)
    call panel_coefficients(t%panel,t%MachNumber,CL,CmLE,Cmc4)
    write(*,*)
    write(*,*) '--------- Inviscid Results -----------'
    write(*,*) '     alpha (deg) = ',t%alpha*180.0/pi
    write(*,*) '              CL = ',CL
    write(*,*) '         Cm(c/4) = ',Cmc4
    write(*,*) '--------------------------------------'
    
    !Scale Airfoil to correct size
    t%surfacePoints(:,:) = t%surfacePoints(:,:) * t%afLength

    !Append to organic grid
!    write(*,*) 'Surface Nodes'
    do i=1,t%nSurfaceNodes
        t%organicGrid(t%nWakeNodes + i,1,:) = t%surfacePoints(i,:)
!        write(*,*) t%surfacePoints(i,:)
    end do
    
end subroutine grid_create_airfoil_surface

!-----------------------------------------------------------------------------------------------------------
subroutine grid_create_wake_surface(t)
    type(grid_t) :: t
    type(dataset_t) :: ds_wake
    real :: beta,ybar,mypercent,x,x0,A
    real :: wakeLength,dist,maxXi
    integer :: i, nPanelWakePoints
    
    write(*,*) 'Computing Wake Profile...'
    !Create Wake surface
    call t%json%get('mesh.wake.length',        wakeLength);         call json_check()
    call t%json%get('mesh.wake.clusterBeta',   beta);    call json_check()
    beta = 1.0 + beta;

    nPanelWakePoints = panel_make_wake_streamline(t%panel,wakeLength*cos(t%alpha),wakeLength)
    call ds_create_from_data(ds_wake,nPanelWakePoints,2,t%afLength*t%panel%Wake(:,:))
    call ds_cubic_setup(ds_wake,0,2,0.0,2,0.0)
!    call ds_print_data(ds_wake)
    maxXi = ds_wake%Xi(ds_wake%datasize)

    if(beta.eq.1.0) then !Calculate best clustering parameter based on cell size of airfoil surface TE
        dist = math_length(2,t%surfacePoints(t%nSurfaceNodes,:),t%surfacePoints(t%nSurfaceNodes-1,:))
        write(*,*) 
        write(*,*) 'Computing beta in wake to match spacing at trailing edge...'
        write(*,*) '   Trailing-Edge Grid Spacing = ',dist
        beta = grid_find_beta(t,maxXi,1.0/real(t%nWakeNodes),dist)
    else
        beta = max(t%minBeta,beta);    
    end if
    
!    write(*,*) 'Wake: '
    do i=1,t%nWakeNodes
        mypercent = real(i)/real(t%nWakeNodes)
        ybar = grid_get_ybar(mypercent,beta)
        call ds_cubic_interpolate(ds_wake,ybar*maxXi,0,t%organicGrid(t%nWakeNodes - (i-1), 1, :))
        t%organicGrid(t%nWakeNodes + t%nSurfaceNodes + i, 1, :) = t%organicGrid(t%nWakeNodes - (i-1), 1, :)
!        write(*,*) mypercent,ybar,ybar*maxXi,t%organicGrid(t%nWakeNodes - (i-1), 1, :)
    end do
        
end subroutine grid_create_wake_surface

!-----------------------------------------------------------------------------------------------------------
real function grid_find_beta(t,length,percent,dist)
    type(grid_t) :: t
    real :: length,percent,dist,beta,y,yprime,dx
    integer :: i

    write(*,*)
    write(*,*) '------------------ Computing Best beta -----------------------------'
    dx =   0.0000000001
    beta = t%minBeta !initial guess
    y = length*grid_get_ybar(percent,beta) - dist ! y = error in solution
    do while(abs(y) > 1.0e-12)
        yprime = 0.5*length*(grid_get_ybar(percent,beta+dx) - grid_get_ybar(percent,beta-dx))/dx
        beta = beta - y/yprime
        if((beta < 1.0) .or. (beta.ne.beta)) then
            beta = t%minBeta
            write(*,*) '---> ERROR! Unable to solve for best beta parameter for clustering. Using beta = ',t%minBeta
            exit
        end if
        y = length*grid_get_ybar(percent,beta) - dist
write(*,*) beta,y
    end do
    write(*,*) '            Computed beta = ',beta
    write(*,*) '      Target grid spacing = ',dist
    write(*,*) '   Resulting grid spacing = ',length*grid_get_ybar(percent,beta)
    write(*,*) '--------------------------------------------------------------------'
    grid_find_beta = beta

end function grid_find_beta

!-----------------------------------------------------------------------------------------------------------
real function grid_get_ybar(percent,beta)
    real :: percent,beta,A

    !Alley, p. 38
    A = ((beta+1.0)/(beta-1.0))**(1.0 - percent)
    grid_get_ybar = 1.0 - beta*(A-1.0)/(A+1.0)
        
end function grid_get_ybar

!-----------------------------------------------------------------------------------------------------------
subroutine grid_create_inner_grids(t)
!Creates both the organic and algebraic grids
    type(grid_t) :: t
    type(atmosphere_t) :: atm
    type(dataset_t) :: layer,farField
    character(len=:),allocatable :: cval
    real :: farFieldpts(9,2)
    real :: beta,dybar,mypercent,prevpercent,ybar,prevybar,dx(2),dxmag,x1(2),x2(2),myxi,A 
    integer :: i,j
    real :: wakeLength, radialLength, growthInfluence
    real :: yplus,Cf,dist,angle,rotate,smoothRatio,slope1(2),slope2(2),slope(2)
    real :: wakeEnd(2),delta(2)

    wakeEnd(:) = t%organicGrid(1,1,:)
    t%algebraicGrid(:,:,:) = t%organicGrid(:,:,:)

    !Create Wake surface
    call t%json%get('mesh.radial.length',            radialLength);           call json_check()
    call t%json%get('mesh.radial.growthInfluence',   growthInfluence);  call json_check()
    call t%json%get('mesh.wake.length',              wakeLength);           call json_check()
    radialLength = radialLength*t%afLength
    wakeLength = wakeLength*t%afLength
    
    !Create FarField
    dist = sqrt(wakeEnd(1)**2 + wakeEnd(2)**2)
    rotate = atan2(wakeEnd(2),wakeEnd(1))
    farFieldpts(1,:) = [ dist,      -radialLength]
    farFieldpts(2,:) = [ t%afLength+0.5*(dist-t%afLength),  -radialLength]
    farFieldpts(3,:) = [ 0.0,                 -radialLength]
    farFieldpts(4,:) = [ -radialLength*sin(45.0*pi/180.0),-radialLength*sin(45.0*pi/180.0)]
    farFieldpts(5,:) = [ -radialLength,        0.0]
    farFieldpts(6,:) = [farFieldpts(4,1), -farFieldpts(4,2)]
    farFieldpts(7,:) = [farFieldpts(3,1), -farFieldpts(3,2)]
    farFieldpts(8,:) = [farFieldpts(2,1), -farFieldpts(2,2)]
    farFieldpts(9,:) = [farFieldpts(1,1), -farFieldpts(1,2)]

    !Rotate farfield by angle "rotate"
    do i=1,9
        dist = sqrt(farFieldpts(i,1)**2 + farFieldpts(i,2)**2)
        angle = atan2(farFieldpts(i,2),farFieldpts(i,1))
        farFieldpts(i,1) = dist*cos(angle+rotate)
        farFieldpts(i,2) = dist*sin(angle+rotate)
    end do

    call ds_create_from_data(farField,9,2,farFieldpts)
    call ds_cubic_setup(farField,0,2,0.0,2,0.0)
    do i=1,t%ni
        call ds_cubic_interpolate(farField,farField%Xi(9)*real(i-1)/real(t%ni-1),0,t%algebraicGrid(i,t%nj,:))
    end do
    
    
    call t%json%get('mesh.radial.spacingDefinition',  cval);  call json_check()

    select case (trim(cval))
        case('Distance')
            call t%json%get('mesh.radial.distance',      dist);  call json_check()
            beta = grid_find_beta(t,radialLength,1.0/real(t%nRadialNodes-1),dist)
        case('yPlus')
            call t%json%get('mesh.radial.yplus',        yplus);  call json_check()
            Cf = 0.027/(t%ReynoldsNumber**(1.0/7.0))
            dist = yplus*t%afLength/t%ReynoldsNumber/sqrt(0.5*Cf)
            write(*,*) 'Computing beta to give correct yplus value on airfoil surface...'
            write(*,*) '     Characteristic Length (m) = ',t%afLength
            write(*,*) '     Freestream Velocity (m/s) = ',t%velocity
            write(*,*) '               Reynolds Number = ',t%ReynoldsNumber
            write(*,*) '                   Mach Number = ',t%MachNumber
            write(*,*) '                 Desired yplus = ',yplus
            write(*,*) '  Desired height of first cell = ',dist
            beta = grid_find_beta(t,radialLength,1.0/real(t%nRadialNodes-1),dist)
        case('ClusterDefinition')
            call t%json%get('mesh.radial.clusterBeta',   beta);  call json_check()
            beta = 1.0 + beta;
            beta = max(t%minBeta,beta);
        case default
            write(*,*) 'Invalid mesh.radial.spacingDefinition. Aborting program.'
            stop
    end select

    !Create Layers
    write(*,*)
    write(*,*) 'Creating Organic and Algebraic Grids...'
    write(*,*) 'beta = ',beta
    t%mergeFactor(:) = 0.0
    do j=2,t%nRadialNodes
        mypercent = real(j-1)/real(t%nRadialNodes-1)
        prevpercent = real(j-2)/real(t%nRadialNodes-1)
        ybar = grid_get_ybar(mypercent,beta)
        prevybar = grid_get_ybar(prevpercent,beta)
        dybar = ybar - prevybar

        t%mergeFactor(j) = ybar
        
        t%organicGrid(1,   j,:) = wakeEnd(:) + ybar*(farFieldpts(1,:) - wakeEnd(:))
        t%organicGrid(t%ni,j,:) = wakeEnd(:) + ybar*(farFieldpts(9,:) - wakeEnd(:))
        
        t%algebraicGrid(1,   j,:) = t%organicGrid(1,   j,:)
        t%algebraicGrid(t%ni,j,:) = t%organicGrid(t%ni,j,:)

        !calculate slopes at ends
        dx(:) = t%organicGrid(1,j,:) - t%organicGrid(1,j-1,:)
        dxmag = sqrt(dx(1)**2 + dx(2)**2)
        slope1(:) = dx(:)/dxmag

        dx(:) = t%organicGrid(t%ni,j,:) - t%organicGrid(t%ni,j-1,:)
        dxmag = sqrt(dx(1)**2 + dx(2)**2)
        slope2(:) = dx(:)/dxmag

        call ds_create_from_data(layer,t%ni,2,t%organicGrid(:,j-1,:))
!        call ds_cubic_setup(layer,0,2,0.0,2,0.0)
       ! call ds_print_data(layer)

        do i=2,t%ni-1
            myxi = layer%Xi(i)
!            call ds_weighted_interpolate(layer,myxi - t%BLInfluence,0.0,x1(:)) !can't use this without cubic setup
!            call ds_weighted_interpolate(layer,myxi + t%BLInfluence,0.0,x2(:))
            call ds_linear_interpolate(layer,myxi - mypercent*growthInfluence,x1(:))
            call ds_linear_interpolate(layer,myxi + mypercent*growthInfluence,x2(:))
            dx(:) = x2(:) - x1(:)
            dxmag = sqrt(dx(1)**2 + dx(2)**2)
            slope(1) = - dx(2)/dxmag
            slope(2) =   dx(1)/dxmag
            if((i < t%nWakeNodes) .and. (t%alpha > 0.0)) then
                smoothRatio = real(i - 1)/real(t%nWakeNodes-1)
                slope(:) = smoothRatio*(slope(:) - slope1(:)) + slope1(:)
            end if
            if((i > (t%ni-t%nWakeNodes)) .and. (t%alpha < 0.0)) then
                smoothRatio = real(t%ni - i)/real(t%nWakeNodes-1)
                slope(:) = smoothRatio*(slope(:) - slope2(:)) + slope2(:)
            end if
            t%organicGrid(  i,j,:) = t%organicGrid(  i,j-1,:) + slope(:)*radialLength*dybar
            t%algebraicGrid(i,j,:) = t%algebraicGrid(i,1,  :) + ybar*(t%algebraicGrid(i,t%nj,:) - t%algebraicGrid(i,1,:))
        end do
        call ds_deallocate(layer)
        write(*,*) '    Layer = ',j,'  Thickness = ',radialLength*dybar, '  Complete( % ) = ',mypercent*100.0
    end do
    write(*,*) 'done.'

end subroutine grid_create_inner_grids

!-----------------------------------------------------------------------------------------------------------
subroutine grid_merge(t)
    type(grid_t) :: t
    integer :: i,j
    real :: relax
    
    do j=1,t%nj
            relax = (real(j)/real(t%nj-1)) !this makes it dependento on the radial clustering
!            relax = t%mergeFactor(j)**(1.0/4.0)          !this is independent of the radial clustering
        do i=1,t%ni
            t%masterGrid(i,j,:) = t%organicGrid(i,j,:) + relax*(t%algebraicGrid(i,j,:) - t%organicGrid(i,j,:))
        end do
    end do

end subroutine grid_merge

!-----------------------------------------------------------------------------------------------------------
subroutine grid_coarsen(t)
    type(grid_t) :: t
    integer :: i,j,ii,jj
    
    !Medium Grid
    j = 0
    do jj=1,t%nj,2
        j = j + 1
        i = 0
        do ii=1,t%ni,2
            i = i + 1
            t%mediumGrid(i,j,:) = t%masterGrid(ii,jj,:)
        end do
    end do

    !Coarse Grid
    j = 0
    do jj=1,t%nj,4
        j = j + 1
        i = 0
        do ii=1,t%ni,4
            i = i + 1
            t%coarseGrid(i,j,:) = t%masterGrid(ii,jj,:)
        end do
    end do

    write(*,*)
    write(*,*) 'Master Grid: ',(t%ni-1)*(t%nj-1),' cells'
    write(*,*) 'Medium Grid: ',(t%ni+1)/2*(t%nj+1)/2,' cells'
    write(*,*) 'Coarse Grid: ',(t%ni+3)/4*(t%nj+3)/4,' cells'
    write(*,*)
end subroutine grid_coarsen


!-----------------------------------------------------------------------------------------------------------
subroutine grid_save_SU2(t,name)
    type(grid_t) :: t
    character(6) :: name
    character(100) :: filename
    character(len=:),allocatable :: cval
    integer :: ios,i,j,k,ii
    integer :: ni,nj,nSurfaceNodes,nWakeNodes,nRadialNodes,offset
    110 FORMAT (2ES25.16, I10)
    120 FORMAT (A, I10)
    
    select case (name)
        case ('master')
            ni = t%ni
            nj = t%nj
            nSurfaceNodes = t%nSurfaceNodes
            nWakeNodes = t%nWakeNodes
            nRadialNodes = t%nRadialNodes
!            filename = trim(t%tagName)//'_A_master.su2'
        case ('medium')
            ni = (t%ni+1)/2
            nj = (t%nj+1)/2
            nSurfaceNodes = (t%nSurfaceNodes+1)/2
            nWakeNodes = t%nWakeNodes/2
            nRadialNodes = (t%nRadialNodes+1)/2
!            filename = trim(t%tagName)//'_B_medium.su2'
        case ('coarse')
            ni = (t%ni+3)/4
            nj = (t%nj+3)/4
            nSurfaceNodes = (t%nSurfaceNodes+3)/4
            nWakeNodes = t%nWakeNodes/4
            nRadialNodes = (t%nRadialNodes+3)/4
!            filename = trim(t%tagName)//'_C_coarse.su2'
        case default
            write(*,*) 'Invalid Save Command'
            return
    end select

    filename = trim(t%master_filename)//'.su2' !Use this for Windows Compilation
!    filename = 'mesh.su2'                       !Use this for hpc Compilation
    write(*,*) 'Saving ',name,' mesh to file: ',trim(filename)
    open(unit = 100, File = trim(filename), action = "write", iostat = ios)
    write(100,120) 'NDIME=',2
    write(100,120) 'NELEM= ',(nj-1)*(ni-1)

    k = 0
    !On the first row, link the wake together so there is no interior surface
    j = 1
    do i=1,nWakeNodes+nSurfaceNodes-2
        write(100,*) 9, i-1, i, ni-nWakeNodes+i-1, ni-nWakeNodes+i-2, k
        k = k + 1
    end do
    i = nWakeNodes+nSurfaceNodes-1
    write(100,*) 9, i-1, nWakeNodes, ni-nWakeNodes+i-1, ni-nWakeNodes+i-2, k !last cell on trailing edge
    k = k + 1
    ii = nWakeNodes
    do i=nWakeNodes+nSurfaceNodes,ni-1 !closed wake
        write(100,*) 9, ii, ii-1, ni-nWakeNodes+i-1, ni-nWakeNodes+i-2, k
        k = k + 1
        ii = ii - 1
    end do

    !Rest of mesh
    offset = nWakeNodes + 1
    do j=2,nj-1
        do i=1,ni-1
            write(100,*) 9, ni*(j-1)+(i-1)-offset, ni*(j-1)+i-offset, ni*(j)+i-offset, ni*(j)+i-1-offset, k
            k = k + 1
        end do
    end do
    
    
    write(100,*)
    write(100,120) 'NPOIN= ',nj*ni-nWakeNodes-1
    k = 0
    select case (name)
        case ('master')
            do i=1,nWakeNodes+nSurfaceNodes-1
                write(100,110) t%masterGrid(i,1,:), k
                k = k + 1
            end do
            k = k + nWakeNodes+1
            do j=2,nj
                do i=1,ni
                    write(100,110) t%masterGrid(i,j,:),k
                    k = k + 1
                end do
            end do
        case ('medium')
            do i=1,nWakeNodes+nSurfaceNodes-1
                write(100,110) t%mediumGrid(i,1,:), k
                k = k + 1
            end do
            k = k + nWakeNodes+1
            do j=2,nj
                do i=1,ni
                    write(100,110) t%mediumGrid(i,j,:),k
                    k = k + 1
                end do
            end do
        case ('coarse')
            do i=1,nWakeNodes+nSurfaceNodes-1
                write(100,110) t%coarseGrid(i,1,:), k
                k = k + 1
            end do
            k = k + nWakeNodes+1
            do j=2,nj
                do i=1,ni
                    write(100,110) t%coarseGrid(i,j,:),k
                    k = k + 1
                end do
            end do
    end select

    write(100,*)
    write(100,120) 'NMARK= 3'
    write(100,120) 'MARKER_TAG= airfoil'
    write(100,120) 'MARKER_ELEMS= ',nSurfaceNodes-1
    do i=1,nSurfaceNodes-2
        write(100,*) 3,nWakeNodes+i-1, nWakeNodes+i
    end do
    write(100,*) 3,nWakeNodes+nSurfaceNodes-2, nWakeNodes
    
    write(100,120) 'MARKER_TAG= farfield'
    write(100,120) 'MARKER_ELEMS= ',ni-1
    do i=1,ni-1
        write(100,*) 3, (nj-1)*(ni)+i-1-offset, (nj-1)*(ni)+i-offset
    end do

    write(100,120) 'MARKER_TAG= outlet'
    write(100,120) 'MARKER_ELEMS= ',(nj-1)*2
    write(100,*) 3, 0, ni-offset
    do j=2,nj-1
        write(100,*) 3, (j-1)*ni-offset, j*ni-offset
    end do

    write(100,*) 3, 0, 2*ni-1-offset
    do j=2,nj-1
        write(100,*) 3, j*ni-1-offset, (j+1)*ni-1-offset
    end do

    close(100)
    
end subroutine grid_save_SU2

!-----------------------------------------------------------------------------------------------------------
subroutine grid_save_VTK(t,name)
    type(grid_t) :: t
    character(6) :: name
    character(100) :: filename
    character(len=:),allocatable :: cval
    integer :: ios,i,j,k,ii
    integer :: ni,nj,nSurfaceNodes,nWakeNodes,nRadialNodes,offset
    130 FORMAT (3ES25.16)
    
    select case (name)
        case ('master')
            ni = t%ni
            nj = t%nj
            nSurfaceNodes = t%nSurfaceNodes
            nWakeNodes = t%nWakeNodes
            nRadialNodes = t%nRadialNodes
!            filename = trim(t%tagName)//'_A_master.vtk'
        case ('medium')
            ni = (t%ni+1)/2
            nj = (t%nj+1)/2
            nSurfaceNodes = (t%nSurfaceNodes+1)/2
            nWakeNodes = t%nWakeNodes/2
            nRadialNodes = (t%nRadialNodes+1)/2
!            filename = trim(t%tagName)//'_B_medium.vtk'
        case ('coarse')
            ni = (t%ni+3)/4
            nj = (t%nj+3)/4
            nSurfaceNodes = (t%nSurfaceNodes+3)/4
            nWakeNodes = t%nWakeNodes/4
            nRadialNodes = (t%nRadialNodes+3)/4
!            filename = trim(t%tagName)//'_C_coarse.vtk'
        case default
            write(*,*) 'Invalid Save Command'
            return
    end select

    filename = trim(t%master_filename)//'.vtk'
    write(*,*) 'Saving ',name,' mesh to file: ',trim(filename)

    write(*,*) 'Saving ',name,' mesh to file: ',trim(filename)
    open(unit = 100, File = trim(filename), action = "write", iostat = ios)
    write(100,'(a)') '# vtk DataFile Version 2.0'
    write(100,'(a)') 'Written by CloudFoil 1.0'
    write(100,'(a)') 'ASCII'
    write(100,'(a)') 'DATASET STRUCTURED_GRID'
    write(100,'(a,3I10)') 'DIMENSIONS ',ni,nj,1
    write(100,'(a,I10,a)') 'POINTS ',ni*nj,' double'

    !Write Points
    select case (name)
        case ('master')
            do j=1,nj
                do i=1,ni
                    write(100,130) t%masterGrid(i,j,:),0.0
                end do
            end do
        case ('medium')
            do j=1,nj
                do i=1,ni
                    write(100,130) t%mediumGrid(i,j,:),0.0
                end do
            end do
        case ('coarse')
            do j=1,nj
                do i=1,ni
                    write(100,130) t%coarseGrid(i,j,:),0.0
                end do
            end do
    end select

    close(100)
    
end subroutine grid_save_VTK

!-----------------------------------------------------------------------------------------------------------
subroutine grid_save_web(t,name)
    type(grid_t) :: t
    character(6) :: name
    character(100) :: filename
    character(len=:),allocatable :: cval
    integer :: ios,i,j,k,ii
    integer :: ni,nj,nSurfaceNodes,nWakeNodes,nRadialNodes,offset
    140 FORMAT (ES14.6,a)
    
    select case (name)
        case ('master')
            ni = t%ni
            nj = t%nj
            nSurfaceNodes = t%nSurfaceNodes
            nWakeNodes = t%nWakeNodes
            nRadialNodes = t%nRadialNodes
!            filename = trim(t%tagName)//'_A_master.vtk'
        case ('medium')
            ni = (t%ni+1)/2
            nj = (t%nj+1)/2
            nSurfaceNodes = (t%nSurfaceNodes+1)/2
            nWakeNodes = t%nWakeNodes/2
            nRadialNodes = (t%nRadialNodes+1)/2
!            filename = trim(t%tagName)//'_B_medium.vtk'
        case ('coarse')
            ni = (t%ni+3)/4
            nj = (t%nj+3)/4
            nSurfaceNodes = (t%nSurfaceNodes+3)/4
            nWakeNodes = t%nWakeNodes/4
            nRadialNodes = (t%nRadialNodes+3)/4
!            filename = trim(t%tagName)//'_C_coarse.vtk'
        case default
            write(*,*) 'Invalid Save Command'
            return
    end select

    filename = trim(t%master_filename)//'_out.json'
    write(*,*) 'Saving ',name,' mesh to file: ',trim(filename)
    
    open(unit = 100, File = trim(filename), action = "write", iostat = ios)


    write(100,'(a)')           '{'
    write(100,'(a)')           '  "airfoil" : {'
    write(100,'(a)')           '      "name" : "'//trim(t%tagName)//'",'
    write(100,'(a)')           '      "mesh" : {'
    write(100,'(a,I10,a)')     '          "ni" : ',ni,','
    write(100,'(a,I10,a)')     '          "nj" : ',nj,','
    write(100,'(a,I10,a)')     '          "npoints" : ',ni*nj,','
    write(100,'(a)')           '          "xpoints" : ['
    !Write x points
    select case (name)
        case ('master')
            do j=1,nj
                do i=1,ni
                    if(i*j<ni*nj) then
                        write(100,140) t%masterGrid(i,j,1),','
                    else
                        write(100,140) t%masterGrid(i,j,1)
                    end if
                end do
            end do
        case ('medium')
            do j=1,nj
                do i=1,ni
                    if(i*j<ni*nj) then
                        write(100,140) t%mediumGrid(i,j,1),','
                    else
                        write(100,140) t%mediumGrid(i,j,1)
                    end if
                end do
            end do
        case ('coarse')
            do j=1,nj
                do i=1,ni
                    if(i*j<ni*nj) then
                        write(100,140) t%coarseGrid(i,j,1),','
                    else
                        write(100,140) t%coarseGrid(i,j,1)
                    end if
                end do
            end do
    end select
    write(100,'(a)')           '                      ],'
    write(100,'(a)')           '          "ypoints" : ['
    !Write y points
    select case (name)
        case ('master')
            do j=1,nj
                do i=1,ni
                    if(i*j<ni*nj) then
                        write(100,140) t%masterGrid(i,j,2),','
                    else
                        write(100,140) t%masterGrid(i,j,2)
                    end if
                end do
            end do
        case ('medium')
            do j=1,nj
                do i=1,ni
                    if(i*j<ni*nj) then
                        write(100,140) t%mediumGrid(i,j,2),','
                    else
                        write(100,140) t%mediumGrid(i,j,2)
                    end if
                end do
            end do
        case ('coarse')
            do j=1,nj
                do i=1,ni
                    if(i*j<ni*nj) then
                        write(100,140) t%coarseGrid(i,j,2),','
                    else
                        write(100,140) t%coarseGrid(i,j,2)
                    end if
                end do
            end do
    end select
    write(100,'(a)')           '                      ]'
    write(100,'(a)')           '      }'
    write(100,'(a)')           '   }'
    write(100,'(a)')           '}'

    close(100)
    
end subroutine grid_save_web

!-----------------------------------------------------------------------------------------------------------
subroutine grid_save_surface(t)
    type(grid_t) :: t
    integer :: i
    do i=1,t%nSurfaceNodes
        write(*,*) t%surfacePoints(i,:)
    end do

end subroutine grid_save_surface

!-----------------------------------------------------------------------------------------------------------
subroutine grid_write(t)
    type(grid_t) :: t
    integer :: i,j

    write(*,*) 'Writing Grid to Screen'
    do j=1,t%nj
        do i=1,t%ni
            write(*,*) t%organicGrid(i,j,:)
        end do
    end do
    
end subroutine grid_write

!-----------------------------------------------------------------------------------------------------------
subroutine json_check()
    if(json_failed()) then
        call print_json_error_message()
        STOP
    end if
end subroutine json_check

!-----------------------------------------------------------------------------------------------------------
subroutine print_json_error_message()
    implicit none
    character(len=:),allocatable :: error_msg
    logical :: status_ok

    !get error message:
    call json_check_for_errors(status_ok, error_msg)

    !print it if there is one:
    if (.not. status_ok) then
        write(*,'(A)') error_msg
        deallocate(error_msg)
        call json_clear_exceptions()
    end if

end subroutine print_json_error_message

!-----------------------------------------------------------------------------------------------------------
subroutine grid_create_deflected_airfoil(t,afdata,lepercent,afdataFlapOut,XiSeg)
    type(grid_t) :: t
    type(dataset_t) :: afdata
    type(dataset_t) :: flapBott, flapTop
    type(dataset_t) :: afdataFlapOut
    real :: lepercent

    real :: myXi
    real :: splitPercent(2), foil(2), foil2(2), newcutpercent(2)
    integer :: i, j
    real, dimension(:,:), allocatable :: flapTemp, afdataFlap
    integer :: pointPerPercent = 2
    integer :: flapDataSize, afdataFlapSize
    real :: XiTempCent, x_hinge, y_hinge, de_rad
    real :: r_round, m_round, th, th1, th2, x_round, y_round, b_round
    real :: startStopPercent(3,2) ![[botflapStart,Botflapstop][foilStart,foilStop][topFlapStart,topFlapStop]]
    integer :: hinge_pts = 31
    real :: dataout(31,2), dataout2(31,2)
    real :: XiSeg(4)
    real :: mindist, xilength
    integer :: mincoord


    x_hinge = 1.0 - t%flapFraction
    y_hinge = t%flapHingeVertPos
    de_rad = t%flapDeflection*pi/180.0

    call grid_find_MinDistance(afdata,x_hinge,y_hinge,0.0,1,splitPercent(1))

    call grid_find_MinDistance(afdata,x_hinge,y_hinge,1.0,0,splitPercent(2))


    flapDataSize = nint(splitPercent(1)*100)*pointPerPercent
    allocate(flapTemp(flapDataSize, 2))

    ! Deflection Bott surface of flap
    XiTempCent = 0.0
    do i=1,flapDataSize
        call ds_cubic_interpolate(afdata,XiTempCent*afdata%Xi(afdata%datasize),0,foil)
        flapTemp(i,1) = x_hinge + (foil(1)-x_hinge)*cos(-de_rad)- &
                        & (foil(2)-y_hinge)*sin(-de_rad)
        flapTemp(i,2) = y_hinge + (foil(1)-x_hinge)*sin(-de_rad)+ &
                        & (foil(2)-y_hinge)*cos(-de_rad)
        XiTempCent = XiTempCent + splitPercent(1)/(flapDataSize-1)
    end do

    call ds_create_from_data(flapBott,flapDataSize,2,flapTemp)
    call ds_cubic_setup(flapBott,0,2,0.0,2,0.0)

    ! Deflection Top surface of flap
    XiTempCent = splitPercent(2)
    do i=1,flapDataSize
        call ds_cubic_interpolate(afdata,XiTempCent*afdata%Xi(afdata%datasize),0,foil)
        flapTemp(i,1) = x_hinge + (foil(1)-x_hinge)*cos(-de_rad)- &
                        & (foil(2)-y_hinge)*sin(-de_rad)
        flapTemp(i,2) = y_hinge + (foil(1)-x_hinge)*sin(-de_rad)+ &
                        & (foil(2)-y_hinge)*cos(-de_rad)
        XiTempCent = XiTempCent + (1.0 - splitPercent(2))/(flapDataSize-1)
    end do
    
    call ds_create_from_data(flapTop,flapDataSize,2,flapTemp(:,:))
    call ds_cubic_setup(flapTop,0,2,0.0,2,0.0)
    
    
    call ds_cubic_interpolate(afdata,splitPercent(1)*afdata%Xi(afdata%datasize),0,foil) !Bottom
    call ds_cubic_interpolate(afdata,splitPercent(2)*afdata%Xi(afdata%datasize),0,foil2) !Top

    startStopPercent(1,1) = 0.0
    startStopPercent(1,2) = 1.0
    startStopPercent(2,1) = splitPercent(1)
    startStopPercent(2,2) = splitPercent(2)
    startStopPercent(3,1) = 0.0
    startStopPercent(3,2) = 1.0

    ! Trimming and filling
    if( (t%flapDeflection >= 0.0 .and. y_hinge > foil2(2)) .or. &
            & (t%flapDeflection < 0.0 .and. y_hinge < foil(2)) ) then
        call grid_flap_trim(t,afdata,flapBott,1,splitPercent(1),x_hinge,y_hinge,dataout,newcutpercent)

        startStopPercent(1,2) = newcutpercent(1)
        startStopPercent(2,1) = newcutpercent(2)
        call grid_flap_trim(t,afdata,flaptop,0,splitPercent(2),x_hinge,y_hinge,dataout2,newcutpercent)
        startStopPercent(2,2) = newcutpercent(2)
        startStopPercent(3,1) = newcutpercent(1)

    else if (t%flapDeflection >= 0.0 .and. y_hinge < foil2(2) .and. y_hinge > foil(2)) then

        call grid_flap_fill(t,afdata,flaptop,0,splitPercent(2),x_hinge,y_hinge,dataout2)
        call grid_flap_trim(t,afdata,flapBott,1,splitPercent(1),x_hinge,y_hinge,dataout,newcutpercent)

        startStopPercent(1,2) = newcutpercent(1)
        startStopPercent(2,1) = newcutpercent(2)
        

    else if (t%flapDeflection <= 0.0 .and. y_hinge < foil2(2) .and. y_hinge > foil(2)) then

        call grid_flap_trim(t,afdata,flaptop,0,splitPercent(2),x_hinge,y_hinge,dataout2,newcutpercent)
        call grid_flap_fill(t,afdata,flapBott,1,splitPercent(1),x_hinge,y_hinge,dataout)
        startStopPercent(2,2) = newcutpercent(2)
        startStopPercent(3,1) = newcutpercent(1)
        
    else if ( (t%flapDeflection >= 0.0 .and. y_hinge < foil(2)) .or. &
            & (t%flapDeflection < 0.0 .and. y_hinge > foil2(2)) ) then
        call grid_flap_fill(t,afdata,flapBott,1,splitPercent(1),x_hinge,y_hinge,dataout)
        call grid_flap_fill(t,afdata,flaptop,0,splitPercent(2),x_hinge,y_hinge,dataout2)  

    else    
        write(*,*)'Bad new for radius'
        Stop

    end if
!
!    open(unit=13,file='dataplot.txt',status='replace',action='write')
!
!    write(13,*) '% toplabel   = "Circle"'
!    write(13,*) '% xlabel     = "x"'
!    write(13,*) '% ylabel     = "y"'
!    write(13,*) '% grid       = False'
!    write(13,*) '% equalscale = True'
!    write(13,*) '% axisscale  = False'
!    write(13,*) '% eyepos.x   = 0.50'
!    write(13,*) '% eyepos.y   = 0.75'
!    write(13,*) '% eyepos.z   = 0.25'
!    write(13,*) '% dlinecolor = 6'
!    write(13,*)
!
!    write(13,*) '% linecolor = 3'
!    myxi = splitPercent(1)
!    do i = 1,afdata%datasize
!        call ds_cubic_interpolate(afdata,myxi*afdata%Xi(afdata%datasize),0,foil)
!        write(13,'(f21.15,f21.15)') foil(1),foil(2)
!        myxi =  myxi  + (splitPercent(2)-splitPercent(1))/real(afdata%datasize-1)
!    end do
!
!
!    write(13,*)
!    write(13,*) '% linecolor = 2'
!        do i = 1,flapBott%datasize
!        myxi = real(i-1)/real(flapBott%datasize-1)
!!        write(*,*) myxi
!        call ds_cubic_interpolate(flapBott, myxi * flapBott%Xi(flapBott%datasize),0,foil)
!        write(13,'(f21.15,f21.15)') foil(1),foil(2)
!    end do
!
!    write(13,*)
!    write(13,*) '% linecolor = 2'
!    do i = 1,flapTop%datasize
!        myxi = real(i-1)/real(flapTop%datasize-1)
!        !        write(*,*) myxi
!        call ds_cubic_interpolate(flapTop, myxi * flapTop%Xi(flapTop%datasize),0,foil)
!        write(13,'(f21.15,f21.15)') foil(1),foil(2)
!    end do
!
!    write(13,*)
!    write(13,*) '% linecolor = 4'
!    do i = 1,30
!        write(13,'(f21.15,f21.15)') dataout(i,1),dataout(i,2)
!    end do
!
!    write(13,*)
!    write(13,*) '% linecolor = 4'
!    do i = 1,30
!        write(13,'(f21.15,f21.15)') dataout2(i,1),dataout2(i,2)
!    end do
!    CALL EXECUTE_COMMAND_LINE('plotmtv dataplot.txt')
!    close(13)
!


        
    afdataFlapSize = flapBott%datasize + afdata%datasize + &
                        & flapTop%datasize + size(dataout(:,2)) + size(dataout2(:,2))

    allocate(afdataFlap(afdataFlapSize,2))

    j = 1
    ! Flap Bottom
    myxi = startStopPercent(1,1)
    do i = 1,flapBott%datasize
        call ds_cubic_interpolate(flapBott, myxi * flapBott%Xi(flapBott%datasize),0,foil)
        afdataFlap(j,:) =  foil(:)
        myxi =  myxi  + (startStopPercent(1,2)-startStopPercent(1,1))/real(flapBott%datasize-1)
        j = j + 1
    end do

    ! Bottom Radius
    do i = 1,hinge_pts
        afdataFlap(j,:) = dataout(i,:)
        j = j + 1
    end do

    ! Foil
    myxi = startStopPercent(2,1)
    do i = 1,afdata%datasize
        call ds_cubic_interpolate(afdata,myxi*afdata%Xi(afdata%datasize),0,foil)
        afdataFlap(j,:) =  foil(:)
        myxi =  myxi  + (startStopPercent(2,2)-startStopPercent(2,1))/real(afdata%datasize-1)
        j = j + 1
    end do 

    ! Top Radius
    do i = 1,hinge_pts
        afdataFlap(j,:) = dataout2(i,:)
        j = j + 1
    end do

    ! Flap Top
    myxi = startStopPercent(3,1)
    do i = 1,flapTop%datasize
        call ds_cubic_interpolate(flapTop, myxi * flapTop%Xi(flapTop%datasize),0,foil)
        afdataFlap(j,:) =  foil(:)
        myxi =  myxi  + (startStopPercent(3,2)-startStopPercent(3,1))/real(flaptop%datasize-1)
        j = j + 1
    end do


    call ds_create_from_data(afdataFlapOut,afdataFlapSize,2,afdataFlap(:,:))
    call ds_cubic_setup(afdataFlapOut,0,2,0.0,2,0.0)
!    call ds_print_data(afdataFlapOut)



    !Find Leading-Edge
    mindist = math_length(2,[0.0,0.0],afdataFlapOut%RawData(1,:))
    mincoord = 1
    do i=1,afdataFlapOut%datasize
        if(math_length(2,[0.0,0.0],afdataFlapOut%RawData(i,:)) < mindist) then
            mindist = math_length(2,[0.0,0.0],afdataFlapOut%RawData(i,:))
            mincoord = i
        end if
    end do
    xilength = afdataFlapOut%Xi(afdataFlapOut%datasize)
    lepercent = afdataFlapOut%Xi(mincoord)/xilength

    XiSeg(1) = afdataflapout%Xi(flapBott%datasize + (hinge_pts+1)/2) ! Center of bottom radius
    XiSeg(2) = lepercent * afdataflapout%Xi(afdataflapout%datasize) ! Leading edge
    XiSeg(3) = afdataflapout%Xi(flapBott%datasize + hinge_pts + afdata%datasize + (hinge_pts+1)/2) ! Center of Top radius
    XiSeg(4) = afdataflapout%Xi(afdataFlapSize) ! Trailing edge approched from top

    write(*,*) flapBott%datasize + (hinge_pts+1)/2
    write(*,*) flapBott%datasize + hinge_pts + (afdata%datasize+1)/2
    write(*,*) flapBott%datasize + hinge_pts + afdata%datasize + (hinge_pts+1)/2
    write(*,*) afdataFlapSize
    

!    open(unit=13,file='dataplot.txt',status='replace',action='write')
!
!    write(13,*) '% toplabel   = "Circle"'
!    write(13,*) '% xlabel     = "x"'
!    write(13,*) '% ylabel     = "y"'
!    write(13,*) '% grid       = False'
!    write(13,*) '% equalscale = True'
!    write(13,*) '% axisscale  = False'
!    write(13,*) '% eyepos.x   = 0.50'
!    write(13,*) '% eyepos.y   = 0.75'
!    write(13,*) '% eyepos.z   = 0.25'
!    write(13,*) '% dlinecolor = 6'
!    write(13,*)
!
!    write(13,*) '% linecolor = 3'
!    do i = 1,afdataFlapSize
!        write(13,'(f21.15,f21.15)') afdataFlap(i,:)
!    end do
!
!    CALL EXECUTE_COMMAND_LINE('plotmtv dataplot.txt')
!    close(13)


    deallocate(flapTemp)
    deallocate(afdataFlap)
end subroutine grid_create_deflected_airfoil

!-----------------------------------------------------------------------------------------------------------
subroutine  grid_flap_trim(t,afdata,flapdata,topOrBott, & 
                & splitPercent,x_hinge,y_hinge,dataout,splitPercent_flap)

    type(grid_t) :: t
    type(dataset_t) :: afdata
    type(dataset_t) :: flapdata
    integer :: topOrBott  ! 1 = Bott, 0 = Top
    real :: splitPercent, x_hinge, y_hinge

    real :: dataout(31,2)
    real :: th_start, th_stop, th, de_rad, r
    real :: foil(2), foil2(2)
    integer :: npts = 31, i

    real :: x_round, y_round, r_round, m_round, b_round
    real :: th1, th2
    real :: splitPercent_flap(2)
    
    real :: trim_radius = .05

    de_rad = t%flapDeflection * pi/180.0

    call ds_cubic_interpolate(afdata,splitPercent*afdata%Xi(afdata%datasize),0,foil)

    if (topOrBott .eq. 1) then
        call ds_cubic_interpolate(flapdata,1.0*flapdata%xi(flapdata%datasize),0,foil2)
    else
        call ds_cubic_interpolate(flapdata,0.0*flapdata%xi(flapdata%datasize),0,foil2)
    end if

    r_round = sqrt((y_hinge-foil(2))**2+(x_hinge-foil(1))**2)+trim_radius

    call grid_find_angle(x_hinge,y_hinge,foil(1),foil(2),th1)
    call grid_find_angle(x_hinge,y_hinge,foil2(1),foil2(2),th2)

    th1 = (th1 + th2)/2.0
    m_round = tan(th1)

    y_round = y_hinge + sin(th1) * r_round
    b_round = y_hinge - m_round * x_hinge
    x_round = (y_round-b_round)/m_round
!    write(*,*) 'hi',th1*180.0/pi,x_round, y_round

    if (topOrBott .eq. 1) then
        call grid_find_MinDistance(flapdata,x_round,y_round,0.0,1,splitPercent_flap(1))
        call grid_find_MinDistance(afdata,x_round,y_round,0.0,1,splitPercent_flap(2))
    else
        call grid_find_MinDistance(flapdata,x_round,y_round,1.0,0,splitPercent_flap(1))
        call grid_find_MinDistance(afdata,x_round,y_round,1.0,0,splitPercent_flap(2))
    end if
    
    call ds_cubic_interpolate(afdata,splitPercent_flap(2)*afdata%Xi(afdata%datasize),0,foil)
    call grid_find_angle(x_round,y_round,foil(1),foil(2),th)

    r = sqrt((y_round-foil(2))**2+(x_round-foil(1))**2)


    if (topOrBott .eq. 1) then !Bott
        th_start = th  - de_rad
        th_stop = th
    else  ! Top
        th_start = th  
        th_stop = th - de_rad 
    end if
    th = th_start    
    do i = 1,npts
        th = th - (th_start-th_stop)/(real(npts+1))
        dataout(i,1) = r*cos(th) + x_round
        dataout(i,2) = r*sin(th) + y_round
    end do

    
end subroutine grid_flap_trim

!-----------------------------------------------------------------------------------------------------------
subroutine  grid_flap_fill(t,afdata,flapdata,topOrBott,splitPercent,x_hinge,y_hinge,dataout)

    type(grid_t) :: t
    type(dataset_t) :: afdata
    type(dataset_t) :: flapdata
    integer :: topOrBott  ! 1 = trim/Bott, 0 = radius/Top
    real :: x_hinge, y_hinge, splitPercent

    real :: dataout(31,2)
    real :: th_start, th_stop, th, de_rad, r
    real :: foil(2), foil2(2)
    integer :: npts = 31, i

    de_rad = t%flapDeflection * pi/180.0

    call ds_cubic_interpolate(afdata,splitPercent*afdata%Xi(afdata%datasize),0,foil)
    call grid_find_angle(x_hinge,y_hinge,foil(1),foil(2),th)
    r = sqrt((y_hinge-foil(2))**2+(x_hinge-foil(1))**2)

    if (topOrBott .eq. 1) then !Bott
        th_start = th - de_rad
        th_stop = th
    else  ! Top
        th_start = th   
        th_stop = th - de_rad
    end if

    th = th_start
    do i = 1,npts
        th = th - (th_start-th_stop)/(real(npts+1))
        dataout(i,1) = r*cos(th) + x_hinge
        dataout(i,2) = r*sin(th) + y_hinge
    end do

end subroutine grid_flap_fill

!-----------------------------------------------------------------------------------------------------------
subroutine grid_find_angle(xc,yc,xp,yp,th)
    real :: xc,yc, xp,yp, th

    th = atan((yc-yp)/(xc-xp))
    if (xp < xc .and. yp < yc) then
        th = th - pi
    elseif (xp < xc .and. yp > yc) then
        th = th + pi
    else    
    end if
end subroutine grid_find_angle

!-----------------------------------------------------------------------------------------------------------
subroutine grid_find_MinDistance(dataset,xo,yo,startPercent,direction,splitPercent)
    type(dataset_t) :: dataset
    real :: xo, yo, startPercent, splitPercent
    integer :: direction ! 1 => +, 0 => -

    real :: intersection(2)
    real :: r, r_1
    real :: convergeInter = 1.0E-15
    real :: step   
    integer :: i

    if (direction .eq. 1) then
        step = 0.05
    else
        step = -0.05
    end if
    
    call ds_cubic_interpolate(dataset,startPercent*dataset%Xi(dataset%datasize),0,intersection)
    r_1 = sqrt((yo-intersection(2))**2+(xo-intersection(1))**2)
    splitPercent = startPercent + step

    i = 0
    do !i = 1, 40
        call ds_cubic_interpolate(dataset,splitPercent*dataset%Xi(dataset%datasize),0,intersection)
        r  = sqrt((yo-intersection(2))**2+(xo-intersection(1))**2)
!        write(*,*) abs(r1-r2)
        if (abs(r - r_1) < convergeInter) then
            Exit
        elseif (r > r_1) then
            step = -step /2.0
        else
!            write(*,*) 'step'
        end if

        splitpercent = splitPercent + step
        r_1 = r
        i = i + 1
    
        if (i > 500) then
            write(*,*)  'I''m Stuck in grid_find_MinDistance: iter =', i
        end if
    end do

end subroutine grid_find_MinDistance

end module grid_m
