        ! for setting up problem
        module Input
        use equations
        use line
        use regions
        use EulerDeclarations
        implicit none
        
        ! Named initial conditions
        integer, parameter :: Gaussian = 1, SquareWave = 2, BurgersShock = 3, explosion = 1
        
        contains
        
        ! read in options, time, space, and equation information
        subroutine InitDomain()
        double precision alpha,beta,radius
        double precision, allocatable :: x(:),WL(:),WR(:),x0(:)
        integer i,j,k,IC
        logical ex

        NAMELIST/ProblemData/ starttime,finaltime,frames,xL,xR,xcells,yL,yR,ycells,zL,zR,zcells,BC,IC,CFL,cov,gamma,tol,maxiter,RiemannSolver,Q,alpha,beta,advectionConst,Eq,nDim,num_threads,lthreading,FluxScheme
        
        ! Make sure problem.data exists. If it does, read in. If not, output error and quit.
        inquire(file='problem.data', exist=ex)
        if (ex) then
         OPEN(UNIT=30, FILE='problem.data', STATUS="OLD")
         READ(30,NML=ProblemData)                               ! read into variables in namelist
         CLOSE(30)
        else
         print *, 'Missing problem.data. Exiting.'
         stop           ! exit if missing
        end if
        
        allocate(x(nDim),x0(nDim),dx(nDim),imom(nDim),iu(nDim))
        
        ! Set array locations for multidimensional problems and density
        ix = 1
        iy = 2
        iz = 3
        irho = 1
        
        ! Set number of variables, location of each (global)
        select case (progNum)
         case (Euler)
          select case (nDim)
           case (1)
            NrVar = 3
            ! Primitives
            iu(ix) = 2
            ip = 3
            ! Conserved
            imom(ix) = 2
            iE = 3
           case (2)
            NrVar = 4
            ! Primitives
            iu(ix) = 2
            iu(iy) = 3
            ip = 4
            ! Conserved
            imom(ix) = 2
            imom(iy) = 3
            iE = 4
           case (3)
            NrVar = 5
            ! Primitives
            iu(ix) = 2
            iu(iy) = 3
            iu(iz) = 4
            ip = 5
            ! Conserved
            imom(ix) = 2
            imom(iy) = 3
            imom(iz) = 4
            iE = 5
          end select
         case (GodunovSingle)
          NrVar = 1
         case (ExactProg)
          NrVar = 3
          ! Primitives
          iu(ix) = 2
          ip = 3
        end select
        
        allocate(WL(NrVar),WR(NrVar))
        
        ! Initialize time and frame and set intercell distances
        time = starttime
        frame = 0
        dx(ix) = (xR-xL)/xcells
        if (nDim > 1) dx(iy) = (yR-yL)/ycells
        if (nDim > 2) dx(iz) = (zR-zL)/zcells
        
        ! Initialize based on selected test (see Toro, Numerical Methods)
        select case (progNum)
               
         case (Euler)
          select case (nDim)
           
           ! 1D case
           case (1)
            ! Allocate domain arrays for conserved variables and intercell flux
            allocate(ProblemRegion%U(-1:xcells+2,NrVar),ProblemRegion%flux(0:xcells,NrVar))
                    
            select case (IC)
             case (1)
              WL(irho) = 1.0
              WR(irho) = 0.125
              WL(iu) = 0.75
              WR(iu) = 0.0
              WL(ip) = 1.0
              WR(ip) = 0.1
              x0(ix) = 0.3
             case (2)
              WL(irho) = 1.0
              WR(irho) = 1.0
              WL(iu) = -2.0
              WR(iu) = 2.0
              WL(ip) = 0.4
              WR(ip) = 0.4
              x0(ix) = 0.5
             case (3)
              WL(irho) = 1.0
              WR(irho) = 1.0
              WL(iu) = 0.0
              WR(iu) = 0.0
              WL(ip) = 1000.0
              WR(ip) = 0.01
              x0(ix) = 0.5
             case (4)
              WL(irho) = 5.99924
              WR(irho) = 5.99242
              WL(iu) = 19.975
              WR(iu) = -6.19633
              WL(ip) = 460.894
              WR(ip) = 46.0950
              x0(ix) = 0.4
             case (5)
              WL(irho) = 1.0
              WR(irho) = 1.0
              WL(iu) = -19.59745
              WR(iu) = -19.59745
              WL(ip) = 1000.0
              WR(ip) = 0.01
              x0(ix) = 0.8
             case (6)
              WL(irho) = 1.4
              WR(irho) = 1.0
              WL(iu) = 0.0
              WR(iu) = 0.0
              WL(ip) = 1.0
              WR(ip) = 1.0
              x0(ix) = 0.5
             case (7)
              WL(irho) = 1.4
              WR(irho) = 1.0
              WL(iu) = 0.1
              WR(iu) = 0.1
              WL(ip) = 1.0
              WR(ip) = 1.0
              x0(ix) = 0.5
            end select
        
            ! Initialize domain using set initial conditions (conserved variables)
            do i = 1, xcells
             x(ix) = (i-0.5)*dx(ix) + xL
             if (x(ix) <= x0(ix)) then
              ProblemRegion%U(i,:) = PrimToCons(WL)
             else
              ProblemRegion%U(i,:) = PrimToCons(WR)
             end if
            end do
           
           ! 2D case
           case (2)
            allocate(ProblemRegion2D%U(-1:xcells+2,-1:ycells+2,NrVar),ProblemRegion%U(-1:max(xcells,ycells)+2,NrVar),ProblemRegion%flux(0:max(xcells,ycells),NrVar))
          
            select case (IC)
             case (explosion)
              ! "Left" is inside circular boundary, "Right" is outside circular boundary
              WL(irho) = 1.0
              WL(iu(ix)) = 0.0
              WL(iu(iy)) = 0.0
              WL(ip) = 1.0
              WR(irho) = 0.125
              WR(iu(ix)) = 0.0
              WR(iu(iy)) = 0.0
              WR(ip) = 0.1
              x0 = (/ 0, 0 /)
              radius = 0.4
              do i = 1, xcells
               x(ix) = (i-0.5)*dx(ix) + xL
               do j = 1, ycells
                x(iy) = (j-0.5)*dx(iy) + yL
                if (sum((x-x0)**2) < radius**2) then
                 ProblemRegion2D%U(i,j,:) = PrimToCons(WL)
                 print *, ProblemRegion2D%U(i,j,:)
                else
                 ProblemRegion2D%U(i,j,:) = PrimToCons(WR)
                 print *, ProblemRegion2D%U(i,j,:)
                end if
               end do
              end do
            end select
           
           ! 3D case          
           case (3)
            allocate(ProblemRegion3D%U(-1:xcells+2,-1:ycells+2,-1:zcells+2,NrVar),ProblemRegion%U(-1:max(xcells,ycells,zcells)+2,NrVar),ProblemRegion%flux(0:max(xcells,ycells,zcells),NrVar))
            
            select case (IC)
             case (explosion)
              ! "Left" is inside spherical boundary, "Right" is outside spherical boundary
              WL(irho) = 1.0
              WL(iu(ix)) = 0.0
              WL(iu(iy)) = 0.0
              WL(iu(iz)) = 0.0
              WL(ip) = 1.0
              WR(irho) = 0.125
              WR(iu(ix)) = 0.0
              WR(iu(iy)) = 0.0
              WR(iu(iz)) = 0.0
              WR(ip) = 0.1
              x0 = (/ 0, 0, 0 /)
              radius = 0.4
              do i = 1, xcells
               x(ix) = (i-0.5)*dx(ix) + xL
               do j = 1, ycells
                x(iy) = (j-0.5)*dx(iy) + yL
                do k = 1, zcells
                 x(iz) = (k-0.5)*dx(iz) + zL
                 if (sum((x-x0)**2) < radius**2) then
                  ProblemRegion3D%U(i,j,k,:) = PrimToCons(WL)
                 else
                  ProblemRegion3D%U(i,j,k,:) = PrimToCons(WR)
                 end if
                end do
               end do
              end do
            end select
          
          end select
          
         case (GodunovSingle)
          allocate(ProblemRegion%u(0:xcells+1,NrVar),ProblemRegion%flux(0:xcells,NrVar))
          
          ! Initialize domain based on selected initial conditions
          select case (IC)
           case (Gaussian)
            do i = 1, xcells
             ProblemRegion%u(i,NrVar) = alpha*exp(-1.*beta*((i-0.5)*dx(ix)+xL)**2)
            end do
           case (SquareWave)
            do i = 1, xcells
            x(ix) = (i-0.5)*dx(ix) + xL
             if (x(ix) <= 0.3) then
              ProblemRegion%u(i,NrVar) = 0
             else if (x(ix) >= 0.7) then
              ProblemRegion%u(i,NrVar) = 0
             else
              ProblemRegion%u(i,NrVar) = 1
             end if
            end do
           case (BurgersShock)
            do i = 1, xcells
            x(ix) = (i-1)*dx(ix) + xL
             if (x(1) <= 0.5) then
              ProblemRegion%u(i,NrVar) = -0.5
             else if (x(ix) >= 1.) then
              ProblemRegion%u(i,NrVar) = 0
             else
              ProblemRegion%u(i,NrVar) = 1
             end if
            end do
          end select
 
         case (ExactProg)
          select case (IC)
           case (1)
            WL(irho) = 1.0
            WR(irho) = 0.125
            WL(iu) = 0.0
            WR(iu) = 0.0
            WL(ip) = 1.0
            WR(ip) = 0.1
           case (2)
            WL(irho) = 1.0
            WR(irho) = 1.0
            WL(iu) = -2.0
            WR(iu) = 2.0
            WL(ip) = 0.4
            WR(ip) = 0.4
           case (3)
            WL(irho) = 1.0
            WR(irho) = 1.0
            WL(iu) = 0.0
            WR(iu) = 0.0
            WL(ip) = 1000.0
            WR(ip) = 0.01
           case (4)
            WL(irho) = 1.0
            WR(irho) = 1.0
            WL(iu) = 0.0
            WR(iu) = 0.0
            WL(ip) = 0.01
            WR(ip) = 100.0
           case (5)
            WL(irho) = 5.99924
            WR(irho) = 5.99242
            WL(iu) = 19.5975
            WR(iu) = -6.19633
            WL(ip) = 460.894
            WR(ip) = 46.0950
          end select
          
          L = region(WL)
          R = region(WR)
        end select
        
        ! Clean up
        if (allocated(x)) deallocate(x)
        if (allocated(x0)) deallocate(x0)
          
        return
        
        end subroutine
        
        end module
