        module sweep
        use line
        use regions
        use equations
        use EulerDeclarations
        use schemes
        implicit none
        
        contains
        
        ! Calculate delta t from CFL and distance between cells
        subroutine TimeStep()
        double precision maxSpeed
        integer i,j,k
        double precision U(NrVar)
        type(region) cell
        
        ! Determine maximum speed in domain
        maxSpeed = 0.
        !!! With implemented boundary conditions, don't need to look at ghost cells - revisit if boundary conditions are ever not a simple equality or negation !!!
        select case (progNum)
         
         ! Euler solver
         case (Euler)
          ! 1D
          if (nDim == 1) then
           do i = 1, xcells
            U = ProblemRegion%U(i,:)
            cell = region(ConsToPrim(U))
            if (cell%c .ne. cell%c) then
             print *, 'NaN in sound speed at time ', time, ' in cell ', i,' in TimeStep.'
             print *, cell%W, cell%c
             stop
            end if
            maxSpeed = max(maxSpeed, abs(cell%W(iu(ix))) + cell%c)
           end do
          ! 2D
          else if (nDim == 2) then
           do i = 1, xcells
            do j = 1, ycells
             U = ProblemRegion2D%U(i,j,:)
             cell = region(ConsToPrim(U))
             if (cell%c .ne. cell%c) then
              print *, 'NaN in sound speed at time ', time, ' in cell ', i,' in TimeStep.'
              stop
             end if
             maxSpeed = max(maxSpeed, abs(cell%W(iu(ix))) + cell%c)
             maxSpeed = max(maxSpeed, abs(cell%W(iu(iy))) + cell%c)
            end do
           end do
          ! 3D
          else if (nDim == 3) then
           do i = 1, xcells
            do j = 1, ycells
             do k = 1, zcells
              U = ProblemRegion3D%U(i,j,k,:)
              cell = region(ConsToPrim(U))
              if (cell%c .ne. cell%c) then
               print *, 'NaN in sound speed at time ', time, ' in cell ', i,' in TimeStep.'
               stop
              end if
              maxSpeed = max(maxSpeed, abs(cell%W(iu(ix))) + cell%c)
              maxSpeed = max(maxSpeed, abs(cell%W(iu(iy))) + cell%c)
              maxSpeed = max(maxSpeed, abs(cell%W(iu(iz))) + cell%c)
             end do
            end do
           end do
          end if
         
         ! Single-variable Godunov solver 
         case (GodunovSingle)
          do i = 0, xcells+1
           maxSpeed = max(maxSpeed, abs(ProblemRegion%u(i,1)) )    ! (See Toro, Ch. 5.6)
          end do
        
        end select
        
        ! Print to console
        print *, "MaxSpeed = ", maxSpeed
        
        ! Calculate delta t
        dt = minval(CFL*dx/maxSpeed)
        
        ! If time step would pass the next output frame, recalculate
        if (time + dt > nexttime) dt = nexttime - time
        
        end subroutine
        
        ! Calculate intercell flux using specified equation
        subroutine CalcFlux(ProblemRegion)
        type(OneDdomain), intent(inout) :: ProblemRegion
     
        select case (FluxScheme)
         case (GodunovUpwind)
          call GodunovUpwindFlux(ProblemRegion)
         case (WAF)
          call WAFFlux(ProblemRegion)
        end select
        
        end subroutine
        
        ! Update velocity array 
        subroutine UpdateCells(ProblemRegion)
        integer i
        type(OneDdomain), intent(inout) :: ProblemRegion
        
        ! Loop through cells using previously calculated flux
        do i = 1, cells
         ProblemRegion%U(i,:) = ProblemRegion%U(i,:) + (dt/dl)*(ProblemRegion%flux(i-1,:) - ProblemRegion%flux(i,:))
         if (any(ProblemRegion%U(i,:) .ne. ProblemRegion%U(i,:))) then
          print *,' NaN at time ', time, 'in cell ',i, ' in UpdateCells.'
          stop
         end if
        end do
        
        end subroutine
        
        end module
