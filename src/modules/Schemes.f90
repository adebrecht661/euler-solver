        module schemes
        use line
        use regions
        use equations
        use EulerDeclarations
        use RiemannSolutions
        implicit none
        
        integer, parameter :: advection = 1, inviscidBurgers = 2
        
        contains
        
        ! Calculate intercell flux using specified equation
        subroutine GodunovUpwindFlux(ProblemRegion)
        integer i
        double precision uL,uR,S0,SL,SR
        double precision U(NrVar)
        type(region) L,R
        type(OneDdomain), intent(inout) :: ProblemRegion
     
        select case (progNum)
         
         ! Euler solver
         case (Euler)
          allocate(L%W(NrVar),R%W(NrVar))
          
          ! Loop through domain using equation set in problem.data
          do i = 0, cells
           
           U = ProblemRegion%U(i,:)
           L = region(ConsToPrim(U))
           if (L%c .ne. L%c) then
            print *, 'NaN in sound speed at time ', time, ' in cell ', i,' in CalcFlux.'
            stop
           end if
           U = ProblemRegion%U(i+1,:)
           R = region(ConsToPrim(U))
           
           ! Some solvers approximate flux rather than fluid state
           select case (RiemannSolver)
            case (HLLC)
             call WaveSpeed(L,R,SL,SR,S0)
             ProblemRegion%flux(i,:) = HLLCFlux(L,R,S0,SL,SR)         ! Flux estimate along t axis
            case (HLL)
             call WaveSpeed(L,R,SL,SR)
             ProblemRegion%flux(i,:) = HLLFlux(L,R,SL,SR)
            case (Roe)
             ProblemRegion%flux(i,:) = RoeFlux(L,R)
            ! If solver approximates fluid state, use it to calculate flux
            case default
             ProblemRegion%flux(i,:) = EulerFlux(RiemannSoln(L,R))
           end select
          end do
          
          ! Clean up
          deallocate(L%W,R%W)
          
         ! Flux for single variable Godunov solver 
         case (GodunovSingle)
          do i = 0, xcells
           uL = ProblemRegion%u(i,NrVar)
           uR = ProblemRegion%u(i+1,NrVar)   
           ! Use selected equation
           select case (Eq)
            case (advection)
             ProblemRegion%flux(i,NrVar) = LinearAdvectionFlux(RiemannSoln(uL,uR))
            case (inviscidBurgers)
             ProblemRegion%flux(i,NrVar) = BurgersFlux(RiemannSoln(uL,uR))
           end select
          end do
          
        end select
          
        return
        
        end subroutine
           
        subroutine WAFFlux(ProblemRegion)
        integer i,k
        double precision S0,SL,SR
        double precision U(NrVar), WAFState(NrVar), Courant(3), temp(NrVar)
        type(region) L,R
        type(OneDdomain), intent(inout) :: ProblemRegion
        
        do i = 0,cells
         U = ProblemRegion%U(i,:)
         L = region(ConsToPrim(U))
         U = ProblemRegion%U(i+1,:)
         R = region(ConsToPrim(U))
         
         call WaveSpeed(L,R,SL,SR,S0)
         
         Courant(1) = dt*SL/dl
         Courant(2) = dt*S0/dl
         Courant(3) = dt*SR/dl
              
         select case (RiemannSolver)
          case (HLLC)
           temp = 0d0
           do k = 1,3
            temp = temp + (sign(5d-1,Courant(k)) - sign(5d-1,-Courant(k)))*WAFLim(Courant(k),DensRatio(L,R,k,Courant(k),i))*(HLLCFlux(L,R,S0,SL,SR,k+1) - HLLCFlux(L,R,S0,SL,SR,k))
           end do
           print *, temp
           ProblemRegion%flux(i,:) = 0.5*(EulerFlux(L%W) + EulerFlux(R%W) - temp)
          
          case default
           temp = 0d0
           do k = 1,3
            temp = temp + (sign(5d-1,Courant(k)) - sign(5d-1,-Courant(k)))*WAFLim(Courant(k),DensRatio(L,R,k,Courant(k),i))*(RiemannSoln(L,R,k+1) - RiemannSoln(L,R,k))
           end do
           WAFState = 0.5*(L%W + R%W - temp)
           ProblemRegion%flux(i,:) = EulerFlux(WAFState)
         end select
        end do
        
        return
        
        end subroutine
        
        double precision function DensRatio(L,R,k,Courant,cell)
        double precision, intent(in) :: Courant
        integer, intent(in) :: k,cell
        type(region), intent(in) :: L,R
        double precision denom
        double precision U(NrVar)
        type(region) Lstar,Rstar,LL,RR,LLstar,RRstar
        
        denom = 0d0
        allocate(Lstar%W(NrVar),Rstar%W(NrVar))
        call StarRegion(L,R,Lstar,Rstar)
        
        select case (k)
         case (1)
          denom = Lstar%W(irho) - L%W(irho)
         case (2)
          denom = Rstar%W(irho) - Lstar%W(irho)
         case (3)
          denom = R%W(irho) - Rstar%W(irho)
        end select
        
        if (Courant < 0) then
         allocate(LL%W(NrVar),LLstar%W(NrVar))
         U = ProblemRegion%U(cell-1,:)
         LL = region(ConsToPrim(U))
         call StarRegion(LL,L,LLstar,Lstar)
         select case (k)
          case (1)
           DensRatio = (LLstar%W(irho) - LL%W(irho))/denom
          case (2)
           DensRatio = (Lstar%W(irho) - LLstar%W(irho))/denom
          case (3)
           DensRatio = (L%W(irho) - Lstar%W(irho))/denom
         end select
         deallocate(LL%W,LLstar%W)
        else
         allocate(RR%W(NrVar),RRstar%W(NrVar))
         U = ProblemRegion%U(cell+2,:)
         RR = region(ConsToPrim(U))
         call StarRegion(R,RR,Rstar,RRstar)
         select case (k)
          case (1)
           DensRatio = (Rstar%W(irho) - R%W(irho))/denom
          case (2)
           DensRatio = (RRstar%W(irho) - Rstar%W(irho))/denom
          case (3)
           DensRatio = (RR%W(irho) - RRstar%W(irho))/denom
         end select
         deallocate(RR%W,RRstar%W)
        end if
        
        deallocate(Lstar%W,Rstar%W)
        
        end function
        
        double precision function WAFLim(Courant,DensRatio)
        double precision, intent(in) :: Courant, DensRatio
        
        if (DensRatio <= 0) then
         WAFLim = 1.
        else if (DensRatio <= 0.5) then
         WAFLim = 1. - 2.*(1. - abs(Courant))*DensRatio
        else if (DensRatio <= 1 .or. DensRatio /= DensRatio) then
         WAFLim = abs(Courant)
        else if (DensRatio <= 2) then
         WAFLim = 1. - (1. - abs(Courant))*DensRatio
        else
         WAFLim = 2.*abs(Courant) - 1.
        end if
        
        end function
        
        end module
