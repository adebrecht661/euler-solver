        ! Uses selected solver to solve or approximate Riemann problem solution
        module RiemannSolutions
        use regions
        use EulerDeclarations
        use ExactSolver
        use ApproximateSolver
        implicit none
        
        ! Selection parameters for Riemann solver
        integer, parameter :: exact = 1, TSRS = 2, ANRS = 3, HLL = 4, HLLC = 5, Roe = 6
        
        private rarefaction
        
        interface sample
         module procedure sampleWithVacuum
         module procedure sampleWithoutVacuum
        end interface
        
        interface RiemannSoln
         module procedure EulerRiemannSoln
         module procedure SingleRiemannSoln
        end interface
        
        contains
        
        ! Solve local Riemann problem RP(L,R) along t axis (x/t = 0)
        function EulerRiemannSoln(L,R,k)
        integer, intent(in), optional :: k
        type(region), intent(in) :: L,R
        double precision EulerRiemannSoln(NrVar)
        type(region) Lstar, Rstar
        
        ! Check for vacuum
        call check_vacuum(L,R)
        
        ! Vacuum present, don't need star region - exact solution is cheap
        if (vacuum /= 3) then
         if (FluxScheme == WAF) then
          print *, 'WAF flux scheme not implemented for vacuum states. Halting.'
          stop
         end if
         EulerRiemannSoln = sample(L,R,0d0)
        ! No vacuum, calculate star region and sample along t axis
        else
         ! Star region may be approximated
         call StarRegion(L,R,Lstar,Rstar)
         if (present(k)) then
          select case (k)
           case (1)
            EulerRiemannSoln = L%W
           case (2)
            EulerRiemannSoln = Lstar%W
           case (3)
            EulerRiemannSoln = Rstar%W
           case (4)
            EulerRiemannSoln = R%W
          end select
         else
          EulerRiemannSoln = sample(L,R,Lstar,Rstar,0d0)
         end if
         deallocate(Lstar%W,Rstar%W)
        end if
        
        end function
        
        ! Solve the local Riemann problem RP(uL,uR) along t axis, velocity only
        double precision function SingleRiemannSoln(uL,uR)
        double precision, intent(in) :: uL,uR
        double precision S
        
        if (uL > uR) then
         S = 0.5*(uL + uR)
         if (S <= 0) then
          SingleRiemannSoln = uR
         else
          SingleRiemannSoln = uL
         end if
        else
         if (uL >= 0) then
          SingleRiemannSoln = uL
         else if (uR <= 0) then
          SingleRiemannSoln = uR
         else
          SingleRiemannSoln = 0
         end if
        end if
        
        end function
        
        ! Solves for primitive variables in interaction region using selected method
        subroutine StarRegion(L,R,Lstar,Rstar)
        type(region), intent(in) :: L,R
        type(region), intent(out) :: Lstar,Rstar
        double precision p0,pmin,pmax
        
        allocate(Lstar%W(NrVar),Rstar%W(NrVar))
        
        ! Passively advected - no math, so no approximations
        if (nDim > 1) then
         Lstar%W(iu(iy)) = L%W(iu(iy))
         Rstar%W(iu(iy)) = R%W(iu(iy))
        end if
        if (nDim > 2) then
         Lstar%W(iu(iz)) = L%W(iu(iz))
         Rstar%W(iu(iz)) = R%W(iu(iz))
        end if
        
        select case (RiemannSolver)
         ! Exact Riemann solver
         case (exact)
          ! Use Newton-Raphson solver to find pressure in star region
          Lstar%W(ip) = exactPressure(L,R)
          ! Calculate velocity and densities in star regions
          Lstar%W(iu(ix)) = exactVelocity(Lstar%W(ip),L,R)
          Rstar%W(ip) = Lstar%W(ip)
          Rstar%W(iu(ix)) = Lstar%W(iu(ix))
          Lstar%W(irho) = exactDensity(Lstar,L)
          Rstar%W(irho) = exactDensity(Rstar,R)
         ! Two-shock approximation
         case (TSRS)
          p0 = max(PrimitiveVariablePressure(L,R),0d0)
          Lstar%W(ip) = TwoShockPressure(L,R,p0)
          Lstar%W(iu(ix)) = TwoShockVelocity(L,R,Lstar,p0)
          Rstar%W(ip) = Lstar%W(ip)
          Rstar%W(iu(ix)) = Lstar%W(iu(ix))
          Lstar%W(irho) = TwoShockDensity(L,Lstar)
          Rstar%W(irho) = TwoShockDensity(R,Rstar)
         ! Adaptive non-iterative method
         case (ANRS)
          p0 = PrimitiveVariablePressure(L,R)
          pmax = max(L%W(ip),R%W(ip))
          pmin = min(L%W(ip),R%W(ip))
          if (pmax/pmin < Q .and. pmin < p0 .and. pmax > p0) then
           Lstar%W(ip) = p0
           Lstar%W(iu(ix)) = PrimitiveVariableVelocity(L,R)
           Rstar%W(ip) = Lstar%W(ip)
           Rstar%W(iu(ix)) = Lstar%W(iu(ix))
           Lstar%W(irho) = PrimitiveVariableDensity(L,Lstar,R)
           Rstar%W(irho) = PrimitiveVariableDensity(R,Rstar,L)
          else if (p0 < pmin) then
           Lstar%W(ip) = TwoRarefactionPressure(L,R)
           Lstar%W(iu(ix)) = TwoRarefactionVelocity(L,R)
           Rstar%W(ip) = Lstar%W(ip)
           Rstar%W(iu(ix)) = Lstar%W(iu(ix))
           Lstar%W(irho) = TwoRarefactionDensity(L,Lstar)
           Rstar%W(irho) = TwoRarefactionDensity(R,Rstar)
          else
           Lstar%W(ip) = TwoShockPressure(L,R,p0)
           Lstar%W(iu(ix)) = TwoShockVelocity(L,R,Lstar,p0)
           Rstar%W(ip) = Lstar%W(ip)
           Rstar%W(iu(ix)) = Lstar%W(iu(ix))
           Lstar%W(irho) = TwoShockDensity(L,Lstar)
           Rstar%W(irho) = TwoShockDensity(R,Rstar)
          end if
         case default   ! if needed for other calculations, select ANRS
          p0 = PrimitiveVariablePressure(L,R)
          pmax = max(L%W(ip),R%W(ip))
          pmin = min(L%W(ip),R%W(ip))
          if (pmax/pmin < Q .and. pmin < p0 .and. pmax > p0) then
           Lstar%W(ip) = p0
           Lstar%W(iu(ix)) = PrimitiveVariableVelocity(L,R)
           Rstar%W(ip) = Lstar%W(ip)
           Rstar%W(iu(ix)) = Lstar%W(iu(ix))
           Lstar%W(irho) = PrimitiveVariableDensity(L,Lstar,R)
           Rstar%W(irho) = PrimitiveVariableDensity(R,Rstar,L)
          else if (p0 < pmin) then
           Lstar%W(ip) = TwoRarefactionPressure(L,R)
           Lstar%W(iu(ix)) = TwoRarefactionVelocity(L,R)
           Rstar%W(ip) = Lstar%W(ip)
           Rstar%W(iu(ix)) = Lstar%W(iu(ix))
           Lstar%W(irho) = TwoRarefactionDensity(L,Lstar)
           Rstar%W(irho) = TwoRarefactionDensity(R,Rstar)
          else
           Lstar%W(ip) = TwoShockPressure(L,R,p0)
           Lstar%W(iu(ix)) = TwoShockVelocity(L,R,Lstar,p0)
           Rstar%W(ip) = Lstar%W(ip)
           Rstar%W(iu(ix)) = Lstar%W(iu(ix))
           Lstar%W(irho) = TwoShockDensity(L,Lstar)
           Rstar%W(irho) = TwoShockDensity(R,Rstar)
          end if
        end select
        
        return
        
        end subroutine
        
        ! Calculate approximate wave speeds in Riemann problem boundary (for flux approximations)
        subroutine WaveSpeed(L,R,SL,SR,S0)
        double precision, intent(out) :: SL, SR
        double precision, intent(out), optional :: S0
        type(region), intent(in) :: L,R
        double precision p0
        
        p0 = exactPressure(L,R) !max(0d0,PrimitiveVariablePressure(L,R))
        ! Approximate shock or rarefaction speeds
        if (p0 > L%W(ip)) then
         SL = L%W(iu(ix)) - L%c*sqrt(1 + (gamma + 1.)*(p0/L%W(ip) - 1.)/(2.*gamma))
        else
         SL = L%W(iu(ix)) - L%c
        end if
        if (p0 > R%W(ip)) then
         SR = R%W(iu(ix)) + R%c*sqrt(1 + (gamma + 1.)*(p0/R%W(ip) - 1.)/(2.*gamma))
        else
         SR = R%W(iu(ix)) + R%c
        end if
        if (present(S0)) S0 = exactVelocity(p0,L,R) !(R%W(ip) - L%W(ip) + L%W(irho)*L%W(iu(ix))*(SL - L%W(iu(ix))) - R%W(irho)*R%W(iu(ix))*(SR - R%W(iu(ix))))/(L%W(irho)*(SL - L%W(iu(ix))) - R%W(irho)*(SR - R%W(iu(ix))))
        
        return
        
        end subroutine
        
        function sampleWithoutVacuum(L,R,Lstar,Rstar,location)
        double precision, intent(in) :: location
        double precision SL,SHL,STL,SR,SHR,STR
        logical lright
        double precision Wfan(NrVar)
        double precision sampleWithoutVacuum(NrVar)
        type(region), intent(in) :: L,R,Lstar,Rstar
        
        ! Calculate shock and rarefaction speeds
        SL = L%W(iu(ix)) - sqrt((Lstar%W(ip) + L%B)/L%A)/L%W(irho)                                          ! Left shock speed
        SHL = L%W(iu(ix)) - L%c                                                                             ! Head speed of left rarefaction
        STL = Lstar%W(iu(ix)) - L%c*(Lstar%W(ip)/L%W(ip))**((gamma - 1)/(2*gamma))                          ! Tail speed of left rarefaction
        SR = R%W(iu(ix)) + sqrt((Rstar%W(ip) + R%B)/R%A)/R%W(irho)                                          ! right
        SHR = R%W(iu(ix)) + R%c                                                                             ! right
        STR = Rstar%W(iu(ix)) + R%c*(Rstar%W(ip)/R%W(ip))**((gamma - 1)/(2*gamma))                          ! right

        if (location <= Lstar%W(iu(ix))) then !left side of contact
         lright = .false.
         if (Lstar%W(ip) > L%W(ip)) then !shock on left side
          if (location <= SL) then !left side of shock
           sampleWithoutVacuum = L%W
          else
           sampleWithoutVacuum = Lstar%W
          end if
         else !rarefaction on left side
          if (location <= SHL) then !left of the rarefaction fan
           sampleWithoutVacuum = L%W
          else if (location >= STL) then !right of the rarefaction fan
           sampleWithoutVacuum = Lstar%W
          else !in rarefaction fan
           call rarefaction(L,R,Wfan,lright,location)
           sampleWithoutVacuum = Wfan
          end if
         end if
        else !right side of contact
         lright = .true.
         if (Rstar%W(ip) > R%W(ip)) then !shock on right side
          if (location >= SR) then !right side of shock
           sampleWithoutVacuum = R%W
          else
           sampleWithoutVacuum = Rstar%W
          end if
         else !rarefaction on right side
          if (location >= SHR) then !right of the rarefaction fan
           sampleWithoutVacuum = R%W
          else if (location <= STR) then !left of the rarefaction fan
           sampleWithoutVacuum = Rstar%W
          else !in rarefaction fan
           call rarefaction(L,R,Wfan,lright,location)
           sampleWithoutVacuum = Wfan
          end if
         end if
        end if
        
        end function
        
        function sampleWithVacuum(L,R,location)
        double precision, intent(in) :: location
        double precision SvacL,SvacR
        logical lright
        double precision Wfan(NrVar)
        double precision sampleWithVacuum(NrVar)
        type(region) L,R
        
        ! Calculate shock and rarefaction speeds
        SvacL = L%W(iu(ix)) + 2.*L%c/(gamma - 1.)                                                           ! Speed of vacuum front (vacuum on right)
        SvacR = R%W(iu(ix)) - 2.*R%c/(gamma - 1.)                                                           ! Speed of vacuum front (vacuum on left)
        
        ! Initial state is vacuum
        if (vacuum == 1) then
         if (R%W(ip) == 0 .and. R%W(irho) == 0) then !vacuum on the right
          lright = .false.
          R%W(iu(ix)) = 0
          if (location <= L%W(iu(ix)) - L%c) then             ! Left of rarefaction fan
           sampleWithVacuum = L%W
          else if (location >= SvacL) then                ! Right of rarefaction fan (in vacuum)
           sampleWithVacuum = R%W
          else                                       ! In rarefaction fan
           call rarefaction(L,R,Wfan,lright,location)
           sampleWithVacuum = Wfan
          end if
         else                                        !vacuum on the left
          lright = .true.
          L%W(iu(ix)) = 0
          if (location >= R%W(iu(ix)) - R%c) then             ! Right of rarefaction fan
           sampleWithVacuum = R%W
          else if (location <= SvacR) then                ! Left of rarefaction fan (in vacuum)
           sampleWithVacuum = L%W
          else
           call rarefaction(L,R,Wfan,lright,location)       ! In rarefaction fan
           sampleWithVacuum = Wfan
          end if
         end if
        ! Initial regions result in vacuum between
        else if (vacuum == 2) then
         if (location <= L%W(iu(ix)) - L%c) then              ! Left of left rarefaction fan
          sampleWithVacuum = L%W
         else if (location <= SvacL) then                 ! In left rarefaction fan
          lright = .false.
          call rarefaction(L,R,Wfan,lright,location)
          sampleWithVacuum = Wfan
         else if (location >= R%W(iu(ix)) + R%c) then         ! Right of right rarefaction fan
          sampleWithVacuum = R%W
         else if (location >= SvacR) then                 ! In right rarefaction fan
          lright = .true.
          call rarefaction(L,R,Wfan,lright,location)
          sampleWithVacuum = Wfan
         else                                        ! In vacuum (between rarefaction fans)
          sampleWithVacuum(irho) = 0
          sampleWithVacuum(iu(ix)) = 0
          sampleWithVacuum(ip) = 0
          if (nDim > 1) sampleWithVacuum(iu(iy)) = 0
          if (nDim > 2) sampleWithVacuum(iu(iz)) = 0
         endif
        end if
        
        end function
        
        ! Calculate density, velocity, and pressure in rarefaction fan - used in output routine when querying rarefaction areas
        subroutine rarefaction(L,R,Wfan,lright,location)
        double precision, intent(in) :: location
        type(region), intent(in) :: L,R
        logical, intent(in) :: lright
        double precision, intent(out) :: Wfan(NrVar)
        double precision c
        
        ! If fluid has covolume constant, iterate to solve for density
        if (cov > 0) then
         if (.not. lright) then
          beta = ((((gamma - 1.)*(L%W(iu(ix)) + 2.*L%c*(1. - cov*L%W(irho))/(gamma - 1.) - location))**2.)/(gamma*L%W(ip)*((1. - cov*L%W(irho))/L%W(irho))**gamma)) ! Constant related to Riemann invariant
          Wfan(irho) = NRSolver(rhoguess(L,R),densityF,densityDf,L,R)
          Wfan(ip) = L%W(ip)*(((1. - cov*L%W(irho))/L%W(irho))*(Wfan(irho)/(1. - cov*Wfan(irho)))**gamma)
          c = sqrt((Wfan(ip)*gamma)/(Wfan(irho)*(1. - cov*Wfan(irho))))                                         ! Speed of sound at current location
          Wfan(iu(ix)) = location + c
          if (nDim > 1) Wfan(iu(iy)) = L%W(iu(iy))
          if (nDim > 2) Wfan(iu(iz)) = L%W(iu(iz))
         else
          beta = (((location - R%W(iu(ix)) + (2.*R%c*(1. - cov*R%W(irho))/(gamma - 1.)))**2.)/(gamma*R%W(ip)*((1. - cov*R%W(irho))/R%W(irho))**gamma))   ! Constant related to Riemann invariant
          Wfan(irho) = NRSolver(rhoguess(L,R),densityF,densityDf,L,R)
          Wfan(ip) = R%W(ip)*(((1. - cov*R%W(irho))/R%W(irho))*(Wfan(irho)/(1. - cov*Wfan(irho)))**gamma)
          c = sqrt((Wfan(ip)*gamma)/(Wfan(irho)*(1. - cov*Wfan(irho))))                                         ! Speed of sound at current location
          Wfan(iu(ix)) = location - c
          if (nDim > 1) Wfan(iu(iy)) = R%W(iu(iy))
          if (nDim > 2) Wfan(iu(iz)) = R%W(iu(iz))
         end if
        ! Otherwise solution is closed-form
        else
         ! Rarefaction fan differs depending on whether it's a right or left rarefaction
         if (.not. lright) then
          Wfan(irho) = L%W(irho)*(2/(gamma + 1) + (gamma - 1)*(L%W(iu(ix)) - location)/((gamma + 1)*L%c))**(2/(gamma - 1))
          Wfan(iu(ix)) = (2/(gamma + 1))*(L%c + (gamma - 1)*L%W(iu(ix))/2 + location)
          Wfan(ip) = L%W(ip)*(2/(gamma + 1) + (gamma - 1)*(L%W(iu(ix)) - location)/((gamma + 1)*L%c))**(2*gamma/(gamma - 1))
          ! Tangential directions unchanged
          if (nDim > 1) Wfan(iu(iy)) = L%W(iu(iy))
          if (nDim > 2) Wfan(iu(iz)) = L%W(iu(iz))
         else
          Wfan(irho) = R%W(irho)*(2/(gamma + 1) - (gamma - 1)*(R%W(iu(ix)) - location)/((gamma + 1)*R%c))**(2/(gamma - 1))
          Wfan(iu(ix)) = (2/(gamma + 1))*(-R%c + (gamma - 1)*R%W(iu(ix))/2 + location)
          Wfan(ip) = R%W(ip)*(2/(gamma + 1) - (gamma - 1)*(R%W(iu(ix)) - location)/((gamma + 1)*R%c))**(2*gamma/(gamma - 1))
          ! Tangential directions unchanged
          if (nDim > 1) Wfan(iu(iy)) = R%W(iu(iy))
          if (nDim > 2) Wfan(iu(iz)) = R%W(iu(iz))
         end if
        end if
        
        return
        
        end subroutine
        
        end module
