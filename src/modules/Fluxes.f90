        ! Intercell fluxes and variable conversion functions
        module equations
        use EulerDeclarations
        use RiemannSolutions
        use regions
        implicit none
        
        contains
        
        function RoeFlux(L,R)
        type(region), intent(in) :: L,R
        integer i
        double precision HL,HR,lambdaLL,lambdaLR,lambdaRL,lambdaRR
        double precision lambda(NrVar),K(NrVar,NrVar),alpha(NrVar),temp(NrVar)
        type(region) RoeAvg,Lstar,Rstar
        double precision RoeFlux(NrVar)
        
        allocate(RoeAvg%W(NrVar))
        
        RoeAvg%W(irho) = sqrt(L%W(irho)*R%W(irho))
        RoeAvg%W(iu) = (sqrt(L%W(irho))*L%W(iu) + sqrt(R%W(irho))*R%W(iu))/(sqrt(L%W(irho)) + sqrt(R%W(irho)))
        HL = (L%W(ip)*(1. - cov*L%W(irho))/(gamma - 1.) + 0.5*L%W(irho)*sum(L%W(iu)**2) + L%W(ip))/L%W(irho)
        HR = (R%W(ip)*(1. - cov*R%W(irho))/(gamma - 1.) + 0.5*R%W(irho)*sum(R%W(iu)**2) + R%W(ip))/R%W(irho)
        RoeAvg%c = sqrt((gamma - 1.)*((sqrt(L%W(irho))*HL + sqrt(R%W(irho))*HR)/(sqrt(L%W(irho)) + sqrt(R%W(irho))) - 0.5*sum(RoeAvg%W(iu)**2)))
        
        lambda = RoeAvg%W(iu(ix))
        lambda(1) = lambda(1) - RoeAvg%c
        lambda(nDim+2) = lambda(nDim+2) + RoeAvg%c
        
        K(2,:) = (/ 1d0, RoeAvg%W(iu), 0.5*sum(RoeAvg%W(iu)**2) /)
        if (nDim > 1) then
         K(1,:) = (/ 1d0, lambda(1), RoeAvg%W(iu(iy:nDim)), (sqrt(L%W(irho))*HL + sqrt(R%W(irho))*HR)/(sqrt(L%W(irho)) + sqrt(R%W(irho))) - RoeAvg%W(iu(ix))*RoeAvg%c /)
         K(nDim+2,:) = (/ 1d0, lambda(nDim+2), RoeAvg%W(iu(iy:nDim)), (sqrt(L%W(irho))*HL + sqrt(R%W(irho))*HR)/(sqrt(L%W(irho)) + sqrt(R%W(irho))) + RoeAvg%W(iu(ix))*RoeAvg%c /)
         if (nDim > 2) then
          K(3,:) = (/ 0d0, 0d0, 1d0, 0d0, RoeAvg%W(iu(iy)) /)
          K(4,:) = (/ 0d0, 0d0, 0d0, 1d0, RoeAvg%W(iu(iz)) /)
         else
          K(3,:) = (/ 0d0, 0d0, 1d0, RoeAvg%W(iu(iy)) /)
         end if
        else
         K(1,:) = (/ 1d0, lambda(1), (sqrt(L%W(irho))*HL + sqrt(R%W(irho))*HR)/(sqrt(L%W(irho)) + sqrt(R%W(irho))) - RoeAvg%W(iu(ix))*RoeAvg%c /)
         K(nDim+2,:) = (/ 1d0, lambda(nDim+2), (sqrt(L%W(irho))*HL + sqrt(R%W(irho))*HR)/(sqrt(L%W(irho)) + sqrt(R%W(irho))) + RoeAvg%W(iu(ix))*RoeAvg%c /)
        end if
        
        alpha(1) = (R%W(ip) - L%W(ip) - RoeAvg%W(irho)*RoeAvg%c*(R%W(iu(ix)) - L%W(iu(ix))))/(2.*RoeAvg%c**2)
        alpha(2) = R%W(irho) - L%W(irho) - (R%W(ip) - L%W(ip))/RoeAvg%c**2
        alpha(nDim+2) = (R%W(ip) - L%W(ip) + RoeAvg%W(irho)*RoeAvg%c*(R%W(iu(ix)) - L%W(iu(ix))))/(2.*RoeAvg%c**2)
        if (nDim > 1) alpha(3) = RoeAvg%W(irho)*(R%W(iu(iy)) - L%W(iu(iy)))
        if (nDim > 2) alpha(4) = RoeAvg%W(irho)*(R%W(iu(iz)) - L%W(iu(iz)))
        
        lambdaLL = L%W(iu(ix)) - L%c
        lambdaRR = R%W(iu(ix)) + R%c
        if (lambdaLL < 0 .or. lambdaRR > 0) then
         RiemannSolver = Exact
         call StarRegion(L,R,Lstar,Rstar)
         RiemannSolver = Roe
         lambdaLR = Lstar%W(iu(ix)) - SoundSpeed(Lstar)
         lambdaRL = Rstar%W(iu(ix)) + SoundSpeed(Rstar)
        end if
        if (lambdaLL < 0 .and. lambdaLR > 0) then
         print *, 'Using entropy fix.'
         lambda(1) = lambdaLL*((lambdaLR - lambda(1))/(lambdaLR - lambdaLL))
        end if
        if (lambdaRL < 0 .and. lambdaRR > 0) then
         print *, 'Using entropy fix.'
         lambda(nDim+2) = lambdaRR*((lambda(nDim+2) - lambdaRL)/(lambdaRR - lambdaRL))
        end if
        
        temp = 0
        do i = 1, NrVar
         temp = temp + alpha(i)*abs(lambda(i))*K(i,:)
        end do
        RoeFlux = 0.5*(EulerFlux(L%W) + EulerFlux(R%W) - temp)
        
        end function
        
        ! Intercell flux for HLLC approximate solver
        function HLLCFlux(L,R,S0,SL,SR,k)
        integer, intent(in), optional :: k
        double precision, intent(in) :: SL,SR,S0
        type(region), intent(in) :: L,R
        double precision U0L(NrVar),U0R(NrVar),UL(NrVar),UR(NrVar)
        double precision HLLCFlux(NrVar)
        
        ! Approximate with flux from left Riemann state
        if ((.not. present(k) .and. 0 <= SL) .or. (present(k) .and. k == 1)) then
         HLLCflux = EulerFlux(L%W)
        ! Approximate with integral average of left star region
        else if ((.not. present(k) .and. SL <= 0 .and. 0 <= S0) .or. (present(k) .and. k == 2)) then
         UL = PrimToCons(L%W)
         U0L(irho) = L%W(irho)*(SL - L%W(iu(ix)))/(SL - S0)
         U0L(imom(ix)) = S0*L%W(irho)*(SL - L%W(iu(ix)))/(SL - S0)
         if (nDim > 1) U0L(imom(iy)) = L%W(iu(iy))*L%W(irho)*(SL - L%W(iu(ix)))/(SL - S0)
         if (nDim > 2) U0L(imom(iz)) = L%W(iu(iz))*L%W(irho)*(SL - L%W(iu(ix)))/(SL - S0)
         U0L(iE) = L%W(irho)*((SL - L%W(iu(ix)))/(SL - S0))*(UL(iE)/L%W(irho) + (S0 - L%W(iu(ix)))*(S0 + L%W(ip)/(L%W(irho)*(SL - L%W(iu(ix))))))
         HLLCflux = EulerFlux(L%W) + SL*(U0L - UL)
        ! Approximate with integral average of right star region
        else if ((.not. present(k) .and. S0 <= 0 .and. 0 <= SR) .or. (present(k) .and. k == 3)) then
         UR = PrimToCons(R%W)
         U0R(irho) = R%W(irho)*(SR - R%W(iu(ix)))/(SR - S0)
         U0R(imom(ix)) = S0*R%W(irho)*(SR - R%W(iu(ix)))/(SR - S0)
         if (nDim > 1) U0R(imom(iy)) = R%W(iu(iy))*R%W(irho)*(SR - R%W(iu(ix)))/(SR - S0)
         if (nDim > 2) U0R(imom(iz)) = R%W(iu(iz))*R%W(irho)*(SR - R%W(iu(ix)))/(SR - S0)
         U0R(iE) = R%W(irho)*((SR - R%W(iu(ix)))/(SR - S0))*(UR(iE)/R%W(irho) + (S0 - R%W(iu(ix)))*(S0 + R%W(ip)/(R%W(irho)*(SR - R%W(iu(ix))))))
         HLLCflux = EulerFlux(R%W) + SR*(U0R - UR)
        ! Approximate with right Riemann state
        else
         HLLCflux = EulerFlux(R%W)
        end if
        
        end function

        ! Intercell flux for HLL approximate solver
        function HLLFlux(L,R,SL,SR)
        double precision, intent(in) :: SL,SR
        type(region), intent(in) :: L,R
        double precision UL(NrVar),UR(NrVar)
        double precision HLLFlux(NrVar)
        
        ! Approximate with left Riemann state
        if (0 <= SL) then
         HLLflux = EulerFlux(L%W)
        ! Approximate with integral average of star region
        else if (SL <= 0 .and. 0 <= SR) then
         UL = PrimToCons(L%W)
         UR = PrimToCons(R%W)
         HLLflux = (SR*EulerFlux(L%W) - SL*EulerFlux(R%W) + SL*SR*(UR - UL))/(SR - SL)
        ! Approximate with right Riemann state
        else
         HLLflux = EulerFlux(R%W)
        end if
        
        end function
        
        ! Intercell flux for Euler equations
        function EulerFlux(W)
        double precision, intent(in) :: W(NrVar)
        double precision EulerFlux(NrVar)
        
        EulerFlux(irho) = W(irho)*W(iu(ix))
        EulerFlux(imom(ix)) = W(irho)*W(iu(ix))**2 + W(ip)
        EulerFlux(iE) = W(iu(ix))*((W(ip)*(1. - cov*W(ip)))/(gamma - 1.) + 0.5*W(irho)*W(iu(ix))**2 + W(ip))
        if (nDim > 1) EulerFlux(imom(iy)) = W(irho)*W(iu(ix))*W(iu(iy))
        if (nDim > 2) EulerFlux(imom(iz)) = W(irho)*W(iu(ix))*W(iu(iz))
        
        end function
        
        ! Equation for intercell flux for invisid Burgers equation
        double precision function BurgersFlux(u)
        double precision, intent(in) :: u        
        
        BurgersFlux = 0.5*u**2
        
        end function

        ! Equation for intercell flux for linear advection equation
        double precision function LinearAdvectionFlux(u)
        double precision, intent(in) :: u

        LinearAdvectionFlux = advectionConst*u

        end function
        
        ! Converts conserved variables (mass, momentum, energy) to primitive variables (mass, velocity, pressure)
        function ConsToPrim(U)
        double precision, intent(in) :: U(NrVar)
        double precision ConsToPrim(NrVar)
        
        ConsToPrim(irho) = U(irho)
        ConsToPrim(iu(ix)) = U(imom(ix))/U(irho)
        ConsToPrim(ip) = (U(iE) - (sum(U(imom)**2))/(2*U(irho)))*(gamma - 1.)/(1. - cov*U(irho))
        ! Tangential directions
        if (nDim > 1) ConsToPrim(iu(iy)) = U(imom(iy))/U(irho)
        if (nDim > 2) ConsToPrim(iu(iz)) = U(imom(iz))/U(irho)

        end function
        
        ! Converts conserved variables (mass, momentum, energy) to primitive variables (mass, velocity, pressure)
        function PrimToCons(W)
        double precision, intent(in) :: W(NrVar)
        double precision PrimToCons(NrVar)
        
        PrimToCons(irho) = W(irho)
        PrimToCons(imom(ix)) = W(iu(ix))*W(irho)
        PrimToCons(iE) = W(ip)*(1. - cov*W(irho))/(gamma - 1.) + 0.5*W(irho)*sum(W(iu)**2)
        ! Tangential directins
        if (nDim > 1) PrimToCons(imom(iy)) = W(iu(iy))*W(irho)
        if (nDim > 2) PrimToCons(imom(iz)) = W(iu(iz))*W(irho)

        end function
        
        double precision function SoundSpeed(reg)
        type(region), intent(in) :: reg
        
        SoundSpeed = sqrt(gamma*reg%W(ip)/(reg%W(irho)*(1. - cov*reg%W(irho))))
        
        end function
        
        double precision function ConsSoundSpeed(U)
        double precision, intent(in) :: U(NrVar)
        
        ConsSoundSpeed = sqrt(gamma*(gamma - 1.)*(U(iE) - (sum(U(imom)**2))/(2*U(irho)))/(U(irho)*(1. - cov*U(irho))**2))
        
        end function
        
        end module
