        ! Functions and rountine for solving the exact Riemann equations with the Newton-Raphson iterative method.
        module ExactSolver
        use regions
        use EulerDeclarations
        implicit none
        
        contains
        
        ! Newton-Raphson iterative solver. Takes guess of root of equation, equation, and derivative and calculates root to tol.
        double precision function NRSolver(guess,f,df,L,R)
        double precision guess
        double precision guess_old, delta
        double precision f,df
        type(region), intent(in) :: L,R
        integer i

        i = 0        
        
        ! do while loop (exits once change is less than tolerance)
        do
         i = i + 1
         guess_old = guess
         guess = guess_old - f(guess_old,L,R)/df(guess_old,L,R)           ! Newton-Raphson formula
         delta = 2.*abs((guess - guess_old)/(guess + guess_old))  ! Relative change in pressure from previous iteration
         if (guess /= guess .or. guess < 0) guess = tol                               ! If we get a negative value (unphysical for our applications), set guess to tol
         if (delta < tol) then                                    ! Exit upon sufficient accuracy
          NRSolver = guess
          exit
         end if
         if (i >= maxiter) then
          print *, 'Maximum number of iterations reached. Halting.'
          stop
         end if
        end do
        
        end function
        
        ! Following functions calculate initial guess for pressure and density in Newton-Raphson solver
        
        !!! should improve these initial guess functions, but I'm not sure how
        ! Makes a guess for pressure in the star region. Currently just average of initial pressures.
        double precision function pguess(L,R)
        type(region), intent(in) :: L,R

        pguess = max(tol,0.5*(L%W(ip) + R%W(ip)))

        end function
        
        ! Makes a guess for density in the rarefaction fan (only used if gas has covolume value). Currently just average of initial densities.
        double precision function rhoguess(L,R)
        type(region), intent(in) :: L,R

        rhoguess = max(tol,0.5*(L%W(irho) + R%W(irho)))

        end function

        ! Following functions are for calculations required for Newton-Raphson solver for pressure

        ! Calculates f as a function of pressure in the star region and the given state
        double precision function f(pstar,reg)
        double precision, intent(in) :: pstar
        type(region), intent(in) :: reg
        
        ! If pressure in star region is greater than pressure in given state, shock forms. Otherwise, rarefaction fan occurs.
        if (pstar > reg%W(ip)) then
         !print *, 'Shock.'
         f = (pstar - reg%W(ip))*sqrt(reg%A/(pstar + reg%B))
        else
         !print *, 'Rarefaction.'
         f =  (2./(gamma - 1.))*reg%c*(1 - cov*reg%W(irho))*((pstar/reg%W(ip))**((gamma - 1.)/(2.*gamma)) - 1.)
        end if
        
        end function
        
        ! Calculates total f for initial states for passing to Newton-Raphson solver
        double precision function covolumeF(pstar,L,R)
        double precision, intent(in) :: pstar
        type(region), intent(in) :: L,R

        covolumeF = f(pstar,L) + f(pstar,R) + R%W(iu(ix)) - L%W(iu(ix))
        
        end function
        
        ! Calculates derivative as function of pressure in star region and given state
        double precision function df(pstar,reg)
        double precision, intent(in) :: pstar
        type(region), intent(in) :: reg
        
        ! Same as above.
        if (pstar > reg%W(ip)) then
         !print *, 'Shock.'
         df = sqrt(reg%A/(pstar + reg%B))*(1. - ((pstar - reg%W(ip))/(2.*(reg%B + pstar))))
        else
         !print *, 'Rarefaction.'
         df = ((1. - cov*reg%W(irho))/(reg%W(irho)*reg%c))*((pstar/reg%W(ip))**(-1.*(gamma + 1.)/(2.*gamma)))
        end if
        
        end function
        
        ! Calculates the derivative of the total function above
        double precision function covolumeDf(pstar,L,R)
        double precision, intent(in) :: pstar
        type(region), intent(in) :: L,R
        
        covolumeDf = df(pstar,L) + df(pstar,R)
        
        end function
        
        ! Following functions are for Newton-Raphson solver for density in rarefaction fan when molecules have finite size
        
        ! Calculates function in rarefaction fan as function of density, position, and initial state for passing to Newton-Raphson solver for covolume densities
        double precision function densityF(rho,L,R)
        double precision, intent(in) :: rho
        type(region), intent(in), optional :: L,R

        densityF = (rho**(gamma - 1.))*(gamma + 1. - 2.*cov*rho)**2. - beta*(1. - cov*rho)**(gamma + 1.)        ! beta is set in output module
        
        end function
        
        ! Calculates derivative of function above
        double precision function densityDf(rho,L,R)
        double precision, intent(in) :: rho
        type(region), intent(in), optional :: L,R
        
        densityDf = (gamma + 1.)*(cov*beta*(1. - cov*rho)**gamma + (gamma + 1. - 2.*cov*rho)*(gamma - 1. - 2*cov*rho)*rho**(gamma - 2.))        ! beta is set in output module
                
        end function
        
        ! Remaining two calculate exact values for velocity and density in star regions
        
        ! Calculates velocity in star region
        double precision function exactVelocity(pstar,L,R)
        double precision, intent(in) :: pstar
        type(region), intent(in) :: L,R
        
        exactVelocity = 0.5*(L%W(iu(ix)) + R%W(iu(ix))) + 0.5*(f(pstar,R) - f(pstar,L)) ! Velocity in star region is the same on either side of CD
        
        end function
        
        ! Calculates densities in star regions
        double precision function exactDensity(starReg,reg)
        type(region), intent(in) :: starReg,reg
        
        ! Density is different for shocks and rarefactions - may differ on either side of CD
        if (starReg%W(ip) > reg%W(ip)) then
         exactDensity = reg%W(irho)*((gamma - 1.)/(gamma + 1.) + (starReg%W(ip)/reg%W(ip)))/(((gamma - 1. + 2.*cov*reg%W(irho))/(gamma + 1.))*(starReg%W(ip)/reg%W(ip)) + ((gamma + 1 - 2.*cov*reg%W(irho))/(gamma + 1)))
        else
         exactDensity = 1./(cov + ((1. - cov*reg%W(irho))*(reg%W(ip)/starReg%W(ip))**(1./gamma))/reg%W(irho))
        end if
        
        end function
        
        double precision function exactPressure(L,R)
        type(region), intent(in) :: L,R
        
        exactPressure = NRSolver(pguess(L,R),covolumeF,covolumeDf,L,R)
        
        end function
        
        end module
