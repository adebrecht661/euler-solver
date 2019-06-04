        ! Functions for various approximate Riemann solvers
        !!! Should look for covolume formulae
        module ApproximateSolver
        use regions
        use EulerDeclarations
        implicit none
        
        contains
        
        ! Functions for primitive variable solver
        double precision function PrimitiveVariablePressure(L,R)
        type(region), intent(in) :: L,R
        
        PrimitiveVariablePressure = 0.5*(L%W(ip) + R%W(ip)) + 0.125*(L%W(iu(ix)) - R%W(iu(ix)))*(L%W(irho) + R%W(irho))*(L%c + R%c)
        
        end function
        
        double precision function PrimitiveVariableVelocity(L,R)
        type(region), intent(in) :: L,R
        
        PrimitiveVariableVelocity = 0.5*(L%W(iu(ix)) + R%W(iu(ix))) + 2.*(L%W(ip) - R%W(ip))/((L%W(irho) + R%W(irho))*(L%c + R%c))
        
        end function
        
        double precision function PrimitiveVariableDensity(reg,starReg,otherReg)
        type(region), intent(in) :: reg,starReg,otherReg
        
        PrimitiveVariableDensity = reg%W(irho) + (reg%W(iu(ix)) - starReg%W(iu(ix)))*(reg%W(irho) + otherReg%W(irho))/(reg%c + otherReg%c)
        
        end function
        
        ! Functions for two-shock solver
        double precision function TwoShockPressure(L,R,p0)
        type(region), intent(in) :: L,R
        double precision, intent(in) :: p0
        
        TwoShockPressure = (sqrt(L%A/(p0 + L%B))*L%W(ip) + sqrt(R%A/(p0 + R%B))*R%W(ip) - R%W(iu(ix)) + L%W(iu(ix)))/(sqrt(L%A/(p0 + L%B)) + sqrt(R%A/(p0 + R%B)))
        
        end function
        
        double precision function TwoShockVelocity(L,R,starReg,p0)
        double precision, intent(in) :: p0
        type(region), intent(in) :: L,R,starReg
        
        TwoShockVelocity = 0.5*(L%W(iu(ix)) + R%W(iu(ix))) + 0.5*((starReg%W(ip) - R%W(ip))*sqrt(R%A/(p0 + R%B)) - (starReg%W(ip) - L%W(ip))*sqrt(L%A/(p0 + L%B)))
        
        end function
        
        double precision function TwoShockDensity(reg,starReg)
        type(region), intent(in) :: reg,starReg
        
        TwoShockDensity = reg%W(irho)*((gamma - 1.)/(gamma + 1.) + (starReg%W(ip)/reg%W(ip)))/(((gamma - 1. + 2.*cov*reg%W(irho))/(gamma + 1.))*(starReg%W(ip)/reg%W(ip)) + ((gamma + 1 - 2.*cov*reg%W(irho))/(gamma + 1)))
        
        end function
        
        ! Functions for two-rarefaction solver
        double precision function TwoRarefactionPressure(L,R)
        type(region), intent(in) :: L,R
        
        TwoRarefactionPressure = ((L%c + R%c - (gamma - 1.)*(R%W(iu(ix)) - L%W(iu(ix)))/2.)/(L%c/(L%W(ip)**((gamma - 1.)/(2.*gamma))) + R%c/(R%W(ip)**((gamma - 1.)/(2.*gamma)))))**(2.*gamma/(gamma - 1.))
        
        end function
        
        double precision function TwoRarefactionVelocity(L,R)
        type(region), intent(in) :: L,R
        double precision PLR
        
        PLR = (L%W(ip)/R%W(ip))**((gamma - 1.)/(2.*gamma))
        
        TwoRarefactionVelocity = (PLR*L%W(iu(ix))/L%c + R%W(iu(ix))/R%c + 2*(PLR - 1.)/(gamma - 1.))/(PLR/L%c + 1./R%c)
        
        end function
        
        double precision function TwoRarefactionDensity(reg,starReg)
        type(region), intent(in) :: reg,starReg
        
        TwoRarefactionDensity = reg%W(irho)*((1. - cov*starReg%W(irho))/(1. - cov*reg%W(irho)))*(starReg%W(ip)/reg%W(ip))**(1./gamma)
        
        end function
        
        end module
