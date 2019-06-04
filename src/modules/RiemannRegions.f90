        ! defines type for region variables and rountine to calculate constants
        module regions
        use EulerDeclarations
        implicit none
        
        ! Derived type for consolidating W vector, A, B constants, and speed of sound in region 
        type region
         sequence
         double precision, allocatable :: W(:)
         double precision A,B,c
        end type region
        
        ! Interface to overload default structure constructor
        interface region
         module procedure ConstructRegion
         module procedure ConstructRegion2
        end interface

        ! Declare regions of global use
        type(region) L,R,Lstar,Rstar
        
        contains
        
        ! Constructor function for region that calculates speed of sound, A and B constants as well
        function ConstructRegion(rho,u,p,v,w)
        double precision, intent(in) :: rho,u,p
        double precision, intent(in), optional :: v,w
        type(region) ConstructRegion
        
        ! Allocate region's state array
        allocate(ConstructRegion%W(NrVar))
        
        ! Fill state
        ConstructRegion%W(irho) = rho
        ConstructRegion%W(iu(ix)) = u
        ConstructRegion%W(ip) = p
        if (nDim > 1) ConstructRegion%W(iu(iy)) = v
        if (nDim > 2) ConstructRegion%W(iu(iz)) = w
        ! Calculate useful constants
        if (.not. (rho == 0 .and. p == 0)) then              ! state is not vacuum
         ConstructRegion%A = 2.*(1. - cov*ConstructRegion%W(irho))/((gamma + 1.)*ConstructRegion%W(irho))
         ConstructRegion%B = (gamma - 1.)*ConstructRegion%W(ip)/(gamma + 1.)
         ConstructRegion%c = sqrt(gamma*ConstructRegion%W(ip)/(ConstructRegion%W(irho)*(1. - cov*ConstructRegion%W(irho))))
        else
         ConstructRegion%c = 0          ! Set speed of sound to zero if state is vacuum (to avoid NaN)          !!! could also set A,B to zero, but since they don't represent physical quantities, not much point
        end if
        
        end function
        
        ! Constructor function for region that calculates speed of sound, A and B constants as well
        function ConstructRegion2(W_in)
        double precision, intent(in), dimension(NrVar) :: W_in
        type(region) ConstructRegion2
        
        ! Allocate region's state array
        allocate(ConstructRegion2%W(NrVar))
        
        ! Fill state
        ConstructRegion2%W(irho) = W_in(irho)
        ConstructRegion2%W(iu(ix)) = W_in(iu(ix))
        ConstructRegion2%W(ip) = W_in(ip)
        if (nDim > 1) ConstructRegion2%W(iu(iy)) = W_in(iu(iy))
        if (nDim > 2) ConstructRegion2%W(iu(iz)) = W_in(iu(iz))
        ! Calculate useful constants
        if (.not. (W_in(irho) == 0 .and. W_in(ip) == 0)) then              ! state is not vacuum
         ConstructRegion2%A = 2.*(1. - cov*ConstructRegion2%W(irho))/((gamma + 1.)*ConstructRegion2%W(irho))
         ConstructRegion2%B = (gamma - 1.)*ConstructRegion2%W(ip)/(gamma + 1.)
         ConstructRegion2%c = sqrt(gamma*ConstructRegion2%W(ip)/(ConstructRegion2%W(irho)*(1. - cov*ConstructRegion2%W(irho))))
        else
         ConstructRegion2%c = 0          ! Set speed of sound to zero if state is vacuum (to avoid NaN)          !!! could also set A,B to zero, but since they don't represent physical quantities, not much point
        end if
        
        end function
        
        ! Validates input for creation of vacuum
        subroutine check_vacuum(L,R)
        double precision critdiff
        type(region), intent(in) :: L,R
        
        critdiff = 2.*L%c/(gamma - 1.) + 2.*R%c/(gamma - 1.)      ! Difference between right and left velocities must be strictly less than this for pressure to remain positive
        
        if ((L%W(ip) == 0 .and. L%W(irho) == 0) .or. (R%W(ip) == 0 .and. R%w(irho) == 0)) then  ! one of the initial states is vacuum
         vacuum = 1
         print *, 'One of initial states is vacuum.'
        else if (critdiff <= R%W(iu(ix)) - L%W(iu(ix))) then             ! vacuum is created by initial states
         vacuum = 2
         print *, 'Vacuum is created by initial states.'
        else            ! no vacuum
         vacuum = 3
        end if
        
        return
        
        end subroutine

        end module
