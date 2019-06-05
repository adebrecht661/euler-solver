        ! defines type for problem domain
        module line
        use EulerDeclarations
        implicit none
        
        ! Selection parameters for boundary conditions
        integer, parameter :: periodic = 1, transparent = 2, reflective = 3
        
        ! Derived type for collecting arrays of flux and conserved variables in the problem domain
        type OneDdomain
         sequence
         double precision, allocatable :: U(:,:), flux(:,:)
        end type
        
        ! 2D and 3D domains don't need flux, since problem is reduced to one dimension for each calculation
        type TwoDdomain
         sequence
         double precision, allocatable :: U(:,:,:)
        end type
        
        type ThreeDdomain
         sequence
         double precision, allocatable :: U(:,:,:,:)
        end type
        
        ! Declare problem domain for global use
        type(OneDdomain) ProblemRegion
        type(TwoDdomain) ProblemRegion2D
        type(ThreeDdomain) ProblemRegion3D

        contains
        
        ! Apply boundary conditions to problem
        subroutine BoundaryConditions(ProblemRegion)
        type(OneDdomain), intent(inout) :: ProblemRegion
        
        select case (BC)
         case(periodic)         ! Wraps around
          ProblemRegion%U(0,1:NrVar) = ProblemRegion%U(cells,1:NrVar)
          ProblemRegion%U(cells+1,1:NrVar) = ProblemRegion%U(1,1:NrVar)
         case(transparent)      ! Wave flows out
          ProblemRegion%U(-1,1:NrVar) = ProblemRegion%U(1,1:NrVar)
          ProblemRegion%U(0,1:NrVar) = ProblemRegion%U(1,1:NrVar)
          ProblemRegion%U(cells+1,1:NrVar) = ProblemRegion%U(cells,1:NrVar)
          ProblemRegion%U(cells+2,1:NrVar) = ProblemRegion%U(cells,1:NrVar)
         case(reflective)       ! Wave bounces off hard boundary
          ProblemRegion%U(-1,1:NrVar) = ProblemRegion%U(2,1:NrVar)
          ProblemRegion%U(-1,imom(ix)) = -ProblemRegion%U(2,imom(ix))  ! Momentum normal to boundary is reversed
          ProblemRegion%U(0,1:NrVar) = ProblemRegion%U(1,1:NrVar)
          ProblemRegion%U(0,imom(ix)) = -ProblemRegion%U(1,imom(ix))
          ProblemRegion%U(cells+1,1:NrVar) = ProblemRegion%U(cells,1:NrVar)
          ProblemRegion%U(cells+1,imom(ix)) = -ProblemRegion%U(cells,imom(ix))
          ProblemRegion%U(cells+2,1:NrVar) = ProblemRegion%U(cells-1,1:NrVar)
          ProblemRegion%U(cells+2,imom(ix)) = -ProblemRegion%U(cells-1,imom(ix))
        end select
        
        end subroutine

        end module
