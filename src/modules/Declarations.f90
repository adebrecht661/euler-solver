        module EulerDeclarations
        implicit none
        
        ! declare variables for global use (mostly set during input)
        double precision xL,xR,yL,yR,zL,zR,starttime,finaltime,CFL,dt,time,nexttime,tol,cov,gamma,beta,Q,advectionConst,dl
        integer frames,xcells,ycells,zcells,frame,irho,ip,iE,vacuum,NrVar,nDim,maxiter,RiemannSolver,progNum,Eq,cells,ix,iy,iz,FluxScheme
        double precision, allocatable :: dx(:)
        integer, allocatable :: imom(:),iu(:)
        ! Initialize some options to defaults
        integer :: BC = 2
        integer :: num_threads = 1
        logical :: lthreading = .false.
        
        ! Flux schemes
        integer, parameter :: GodunovUpwind = 1, WAF = 2
        
        ! Program control - should be set in main source file
        integer, parameter :: Euler = 3, GodunovSingle = 2, ExactProg = 1

        end module
