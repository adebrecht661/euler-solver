        ! This program takes two fluid states, assumed to meet at the origin, and calculates the density, pressure, and velocity over time of the resulting states.
        program OneDimensionalRiemannSolver
        use Input
        use output
        use regions
        use ExactSolver
        implicit none
        
        ! Executable for exact Riemann solver
        progNum = 1
        
        ! Get values from problem.data, set up Riemann problem
        call InitDomain
        
        ! Check to see if there's vacuum in the simulation
        call check_vacuum(L,R)

        ! If vacuum is present, just write out.
        if (vacuum /= 3) then
         call WriteOut()  ! no star region with vacuum (just vacuum and rarefaction fans)
         deallocate (L%W,R%W)
        else
         ! Force exact Riemann solver
         RiemannSolver = 1
         call StarRegion(L,R,Lstar,Rstar)
        
         ! Print pressure in star regions (mostly for comparison to tests)
         print *, 'Final pstar = ', Lstar%W(ip)
        
         ! Write pressure, velocity, and density to output files (one per frame)
         print *, 'Writing out.'
         call WriteOut()
        
         ! Clean up
         deallocate(L%W,R%W,Lstar%W,Rstar%W)
        end if
        
        ! Clean up from input (other programs)
        deallocate(iu,imom,dx)
        
        stop
        
        end program
        
        
        
        
        
        
        
        
        
        
        
        
