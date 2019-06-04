        ! This program takes a domain, time, boundary conditions, and initial conditions, and calculates the requested state of the fluid using Godunov's method.
        program OneDimensionalGodunovSolver
        use Input
        use line
        use Godunov
        use Output
        implicit none
        
        ! Create and initialize time step counter
        integer :: loops = 1
        
        ! Executable for single variable Godunov solver
        progNum = 2
        
        ! Get values from problem.data
        call InitDomain

        ! Calculate next output time
        nexttime = starttime + (finaltime - starttime)/frames
        
        ! At first run, output initial conditions
        if (frame == 0) then
         print *, "Writing out frame ", frame
         call WriteOut
         frame = frame + 1
        end if
        
        ! Start time-dependent calculation
        do
         ! Print information to console
         print *, "Time step ", loops
         print *, "Time = ", time
         
         ! Set cell size
         dl = dx(ix)            ! dl is global for 1D routines
         cells = xcells
         
         ! Apply boundary conditions
         call BoundaryConditions(ProblemRegion)
         
         ! Calculate delta t for this step
         call TimeStep
         ! Print to console
         print *, "dt = ", dt
         
         ! Calculate flux between cells and update velocity
         call CalcFlux(ProblemRegion)
         call UpdateCells(ProblemRegion)
         
         ! Time of cells is now advanced
         time = time + dt
         
         ! If we're at an output time, write to file
         if (time == nexttime) then
          print *, "Writing out frame ", frame
          call WriteOut
          nexttime = nexttime + (finaltime-starttime)/frames    ! Set next output time
          frame = frame + 1
         end if
         
         ! Stop after final time !!!(still a bit buggy - get an extra time step at the end [rounding?])
         if (time >= finaltime) exit
         
         ! Count number of time steps performed
         loops = loops + 1
        
        end do
        
        ! Clean up
        deallocate(ProblemRegion%U,ProblemRegion%flux,iu,imom,dx)
        
        stop
        
        end program
