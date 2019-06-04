        ! This program takes a domain, time, boundary conditions, and test number, and calculates the requested state of the fluid using Godunov's first order upwind method.
        program MultiDimensionalEulerSolver
        use Input
        use line
        use sweep
        use Output
        use EulerDeclarations
        implicit none
        
        ! Create loop counters and initialize time step counter
        integer :: loops = 1
        integer i,j,k,time1,time2,count_rate
        
        ! Executable for Euler solver using Godunov method
        progNum = 3
        
        ! Get values from problem.data and initialize problem domain
        call InitDomain

        ! Calculate next output time
        nexttime = starttime + (finaltime - starttime)/frames
        
        ! At first run, output initial conditions
        if (frame == 0) then
         call system_clock(time1)
         print *, "Writing out frame ", frame
         call WriteOut
         call system_clock(time2,count_rate)
         print *, 'Time to output frame ', frame, ' = ', (time2-time1)/dble(count_rate),'s.'
         frame = frame + 1
        end if
        
        ! Start time-dependent calculation
        do
        
         ! Print information to console
         print *, "Time step ", loops
         print *, "Time = ", time
         
         ! Calculate delta t for this step (uses all dimensions present)
         call TimeStep
         ! Print to console
         print *, "dt = ", dt
         
         ! Do a sweep through each dimension, depending on how many there are
         select case (nDim)
           
          ! 1D   
          case (1)
           cells = xcells               ! cells is global for 1D solver
           dl = dx(ix)                  ! dl similarly for cell size
           ! Apply boundary conditions
           call BoundaryConditions(ProblemRegion)
           ! Calculate flux between cells and update velocity
           call CalcFlux(ProblemRegion)
           call UpdateCells(ProblemRegion)
          
          ! 2D
          case (2)
           
           ! x sweep
           cells = xcells               ! cells is global for Godunov 1D solver
           dl = dx(ix)                  ! dl similarly for cell size
           
           ! For each strip in x direction, copy to 1D region, then solve
           do j = 1, ycells
            ProblemRegion%U(:,:) = ProblemRegion2D%U(:,j,:)
            
            ! Apply boundary conditions
            call BoundaryConditions(ProblemRegion)
            
            ! Calculate flux between cells and update variables
            call CalcFlux(ProblemRegion)
            call UpdateCells(ProblemRegion)
            
            ! Update 2D grid from 1D solution
            ProblemRegion2D%U(:,j,:) = ProblemRegion%U(:,:)
           end do
           
           ! y sweep
           cells = ycells               ! cells is global for Godunov 1D solver
           dl = dx(iy)                  ! dl similarly for cell size
           
           ! For each strip in y direction, copy to 1D region, then solve
           do i = 1, xcells
            ProblemRegion%U(:,:) = ProblemRegion2D%U(i,:,:)
            ! x is considered the direction normal to the cells in 1D routines - so switch v to u for y sweep
            ProblemRegion%U(:,imom(ix)) = ProblemRegion2D%U(i,:,imom(iy))
            ProblemRegion%U(:,imom(iy)) = ProblemRegion2D%U(i,:,imom(ix))
            
            ! Apply boundary conditions
            call BoundaryConditions(ProblemRegion)
            
            ! Calculate flux between cells and update velocity
            call CalcFlux(ProblemRegion)
            call UpdateCells(ProblemRegion)
            
            ! Update 2D grid
            ProblemRegion2D%U(i,:,:) = ProblemRegion%U(:,:)
            ! Switch speeds back
            ProblemRegion2D%U(i,:,imom(iy)) = ProblemRegion%U(:,imom(ix))
            ProblemRegion2D%U(i,:,imom(ix)) = ProblemRegion%U(:,imom(iy))
           end do
          
          ! 3D
          case (3)
           ! x sweep
           cells = xcells               ! cells is global for Godunov 1D solver
           dl = dx(ix)                  ! dl similarly for cell size
           
           ! For each strip in x direction, copy to 1D region, then solve
           do j = 1, ycells
            do k = 1, zcells
             ProblemRegion%U(:,:) = ProblemRegion3D%U(:,j,k,:)
             
             ! Apply boundary conditions
             call BoundaryConditions(ProblemRegion)
             
             ! Calculate flux between cells and update velocity
             call CalcFlux(ProblemRegion)
             call UpdateCells(ProblemRegion)
             
             ! Update 3D grid
             ProblemRegion3D%U(:,j,k,:) = ProblemRegion%U(:,:)
            end do
           end do
           
           ! y sweep
           cells = ycells               ! cells is global for Godunov 1D solver
           dl = dx(iy)                  ! dl similarly for cell size
           
           ! For each strip in y direction, copy to 1D region, then solve
           do i = 1, xcells
            do k = 1, zcells
             ProblemRegion%U(:,:) = ProblemRegion3D%U(i,:,k,:)
             ! x is considered the direction normal to the cells in 1D routines - so switch v to u for y sweep (w doesn't change - still passively advected)
             ProblemRegion%U(:,imom(ix)) = ProblemRegion3D%U(i,:,k,imom(iy))
             ProblemRegion%U(:,imom(iy)) = ProblemRegion3D%U(i,:,k,imom(ix))
             
             ! Apply boundary conditions
             call BoundaryConditions(ProblemRegion)

             ! Calculate flux between cells and update velocity
             call CalcFlux(ProblemRegion)
             call UpdateCells(ProblemRegion)

             ! Update 3D grid
             ProblemRegion3D%U(i,:,k,:) = ProblemRegion%U(:,:)
             ! switch speeds back
             ProblemRegion3D%U(i,:,k,imom(ix)) = ProblemRegion%U(:,imom(iy))
             ProblemRegion3D%U(i,:,k,imom(iy)) = ProblemRegion%U(:,imom(ix))
            end do
           end do
           
           ! z sweep
           cells = zcells               ! cells is global for Godunov 1D solver
           dl = dx(iz)                  ! dl similarly for cell size
           ! For each strip in z direction, copy to 1D region, then solve
           do i = 1, xcells
            do j = 1, ycells
             
             ProblemRegion%U(:,:) = ProblemRegion3D%U(i,j,:,:)
             ! x is considered the direction normal to the cells for 1D routines - so switch w to u for y sweep (v doesn't change - still passively advected)
             ProblemRegion%U(:,imom(ix)) = ProblemRegion3D%U(i,j,:,imom(iz))
             ProblemRegion%U(:,imom(iz)) = ProblemRegion3D%U(i,j,:,imom(ix))
             
             ! Apply boundary conditions
             call BoundaryConditions(ProblemRegion)
             
             ! Calculate flux between cells and update velocity
             call CalcFlux(ProblemRegion)
             call UpdateCells(ProblemRegion)
             
             ! Update 3D grid
             ProblemRegion3D%U(i,j,:,:) = ProblemRegion%U(:,:)
             ! switch speeds back
             ProblemRegion3D%U(i,j,:,imom(ix)) = ProblemRegion%U(:,imom(iz))
             ProblemRegion3D%U(i,j,:,imom(iz)) = ProblemRegion%U(:,imom(ix))
             
            end do
           end do
             
         end select
         
         ! Grid is now at next time state
         time = time + dt
           
         ! If we're at an output time, write to file
         if (time == nexttime) then
          print *, "Writing out frame ", frame
          call system_clock(time1)
          call WriteOut
          nexttime = nexttime + (finaltime-starttime)/frames    ! Set next output time
          call system_clock(time2,count_rate)
          print *, 'Time to output frame ', frame, ' = ', (time2-time1)/dble(count_rate),'s.'
          frame = frame + 1
         end if
 
         ! Stop after final time !!! still a bit buggy - get an extra time step at the end (rounding?) !!!
         if (time >= finaltime) exit
        
         ! Count number of time steps performed
         loops = loops + 1
         !if (loops > 1) stop
        
        end do
        
        deallocate(ProblemRegion%U,ProblemRegion%flux,imom,iu,dx)
        if (allocated(ProblemRegion2D%U)) deallocate(ProblemRegion2D%U)
        if (allocated(ProblemRegion3D%U)) deallocate(ProblemRegion3D%U)
        
        stop
        
        end program
