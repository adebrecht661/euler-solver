        ! Writes out to files (one per timestate)
        module Output
        use line
        use equations
        use EulerDeclarations
        use RiemannSolutions
        implicit none
        
        contains
        
        ! Write requested information to files
        subroutine WriteOut()
        integer stat,i,j,k
        character(len=20) :: fileformat,filename
        double precision U(NrVar),x(nDim)
        
        fileformat = "(A,I5.5,A)"       ! Filename format - characters, integer with 5 leading zeros, characters
        write(filename,fileformat) "out/result",frame,".dat"   ! filename to output - goes to out/resultxxxxx.dat
        open(unit=40,file=filename,iostat=stat,status='replace')
        
        ! Depending on which executable was compiled, behave differently
        select case (progNum)

         ! Euler solver with Godunov method
         case (Euler) 
          if (nDim == 1) then
           write (40,*) 'x rho u p'                               ! write column headers
           do i = 1, xcells
            x(ix) = xL + (i-0.5)*dx(ix)                    ! Data is located at center of each cell
            U = ProblemRegion%U(i,:)
            !print *, x, ConsToPrim(U)
            write (40,*) x, ConsToPrim(U)
           end do
          else if (nDim == 2) then
           write (40,*) 'x y rho u v p'
           do i = 1, xcells
            x(ix) = xL + (i-0.5)*dx(ix)                    ! Data is located at center of each cell
            do j = 1, ycells
             x(iy) = yL + (j-0.5)*dx(iy)                    ! Data is located at center of each cell
             U = ProblemRegion2D%U(i,j,:)
             write (40,*) x, ConsToPrim(U)
            end do
           end do
          else if (nDim == 3) then
           write (40,*) 'x y z rho u v w p'
           do i = 1, xcells
            x(1) = xL + (i-0.5)*dx(1)                    ! Data is located at center of each cell
            do j = 1, ycells
             x(2) = yL + (j-0.5)*dx(2)                    ! Data is located at center of each cell
             do k = 1, zcells
              x(3) = zL + (k-0.5)*dx(3)                    ! Data is located at center of each cell
              U = ProblemRegion3D%U(i,j,k,:)
              write (40,*) x, ConsToPrim(U)
             end do
            end do
           end do
          end if
          
         ! Godunov method for one variable
         case (GodunovSingle)
          write (40,*) 'x u'                               ! write column headers
          do i = 1, xcells
           x(1) = xL + (i-0.5)*dx(1)                    ! Data is located at center of each cell
           write (40,*) x, ProblemRegion%u(i,1)
          end do
         
         ! Exact Riemann solver
         case (ExactProg)
          dt = (finaltime-starttime)/frames
          ! Loop through requested times
          do while (time <= finaltime)
           write(filename,fileformat) "out/result",frame,".dat"   ! filename to output - goes to out/resultxxxxx.dat
           open(unit=40,file=filename,iostat=stat,status='replace')
           write (40,*) 'x rho u p'                               ! write column headers
           x(ix) = xL
           ! Calculate values for each position - similar to time above
           do while (x(ix) <= xR)
            if (time == 0) then      ! initial state (so we don't divide by zero)
             if (x(ix) < 0) write (40,*) x, L%W
             if (x(ix) > 0) write (40,*) x, R%W
             if (x(ix) == 0) write (40,*) x, 0, 0, 0
            else
             if (vacuum /= 3) then
              write (40, *) x, sample(L,R,x(ix)/time)
             else
              write (40,*) x, sample(L,R,Lstar,Rstar,x(ix)/time)
             end if
            end if
            x(ix) = x(ix) + dx(ix)            ! increment x to next position
           end do
           time = time + dt             ! increment t to next time
           frame = frame + 1      ! update frame counter
          end do
          
        end select
         
        return
        
        end subroutine
        
        end module
