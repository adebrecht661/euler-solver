fc = gfortran
fflags += -O3 -ffree-line-length-none
debugflags += -g -fcheck=all -Wall -ffree-line-length-none -fbacktrace
OMPflags += -fopenmp

#modules must be in appropriate order
modules = src/modules/Declarations.f90 src/modules/RiemannRegions.f90 src/modules/Domain.f90 src/modules/ApproxRiemannEqns.f90 src/modules/ExactRiemannEqns.f90 src/modules/RiemannSolutions.f90 src/modules/Fluxes.f90 src/modules/Input.f90 src/modules/Output.f90 src/modules/Schemes.f90 src/modules/Sweep.f90

ExactExec = src/ExactRiemann.f90

GodunovExec = src/GodunovSolver.f90

EulerExec = src/EulerSolver.f90

EulerExec_parallel = src/EulerSolver_parallel.f90

Euler : $(modules) $(EulerExec)
	$(fc) $(fflags) $(modules) $(EulerExec) -o EulerSolver
	make clean

parallel_euler : $(modules) $(EulerExec)
	$(fc) $(fflags) $(OMPflags) $(modules) $(EulerExec_parallel) -o EulerSolver
	make clean

Godunov : $(modules) $(GodunovExec)
	$(fc) $(fflags) $(modules) $(GodunovExec) -o GodunovSolver
	make clean

Exact : $(modules) $(ExactExec)
	$(fc) $(fflags) $(modules) $(ExactExec) -o ExactRiemannSolver
	make clean

debug_euler : $(modules) $(EulerExec)
	$(fc) $(debugflags) $(modules) $(EulerExec) -o EulerSolver_debug
	
debug_parallel_euler : $(modules) $(EulerExec_parallel)
	$(fc) $(debugflags) $(OMPflags) $(modules) $(EulerExec) -o EulerSolver_debug

debug_godunov : $(modules) $(GodunovExec)
	$(fc) $(debugflags) $(modules) $(GodunovExec) -o GodunovSolver_debug

debug_exact : $(modules) $(ExactExec)
	$(fc) $(debugflags) $(modules) $(ExactExec) -o ExactRiemannSolver_debug

clean : 
	rm *.mod

delete :
	-rm ExactRiemannSolver
	-rm ExactRiemannSolver_debug
	-rm GodunovSolver
	-rm GodunovSolver_debug
	-rm EulerSolver
	-rm EulerSolver_debug

