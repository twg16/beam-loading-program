CXX = g++
MPIXX = mpicxx 
CXXFLAGS = -std=c++11	-O2
HDRS = functions.h
OBJS = main.o  builders.o 
LDLIBS = -llapack -lblas
MPILDLIBS = -lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -llapack -lblas

# problem parameters
beamLength = 10
numberOfNodes = 25
crossSecArea = 0.012
2momentArea = 0.0000144
yMod = 210000000000
time = 1
timesteps = 10000
rho = 7850

default: compile task1 clean

all: compile task1 task2 task3 task4 task5 clean

%.o : %.cpp $(HDRS)
	$(MPIXX) $(CXXFLAGS)	-o   $@	-c $<

compile: $(OBJS)
	$(MPIXX) -o $@ $^ $(MPILDLIBS)

task1:	$(OBJS)
	time -p mpiexec -np 1 ./compile 1 $(beamLength) $(numberOfNodes) $(crossSecArea) $(2momentArea) $(yMod) $(time) $(timesteps) $(rho)

task2: $(OBJS)
	time -p mpiexec -np 1 ./compile 2 $(beamLength) $(numberOfNodes) $(crossSecArea) $(2momentArea) $(yMod) $(time) $(timesteps) $(rho)

task3: $(OBJS)
	time -p mpiexec -np 1 ./compile 3 $(beamLength) $(numberOfNodes) $(crossSecArea) $(2momentArea) $(yMod) $(time) $(timesteps) $(rho)

task4: $(OBJS)
	time -p mpiexec -np 2 ./compile 4 $(beamLength) $(numberOfNodes) $(crossSecArea) $(2momentArea) $(yMod) $(time) $(timesteps) $(rho)

task5: $(OBJS)
	time -p mpiexec -np 2 ./compile 5 $(beamLength) $(numberOfNodes) $(crossSecArea) $(2momentArea) $(yMod) $(time) $(timesteps) $(rho)

clean:
	rm *.o && rm compile



