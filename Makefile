#  Make File for calculate_flux project
#

OBJ  = time.o declaration.o QC.o calculate_flux.o 
		
F  = pgf90

#-------------------------------------------------------------------

all: $(OBJ) 
	$(F)  -o calculate_flux -Kieee $(OBJ) 
	$(F)  -o read_data src/read_data.f90
	$(F)  -o split_data src/split_data.f90

	rm *.o *.mod

#-------------------------------------------------------------------

calculate_flux:$(OBJ)
	$(F)  -o calculate_flux $(OBJ)
	rm *.o *.mod
#-------------------------------------------------------------------

time.o:
	$(F)  -c -o time.o src/time.f90

declaration.o:
	$(F) -c -o declaration.o src/declaration.f90

QC.o:
	$(F)  -c -o QC.o src/QC.f90

calculate_flux.o: 
	$(F)  -c -o calculate_flux.o -Kieee src/calculate_flux.f90

#-------------------------------------------------------------------
clean:
	rm calculate_flux
	rm read_data
	rm split_data
