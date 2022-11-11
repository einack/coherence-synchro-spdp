SRC_DIR=src
OBJ_DIR=objects
BIN_DIR=bin
RESULTS_DIR=results 
FC=gfortran

FFLAGS= -O3 -std=f2003 -finline-functions #-fcheck=all -Wtabs -Wall -Wunused-variable  #-DDEBUG 

##############################################################################

stdp.x: setup $(BIN_DIR)/stdp.x 

$(BIN_DIR)/stdp.x: $(OBJ_DIR)/stdp.o  
	$(FC)  -o $@ $^

##############################################################################

stdp: setup $(OBJ_DIR)/stdp.o 

$(OBJ_DIR)/stdp.o: $(SRC_DIR)/Source_Code_Github.f95 
	$(FC) $(FFLAGS)  -c $< -o $@ 

##############################################################################

setup:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)	

run:
	@mkdir -p $(RESULTS_DIR)
	make setup
	make stdp.x 
	./$(BIN_DIR)/stdp.x 
	 
##############################################################################

clean:
	rm -f *.o *.out *.mod *.dat *.txt fort*
	rm -f $(SRC_DIR)/*.mod 
	rm -f -r $(OBJ_DIR) $(BIN_DIR) 