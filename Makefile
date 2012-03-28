.PHONY: com
#
# Main Makefile
#

# Directories and general parameters:
BASE			=$(shell pwd)
BIN				=$(BASE)/bin
REP_DIR			=$(BASE)/reports
BL_SCALER_FILE	=$(BASE)/dat/bl_scaler.txt
MUT_MODEL_FILE	=$(BASE)/dat/mutation_model.tab
RUN_DIR			=/tmp

# Calibration parameters:
DESIRED_MUTRATE	= 0.05	# Target mutation rate.
NR_CAL_CYCLES	= 10	# Number of simulated calibration cycles.
CAL_SIM			= 5000	# Number of simulated calibration experiments.

# Commit all changes:
com:
	git commit -a

# Push to origin:
push:
	git push --all

# Run a test simulation:
t:
	@bin/sim_exp -vm 300 -vk 80 -vi 400 -vd 0.1 -n $(MUT_MODEL_FILE) -b $(BL_SCALER_FILE) -f dat/dmel_eater.fas -i 5000 -e 0.75 -cm 15 -dm 4900 -cc 20 -dc 16000000 -cf 30 -P $(BIN) -R ./ -S $(SIMNGS_RUNFILE) -L $(READ_LENGTH) -I $(INSERT_SIZE) -t 4000
	@cat dmel_eater.out; rm dmel_eater.out

# Calculate branch length scaling factor:	
$(BL_SCALER_FILE): $(MUT_MODEL_FILE)
	bin/calibrate_mut -q $(MUT_MODEL_FILE) -o $(BL_SCALER_FILE) -g $(CAL_SIM) -c $(NR_CAL_CYCLES)  -m $(DESIRED_MUTRATE) -f dat/MH22.fas -r $(REP_DIR)/calibration_report.pdf -P $(BIN) 

calibration: $(BL_SCALER_FILE)
	@echo -n

include simulations.mk
