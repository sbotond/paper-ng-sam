#
# Simulations Makefile
#

.PHONY:	seq_sim dil_sim

# General parameters:
LSF_QUEUE		= research-rh6

SIMNGS_RUNFILE  =$(BASE)/dat/s_8_4x.runfile
READ_LENGTH     =101
INSERT_SIZE     =400

INIT_POPSIZE    = 5000
PCR_EFFICIENCY  = 0.75
CYCLES_MUT      = 20
DILF_MUT        = 70000
CYCLES_CLEAN    = 20
DILF_CLEAN      = 16000000
CYCLES_FINAL    = 30
SAMPLE_SIZE     = 40
VMIN_CTGL       = 400
VKMER_LENGTH    = 90
VMAX_DIV        = 0.1
TOTAL_COV       = 4000

# seq_sim specific parameters:
SEQ_OUT_DIR		=$(BASE)/seq_sim/out
SEQ_TARGET_DIR	= $(BASE)/seq_sim/targets
NR_REPS			= 5				# Number of replications for a given unit length/number combination.
UNIT_NR_RANGE	= 4:100:4		# Unit number range and step size.
UNIT_LEN_RANGE	= 4:4000:5		# Unit length range and step size.
MIN_TLEN		= 500			# Minimum target sequence length.
MAX_TLEN		= 30000			# Maximum target sequence length.

SEQ_SIM_PARAMS  = "-n $(MUT_MODEL_FILE) -b $(BL_SCALER_FILE) -i $(INIT_POPSIZE) -e $(PCR_EFFICIENCY) -cm $(CYCLES_MUT) -dm $(DILF_MUT)  \
-cc $(CYCLES_CLEAN) -dc $(DILF_CLEAN) -cf $(CYCLES_FINAL) -ss $(SAMPLE_SIZE) -vm $(VMIN_CTGL) -vk $(VKMER_LENGTH) -vi $(INSERT_SIZE) \
-vd $(VMAX_DIV) -P $(BIN) -S $(SIMNGS_RUNFILE) -L $(READ_LENGTH) -I $(INSERT_SIZE) -t $(TOTAL_COV)"

# Simulate NG-SAM experiments on random targets with varying repetitive structure:
seq_sim:
	@bin/run_seq_sim -m $(MIN_TLEN) -M $(MAX_TLEN) -Q $(LSF_QUEUE) -X '$(SEQ_SIM_PARAMS)' -R $(RUN_DIR) -T $(SEQ_TARGET_DIR) -n $(NR_REPS) -u $(UNIT_NR_RANGE) -l $(UNIT_LEN_RANGE) -o $(SEQ_OUT_DIR)

# Visualise the results of seq_sim:
plot_seq_res:
	@bin/plot_seq_res -i $(SEQ_OUT_DIR) -r $(REP_DIR)/seq_sim.pdf -g 27

# dil_sim spacific parameters:
DIL_OUT_DIR		=$(BASE)/dil_sim/out
DIL_TARGET_DIR	= $(BASE)/dil_sim/targets
NR_REPS			= 5							# Number of replications for a given dm/dc combination.
D1_RANGE		= 17000:280000:20000		# First dilution range.
D2_RANGE		= 2000000:128000000:200000	# Second dilution range.
DIL_TARGET		= $(BASE)/dat/eater_root.fas

DIL_SIM_PARAMS  = "-n $(MUT_MODEL_FILE) -b $(BL_SCALER_FILE) -i $(INIT_POPSIZE) -e $(PCR_EFFICIENCY) -cm $(CYCLES_MUT) \
-cc $(CYCLES_CLEAN) -cf $(CYCLES_FINAL) -ss $(SAMPLE_SIZE) -vm $(VMIN_CTGL) -vk $(VKMER_LENGTH) -vi $(INSERT_SIZE) \
-vd $(VMAX_DIV) -P $(BIN) -S $(SIMNGS_RUNFILE) -L $(READ_LENGTH) -I $(INSERT_SIZE) -t $(TOTAL_COV)"

# Simulate NG-SAM experiments under a range of dilution factors.
dil_sim:
	@bin/run_dil_sim -t $(DIL_TARGET) -m $(MIN_TLEN) -M $(MAX_TLEN) -Q $(LSF_QUEUE) -X '$(DIL_SIM_PARAMS)' -R $(RUN_DIR) -n $(NR_REPS) -o $(DIL_OUT_DIR) -d1 $(D1_RANGE) -d2 $(D2_RANGE)

# Visualise the results of dil_sim:
plot_dil_res:
	@bin/plot_dil_res -i $(DIL_OUT_DIR) -r $(REP_DIR)/dil_sim.pdf -g 13

