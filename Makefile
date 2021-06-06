# Aligning BDDs: the main makefile
# Contains all key instructions to reproduce the paper
# results (figures)
#
# (c) Alexey Bochkraev, Clemson University, 2020
# abochka@clemson.edu

######################################################################
## Files and directories
PREF=./
FIGS=./figures
INST=$(PREF)/instances
LOGS=$(PREF)/run_logs
PP=./post_processing

######################################################################
## Constants
JUFL_P=0.3
JUFL_N=20
TUFL_N=20
TUFL_P=0.3
RND_N=15
RND_P=0.6
HEU_BM_NO_INST=10

## Technical constants
# how many jobs to run in parallel?
PARFLAG=2

## Key targets
figures: \
	$(FIGS)/fig_simpl_heuristics.eps \
	$(FIGS)/LB.eps \
	$(FIGS)/.opts_stats \

######################################################################
## Calculated values
HEU_BM_K=$(shell expr $(HEU_BM_NO_INST) / $(PARFLAG) + 1)

######################################################################
## Figures description

$(FIGS)/simpl_heuristics.eps: $(LOGS)/simplified_problem.csv $(PP)/fig_simpl_heuristics.R
				Rscript $(PP)/fig_simpl_heuristics.R -i $< -o $@

$(FIGS)/LB.eps: $(LOGS)/simpl_LB.csv $(PP)/fig_LBs.R
				Rscript $(PP)/fig_LBs.R -i $< -o $@

$(FIGS)/orig_obj_histograms.eps: $(LOGS)/original_problem.csv $(PP)/fig_obj_hist.R
				Rscript $(PP)/fig_obj_hist.R -i $< -o $@

$(FIGS)/orig_runtimes.eps: $(LOGS)/orig_scal.csv $(PP)/fig_scal.R
				Rscript $(PP)/fig_scal.R -i $< -o $@

$(FIGS)/no_opts.eps $(FIGS)/opts_diam.eps $(FIGS)/heuristic_simscore.eps $(FIGS)/heuristic_simscore_vs_AB_simscore.eps: .opts_stats

$(FIGS)/.opts_stats: $(LOGS)/simpl_sol_struct.csv $(PP)/figs_simpl_opt_struct.R
				Rscript $(PP)/figs_simpl_opt_struct.R --outdir $< && \
				touch $(FIGS)/.opts_stats

$(FIGS)/tUFLP_runtimes_overview.eps: $(LOGS)/tUFLP_runtimes_scal.csv $(PP)/fig_tUFLP_runtimes_scal.R
				Rscript $(PP)/fig_tUFLP_runtimes_scal.R -i $< -o $@

$(FIGS)/various_simpl_vs_min.eps: $(PP)/fig_simpl_vs_min.R $(LOGS)/heu_bm/jUFLP.csv $(LOGS)/heu_bm/tUFLP_nat.csv $(LOGS)/heu_bm/tUFLP_rnd.csv $(LOGS)/heu_bm/rnd_dia.csv
				Rscript post_processing/fig_simpl_vs_min.R --indir $(LOGS)/heu_bm --out $@

######################################################################
## Logfiles generation

$(LOGS)/heu_bm/jUFLP.csv: experiments/jUFL_hist_sizes.py
				python -m experiments.jUFL_hist_sizes -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.jUFL_hist_sizes -n $(JUFL_N) -K $(HEU_BM_K) -p $(JUFL_P) -P {} -l $(INST)/jUFLP_inst.tmp.{} ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				cat $(INST)/jUFLP_inst.tmp* > $(INST)/jUFLP_instances.json && \
				rm $(INST)/jUFLP_inst.tmp*


$(LOGS)/heu_bm/tUFLP_nat.csv: experiments/tUFL_hist_sizes_control.py
				python -m experiments.tUFL_hist_sizes_control -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.tUFL_hist_sizes_control -n $(TUFL_N) -K $(HEU_BM_K) -p $(TUFL_P) -P {} -l $(INST)/tUFLP_nat_inst.tmp.{} ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				cat $(INST)/tUFLP_nat_inst.tmp* > $(INST)/tUFLP_nat_instances.json && \
				rm $(INST)/tUFLP_nat_inst.tmp*


$(LOGS)/heu_bm/tUFLP_rnd.csv: experiments/tUFL_hist_sizes_control.py
				python -m experiments.tUFL_hist_sizes_control -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.tUFL_hist_sizes_control -R -n $(TUFL_N) -K $(HEU_BM_K) -p $(TUFL_P) -P {} -l $(INST)/tUFLP_rnd_inst.tmp.{} ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				cat $(INST)/tUFLP_rnd_inst.tmp* > $(INST)/tUFLP_rnd_instances.json && \
				rm $(INST)/tUFLP_rnd_inst.tmp*


$(LOGS)/heu_bm/rnd_dia.csv: experiments/rnd_dia_hist_sizes_control.py
				rm -rf $(INST)/bm_heu_inst && \
				mkdir -p $(INST)/bm_heu_inst && \
				python -m experiments.rnd_dia_hist_sizes_control -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.rnd_dia_hist_sizes_control -n $(RND_N) -K $(HEU_BM_K) -p $(RND_P) -P {} -l $(INST)/bm_heu_inst ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				tar czf $(INST)/bm_heu_random_diagrams.tar.gz -C $(INST)/bm_heu_inst . --remove-files


######################################################################
## Files and directories
# INST=$(PREF)/instances/raw
# ARC=./data_archive
# FIGS=./figures
# SOLVE=./solve_inst.py
# BBLOG=./log_bb.py
# SCAL=./scal_test.py
# STATS=./gen_lsizes_stats.py
# CMP_LBs=./compare_simpl_LBs.py
# PSCAL=./par_scal_test.py
# ######################################################################
# ## Numerical parameters
# PAR_SOL=8

# ### random dataset stats
# LW_n=10000
# LW_N=15
# LW_Ps=0.2 0.4 0.6 0.8

# ### dataset generation parameters:
# p=0.6# dataset generation parameter
# n=10000# number of instances
# N=15# number of variables per instance
# n_LBs=1000

# ### scalability figure
# SCAL_N=5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 #26 27 28 29 30
# SCAL_K=100
# SCAL_P=$(p)
# SCAL_R=R

# ######################################################################
# ## calculated vars/parameters
# SCAL_FILES=$(addsuffix .log, $(addprefix $(LOGS)/scal_$(SCAL_R),$(SCAL_N)))
# LW_FILES_R=$(addsuffix .log, $(addprefix $(LOGS)/lwidths_Rp,$(LW_Ps)))
# LW_FILES_N=$(addsuffix .log, $(addprefix $(LOGS)/lwidths_Np,$(LW_Ps)))

# MAX_I=$(shell expr $(PAR_SOL) - 1 )
# SOL_FILES_R = $(addsuffix .log, $(addprefix $(LOGS)/part.solved_R., $(shell seq -f %02g 0 $(MAX_I))))
# SOL_FILES_N = $(addsuffix .log, $(addprefix $(LOGS)/part.solved_N., $(shell seq -f %02g 0 $(MAX_I))))

# BB_FILES_R = $(addsuffix .log, $(addprefix $(LOGS)/part.BB_bounds_R., $(shell seq -f %02g 0 $(MAX_I))))
# BB_FILES_N = $(addsuffix .log, $(addprefix $(LOGS)/part.BB_bounds_N., $(shell seq -f %02g 0 $(MAX_I))))

# PSCAL_FILES = $(addsuffix .log, $(addprefix $(LOGS)/part.scal., $(shell seq -f %02g 0 $(MAX_I))))

# DTE=$(shell date +%F)
# TS=$(shell date +%H-%M-%S)

# MAX_I=$(shell expr $(PAR_SOL) - 1 )

# .PHONY: all figures clean-raw-inst clean-insts move-logs figures/sample_BB_tree.png

# ######################################################################
# ## Ad-hoc figures
# figures/fig_jUFL_simpl_eff.png: post_processing/fig_jUFL_simpl_eff.R run_logs/jUFL/logfile.csv
# 				Rscript post_processing/fig_jUFL_simpl_eff.R -i run_logs/jUFL/logfile.csv -o $@

# figures/fig_tUFL_simpl_eff_nat.png: post_processing/fig_jUFL_simpl_eff.R run_logs/jUFL/logfile_tUFL_nat.csv
# 				Rscript post_processing/fig_jUFL_simpl_eff.R -i run_logs/jUFL/logfile_tUFL_nat.csv -o $@ -p "Typed UFL (natural order)"

# figures/fig_tUFL_simpl_eff_rnd.png: post_processing/fig_jUFL_simpl_eff.R run_logs/jUFL/logfile_tUFL_rnd.csv
# 				Rscript post_processing/fig_jUFL_simpl_eff.R -i run_logs/jUFL/logfile_tUFL_rnd.csv -o $@ -p "Typed UFL (random order)"
# # run_logs/jUFL/logfile.csv: experiments/jUFL_hist_sizes.py
# #				./get_jUFL_hist.sh 4 # or qsub pbs/jUFL_hists.pbs

# figures/fig_simpl_vs_min.png: post_processing/fig_simpl_vs_min.R run_logs/jUFL/logfile.csv run_logs/jUFL/logfile_tUFL_nat.csv run_logs/jUFL/logfile_tUFL_rnd.csv run_logs/jUFL/logfile_rnd_dia.csv
# 				Rscript post_processing/fig_simpl_vs_min.R -o $@

# figures/fig_tUFLP_runtimes_breakdown.eps: run_logs/cUFL_runtimes.csv post_processing/fig_cUFL_runtimes_breakdown.R
# 				Rscript post_processing/fig_cUFL_runtimes_breakdown.R -i $< -o $@

# ######################################################################
# ## High-level recipes
# .SECONDARY: # keep all the intermediary logfiles (will not work o/w)

# all:

# figures/sample_BB_tree.png: ./sample_BB_tree.py
# 	python ./sample_BB_tree.py -v -V 8 -n 10 -o ./run_logs/sample_BB_tree.dot && \
# 	dot -Tpng ./run_logs/sample_BB_tree.dot > ./figures/sample_BB_tree.png

# figures: $(FIGS)/fig_sol_guessing_R.eps $(FIGS)/fig_sol_fireplace_R.eps $(FIGS)/fig_BB_gaps_R.eps $(FIGS)/LB.eps $(FIGS)/fig_sol_obj_hist_R.eps $(FIGS)/fig_sol_obj_int_R.eps $(FIGS)/guessing.tex

# #figures_R figures_N

# scalfig: $(FIGS)/fig_scal.eps

# summary_figs: $(FIGS)/fig_summary_R.eps $(FIGS)/fig_summary_N.eps
# 	touch summary_figs

# figures_R: $(FIGS)/fig_sol_guessing_R.eps $(FIGS)/fig_BB_gaps_R.eps $(FIGS)/fig_sol_fireplace_R.eps $(FIGS)/fig_sol_obj_hist_R.eps $(FIGS)/fig_sol_obj_int_R.eps
# figures_N: $(FIGS)/fig_sol_guessing_N.eps $(FIGS)/fig_BB_gaps_N.eps $(FIGS)/fig_sol_fireplace_N.eps $(FIGS)/fig_sol_obj_hist_N.eps $(FIGS)/fig_sol_obj_int_N.eps

# hists: $(FIGS)/fig_sol_obj_hist_R.eps $(FIGS)/fig_sol_obj_hist_N.eps
# 	touch hists

# ######################################################################
# ## Figure recipes
# SOL_RECIPE = Rscript $(PP)/fig_$*.R -i $< -o $@
# BB_RECIPE = Rscript $(PP)/fig_BB_$*.R -i $< -o $@

# $(FIGS)/fig_sol_%_R.eps: $(LOGS)/solved_R.log $(PP)/fig_%.R
# 	$(SOL_RECIPE)

# $(FIGS)/fig_sol_%_N.eps: $(LOGS)/solved_N.log $(PP)/fig_%.R
# 	$(SOL_RECIPE)

# $(FIGS)/fig_BB_%_R.eps: $(LOGS)/BB_bounds_R.log $(PP)/fig_BB_%.R
# 	$(BB_RECIPE)

# $(FIGS)/fig_BB_%_N.eps: $(LOGS)/BB_bounds_N.log $(PP)/fig_BB_%.R
# 	$(BB_RECIPE)

# $(FIGS)/fig_scal.eps: $(LOGS)/scal_par.log $(PP)/fig_scal.R
# 	Rscript $(PP)/fig_scal.R -i $< -o $@

# $(FIGS)/fig_summary_R.eps: DS_FLAG=R
# $(FIGS)/fig_summary_N.eps: DS_FLAG=N

# $(FIGS)/fig_summary_%.eps: $(LOGS)/lwidths_%.log $(PP)/fig_summary.R
# 	Rscript $(PP)/fig_summary.R -i $< -o $@

# $(FIGS)/guessing.tex: $(LOGS)/solved_R.log $(PP)/tab_guessing.R
# 	Rscript $(PP)/tab_guessing.R -i $< -o $@

# ######################################################################
# ## scalability figure -- parallel implementation
# $(LOGS)/scal_par.log: $(PSCAL_FILES)
# 	python $(PSCAL) --header > $@ && \
# 	tail -qn +2 $(PSCAL_FILES) >> $@ && \
# 	rm -rf $(PSCAL_FILES)

# $(LOGS)/part.scal.%.log: $(INST)/scal/instances.list
# 	python $(PSCAL) -i $<.$* -o $@ -d $(INST)/scal/

# $(INST)/scal/instances.list: ./gen_BDD_pair.py
# 	@echo Preparing instances for the scalability set...
# 	i=0 && \
# 	rm -rf $(INST)/scal && \
# 	mkdir -p $(INST)/scal && \
# 	for S in $(SCAL_N); do \
# 		python ./gen_BDD_pair.py -s $$(( ${SCAL_K} * $${i} )) -n $(SCAL_K) -v $$S -p $p -RU $(INST)/scal/ > /dev/null &&\
# 		i=$$(( $${i} + 1 ));\
# 	done && \
# 	ls $(INST)/scal | grep -Po "A\\K[^\\.]*" > $@ && \
# 	shuf $@ > $@.shuffled && \
# 	split -d -nl/$(PAR_SOL) $@.shuffled $(INST)/scal/instances.list. && \
# 	rm $@.shuffled

# # was: ls $(INST)/$*/A*.bdd | grep -Po "$(INST)/$*/A\\K[^\\.]*" > $@ &&\
# # basename -a `ls $(INST)/$*/A*.bdd` | sed -n -e's/^A//p' | sed -n -e's/\.bdd//p' > $@ && \
# #
# ######################################################################
# ## Main calculations (creating .log-files)
# .SECONDEXPANSION:
# $(LOGS)/solved_%.log: $$(SOL_FILES_%)
# 	python $(SOLVE) --header > $@ && \
# 	tail -qn +2 $(SOL_FILES_$*) >> $@


# $(LOGS)/part.solved_R.%.log: $(INST)/R/instances.list $(SOLVE)
# 	python $(SOLVE) -i $<.$* -o $@ -d $(INST)/R/
# ######################################################################
# $(FIGS)/LB.eps: $(LOGS)/LBs.log $(PP)/fig_LBs.R
# 	Rscript $(PP)/fig_LBs.R -i $< -o $@

# $(LOGS)/LBs.log: $(CMP_LBs)
# 	mkdir -p $(INST)/LBs && \
# 	python ./gen_BDD_pair.py -n $(n_LBs) -v $N -p $p -RU $(INST)/LBs/ > $(LOGS)/$(DTE)_gen_$(N)var_LBs.log && \
# 	ls $(INST)/LBs | grep -Po "A\\K[^\\.]*" > $(INST)/LBs/instances.list && \
# 	python $(CMP_LBs) -d $(INST)/LBs/ -l $(INST)/LBs/instances.list > $@

# ## Generating and solving instances
# $(INST)/%/instances.list:
# 	@echo Preparing the $*-instances dataset...
# 	mkdir -p $(INST)/$* && \
# 	if [ -f $(ARC)/dataset_$*.tar.gz ]; then tar -zxmf $(ARC)/dataset_$*.tar.gz; \
# 	else \
# 	python ./gen_BDD_pair.py -n $n -v $N -p $p -$*U $(INST)/$*/ > $(LOGS)/$(DTE)_gen_$(N)var_R.log; fi && \
# 	ls $(INST)/$* | grep -Po "A\\K[^\\.]*" > $@ && \
# 	split -d -nl/$(PAR_SOL) $@ $(INST)/$*/instances.list.

# # was: ls $(INST)/$*/A*.bdd | grep -Po "$(INST)/$*/A\\K[^\\.]*" > $@ &&\
# # basename -a `ls $(INST)/$*/A*.bdd` | sed -n -e's/^A//p' | sed -n -e's/\.bdd//p' > $@ && \
# #
# ######################################################################
# ## Main calculations (creating .log-files)
# .SECONDEXPANSION:
# $(LOGS)/solved_%.log: $$(SOL_FILES_%)
# 	python $(SOLVE) --header > $@ && \
# 	tail -qn +2 $(SOL_FILES_$*) >> $@


# $(LOGS)/part.solved_R.%.log: $(INST)/R/instances.list $(SOLVE)
# 	python $(SOLVE) -i $<.$* -o $@ -d $(INST)/R/

# $(LOGS)/part.solved_N.%.log: $(INST)/N/instances.list $(SOLVE)
# 	python $(SOLVE) -i $<.$* -o $@ -d $(INST)/N/

# $(LOGS)/BB_bounds_%.log: $$(BB_FILES_%)
# 	python $(BBLOG) --header > $@ && \
# 	tail -qn +2 $(BB_FILES_$*) >> $@

# $(LOGS)/part.BB_bounds_R.%.log: $(INST)/R/instances.list $(SOLVE)
# 	python $(BBLOG) -i $<.$* -o $@ -d $(INST)/R/

# $(LOGS)/part.BB_bounds_N.%.log: $(INST)/N/instances.list $(SOLVE)
# 	python $(BBLOG) -i $<.$* -o $@ -d $(INST)/N/

# ## scalability experiment
# $(LOGS)/scal_$(SCAL_R)%.log: $(SCAL)
# 	mkdir -p $(INST)/sc$* && \
# 	python ./gen_BDD_pair.py -n $(SCAL_K) -v $* -p $(SCAL_P) -U$(SCAL_R) $(INST)/sc$* > $(INST)/sc$*/gen.log && \
# 	basename -a -s .bdd $(INST)/sc$*/A*.bdd | sed 's/A//' > $(INST)/sc$*/inst.list && \
# 	python $(SCAL) -l $(INST)/sc$*/inst.list -d $(INST)/sc$*/ > $@

# $(LOGS)/scalability.log: $(SCAL) $(SCAL_FILES)
# 	python $(SCAL) --header > $@ && \
# 	cat $(SCAL_FILES) >> $@ && \
# 	tar --remove-files -czf $(ARC)/scalability_p$(SCAL_P)$(SCAL_R)x$(SCAL_K).tar.gz $(INST)/sc* && \
# 	rm -rf $(INST)/sc$* && \
# 	rm $(LOGS)/scal_N*.log
# ## end of scalability experiment

# ## random dataset stats
# GEN_LW_LOGS =	mkdir -p $(INST)/ds_stats/$(DS_FLAG)p$* && \
# 	python ./gen_BDD_pair.py -n $(LW_n) -v $(LW_N) -p $* -$(DS_FLAG)U $(INST)/ds_stats/$(DS_FLAG)p$* > /dev/null && \
# 	ls $(INST)/ds_stats/$(DS_FLAG)p$*/*.bdd > $(INST)/ds_stats/$(DS_FLAG)p$*/lwidths_$(DS_FLAG)p$*.list && \
# 	python $(STATS) -s $* $(INST)/ds_stats/$(DS_FLAG)p$*/lwidths_$(DS_FLAG)p$*.list > $@ && \
# 	rm -rf $(INST)/ds_stats/$(DS_FLAG)p$*

# $(LOGS)/lwidths_Np%.log: $(STATS)
# 	$(GEN_LW_LOGS)

# $(LOGS)/lwidths_Rp%.log: $(STATS)
# 	$(GEN_LW_LOGS)

# .SECONDEXPANSION:
# $(LOGS)/lwidths_%.log: $$(LW_FILES_%)
# 	ls $(LOGS)/lwidths_$(DS_FLAG)p*.log && \
# 	tail -qn +2 $(LOGS)/lwidths_$(DS_FLAG)p*.log > $@

# ## end of random dataset stats

# ######################################################################
# # auxiliary recipes

# prep_dirs:
# 	mkdir -p $(INST)/R $(INST)/N
# 	mkdir -p $(LOGS)

# check_src:
# 	egrep -nr --color 'TODO|FIXME|BUG|NOTE'

# install_R_pkg:
# 	Rscript ./aux/R_install_packages.R

# move-logs:
# 	@echo moving logs away from $(LOGS) to $(ARC)
# 	rm -f $(LOGS)/part.* && \
# 	tar --remove-files -czf $(ARC)/$(DTE)_$(TS)_logs.tar.gz -C $(LOGS) .

# # clean recipes
# clean-raw-inst:
# 	@echo Cleaning up raw instances in $(INST)...
# 	rm -f $(INST)/R/*.bdd
# 	rm -f $(INST)/R/*.list*
# 	rm -f $(INST)/N/*.bdd
# 	rm -f $(INST)/N/*.list*
# 	rm -rf $(INST)/ds_stats
# 	rm -f $(INST)/*.bdd

# clean-archive:
# 	@echo Cleaning up instance archives...
# 	rm -f $(ARC)/*.tar.gz

# clean-logs:
# 	@echo Cleaning up log files...
# 	rm -f $(LOGS)/*.log

# clean-figures:
# 	@echo Cleaning up figures...
# 	rm -f $(FIGS)/*.eps

# clean: clean-raw-inst clean-logs clean-figures
# clean-tmp: clean-raw-inst clean-logs

# test: $(LOGS)/BB.test
# 	@echo Running tests...

# $(LOGS)/BB.test: $(LOGS)/BB_bounds_R.log $(LOGS)/solved_R.log tests/BB_log_correct.R
# 	Rscript tests/BB_log_correct.R -b $(LOGS)/BB_bounds_R.log -s $(LOGS)/solved_R.log > $(LOGS)/BB.test

# qtest:
# 	python -m pytest tests/UFL_test.py
