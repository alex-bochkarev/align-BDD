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

STATS=python -m experiments.gen_lsizes_stats
######################################################################
## Constants
JUFL_P=0.3
JUFL_N=20
TUFL_N=20
TUFL_P=0.3
RND_N=15
RND_P=0.6
HEU_BM_NO_INST=2000

# Number of vars of the main run experiment
ORIG_N=15
ORIG_K_TOTAL=10000

STAT_K_TOTAL = 1000
STAT_PS = 0.2 0.5 0.8

# Scaling figure: how many instances per size?
SCAL_K=200

## Technical constants
# how many jobs to run in parallel?
PARFLAG=2

## Key targets
figures: \
	$(FIGS)/guessing.tex \
	$(FIGS)/simpl_heuristics.eps \
	$(FIGS)/LB.eps \
	$(FIGS)/orig_obj_histograms.eps \
	$(FIGS)/orig_runtimes.eps \
	$(FIGS)/.opts_stats \
	$(FIGS)/tUFLP_runtimes_overview.eps \
	$(FIGS)/various_simpl_vs_min.eps \
	$(FIGS)/orig_lwidth_stats.eps \
	$(FIGS)/tUFLP_runtimes_breakdown.eps

######################################################################
## Calculated values
HEU_BM_K=$(shell expr $(HEU_BM_NO_INST) / $(PARFLAG) + 1)
STAT_K=$(shell expr $(STAT_K_TOTAL) / $(PARFLAG) + 1)
ORIG_K=$(shell expr $(ORIG_K_TOTAL) / $(PARFLAG) + 1)


######################################################################
## Figures description

## The table (simplified problem heuristics)

$(FIGS)/guessing.tex: $(LOGS)/main_rnd_run.csv $(PP)/tab_guessing.R
				Rscript $(PP)/tab_guessing.R -i $< -o $@


## Comparison of the heuristics to the *simplified* problem

$(FIGS)/simpl_heuristics.eps: $(LOGS)/main_rnd_run.csv $(PP)/fig_simpl_heuristics.R
				Rscript $(PP)/fig_simpl_heuristics.R -i $< -o $@


## Different approaches to lower bound
$(FIGS)/LB.eps: $(LOGS)/simpl_LB.csv $(PP)/fig_LBs.R
				Rscript $(PP)/fig_LBs.R -i $< -o $@

$(LOGS)/simpl_LB.csv: $(INST)/orig_problem.list experiments/compare_simpl_LBs.py
				rm -f $<.tmp.* && \
				split -d -nl/$(PARFLAG) $< $<.LBs.tmp. && \
				parallel -j $(PARFLAG) python -m experiments.compare_simpl_LBs -l {} ">" $@.tmp.{#} -d $(INST)/orig_problem/ ::: $$(ls $<.LBs.tmp.*) && \
				python -m experiments.compare_simpl_LBs --header > $@ && \
				cat $@.tmp.* >> $@ && \
				rm $@.tmp.* && \
				rm $<.LBs.tmp.* && \
				if [ "$(PREF)" != "./" ]; then cp $@ ./run_logs/simpl_LB.csv; fi



## Histograms of the original objective value
$(FIGS)/orig_obj_histograms.eps: $(LOGS)/main_rnd_run.csv $(PP)/fig_obj_hist.R
				Rscript $(PP)/fig_obj_hist.R -i $< -o $@

$(LOGS)/main_rnd_run.csv: $(INST)/orig_problem.list
				rm -f $<.tmp.* && \
				split -d -nl/$(PARFLAG) $< $<.tmp. && \
				parallel -j $(PARFLAG) python -m solve_inst -i {} -o $@.tmp.{#} -d $(INST)/orig_problem/ ::: $$(ls $<.tmp.*) && \
				python -m solve_inst --header > $@ && \
				cat $@.tmp.* >> $@ && \
				rm $@.tmp.* && \
				rm $<.tmp.* && \
				if [ "$(PREF)" != "./" ]; then \
					cp $@ ./run_logs/main_rnd_run.csv ; \
					cp -r $(INST)/orig_problem ./instances/ ; \
				fi

$(INST)/orig_problem.list: gen_BDD_pair.py
				rm -rf $(INST)/orig_problem && \
				mkdir $(INST)/orig_problem && \
				parallel -j $(PARFLAG) 'python -m gen_BDD_pair -s $$(( $(ORIG_K) * {#} - $(ORIG_K)))' -K $(ORIG_K) -v $(ORIG_N) -p $(RND_P) -R -U $(INST)/orig_problem/ --quiet ::: $(shell seq 1 $(PARFLAG)) && \
				ls $(INST)/orig_problem | grep -Po "A\\K[^\\.]*" > $@ && \
				if [ "$(PREF)" != "./" ]; then \
					cp -r $(INST)/orig_problem ./instances/; \
					cp $@ ./instances/; \
				fi

## Original problem (rnd) scaling figure
$(FIGS)/orig_runtimes.eps: $(LOGS)/orig_scal.csv $(PP)/fig_scal.R
				Rscript $(PP)/fig_scal.R -i $< -o $@


$(LOGS)/orig_scal.csv: $(INST)/scal/instances.list experiments/par_scal_test.py
				rm -f $(INST)/scal/instances.list.* && \
				split -d -nl/$(PARFLAG) $< $(INST)/scal/instances.list. && \
				parallel -j $(PARFLAG) python -m experiments.par_scal_test -i {} -o $@.tmp.{#} -d $(INST)/scal/ ::: $$(ls $(INST)/scal/instances.list.*) && \
				python -m experiments.par_scal_test --header > $@ && \
				cat $@.tmp.* >> $@ && \
				rm $@.tmp.* && \
				tar czf $(INST)/orig_scal.tar.gz -C $(INST)/scal . --remove-files && \
				if [ "$(PREF)" != "./" ]; then \
					cp $(INST)/orig_scal.tar.gz ./instances/; \
					cp $@ ./run_logs/; \
					rm -f ./instances/scal; \
				fi

$(INST)/scal/instances.list: gen_BDD_pair.py
				mkdir -p $(INST)/scal && \
				parallel -j $(PARFLAG) 'python -m gen_BDD_pair -s $$(( $(SCAL_K) * {#} - $(SCAL_K)))' -K $(SCAL_K) -v {} -p $(RND_P) -R -U $(INST)/scal/ --quiet ::: $(shell seq 5 25 | shuf) && \
				ls $(INST)/scal | grep -Po "A\\K[^\\.]*" > $@ && \
				if [ "$(PREF)" != "./" ]; then \
					cp -r $(INST)/scal ./instances/; \
				fi


## Optima stats for the simplified problem
$(FIGS)/no_opts.eps $(FIGS)/opts_diam.eps $(FIGS)/heuristic_simscore.eps $(FIGS)/heuristic_simscore_vs_AB_simscore.eps: $(FIGS)/.opts_stats

$(FIGS)/.opts_stats: $(LOGS)/simpl_sol_struct.csv $(PP)/figs_simpl_opt_struct.R
				Rscript $(PP)/figs_simpl_opt_struct.R -i $< --outdir $(FIGS) && \
				touch $(FIGS)/.opts_stats

$(LOGS)/simpl_sol_struct.csv: experiments/heu_sol_struct.py
				mkdir -p $(INST)/simpl_sol_struct && \
				python -m experiments.heu_sol_struct -k 1000 -n 7 -l $(INST)/simpl_sol_struct > $@ && \
				tar czf $(INST)/heu_sol_struct.tar.gz -C $(INST)/simpl_sol_struct . --remove-files ;  \
				if [ "$(PREF)" != "./" ]; then \
					cp $@ ./run_logs/; \
					mv $(INST)/heu_sol_struct.tar.gz ./instances/ ; \
				fi

## t-UFLP runtimes overview (scaling)
$(FIGS)/tUFLP_runtimes_overview.eps: $(LOGS)/tUFLP_runtimes_scal.csv $(PP)/fig_tUFLP_runtimes_scal.R
				Rscript $(PP)/fig_tUFLP_runtimes_scal.R -i $< -o $@

$(LOGS)/tUFLP_runtimes_scal.csv: experiments/tUFLP_runtimes.py
				parallel -j $(PARFLAG) python -m experiments.tUFLP_runtimes -K 200 -n {} --with-gsifts -l $(INST)/tUFLP_scal_inst.json.tmp.{} ">" $@.tmp.{} ::: $(shell seq 5 15 | shuf) && \
				python -m experiments.tUFLP_runtimes -H > $@ && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				:> $(INST)/tUFLP_scal_inst.json && \
				cat $(INST)/tUFLP_scal_inst.json.tmp.* >> $(INST)/tUFLP_scal_inst.json && rm $(INST)/tUFLP_scal_inst.json.tmp.* && \
				if [ "$(PREF)" != "./" ]; then \
					cp $@ ./run_logs/; \
					cp $(INST)/tUFLP_scal_inst.json ./instances/ ; \
				fi

$(FIGS)/various_simpl_vs_min.eps: $(PP)/fig_simpl_vs_min.R $(LOGS)/heu_bm/jUFLP.csv $(LOGS)/heu_bm/tUFLP_nat.csv $(LOGS)/heu_bm/tUFLP_rnd.csv $(LOGS)/heu_bm/rnd_dia.csv
				Rscript post_processing/fig_simpl_vs_min.R --indir $(LOGS)/heu_bm --out $@

######################################################################
## Logfiles generation

$(LOGS)/heu_bm/jUFLP.csv: experiments/jUFL_hist_sizes.py
				mkdir -p $(LOGS)/heu_bm && \
				python -m experiments.jUFL_hist_sizes -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.jUFL_hist_sizes -n $(JUFL_N) -K $(HEU_BM_K) -p $(JUFL_P) -P {} -l $(INST)/jUFLP_inst.tmp.{} ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				cat $(INST)/jUFLP_inst.tmp* > $(INST)/jUFLP_instances.json && \
				rm $(INST)/jUFLP_inst.tmp* && \
				if [ "$(PREF)" != "./" ]; then \
					mkdir -p ./run_logs/heu_bm ; \
					cp $@ ./run_logs/heu_bm/; \
					cp $(INST)/jUFLP_instances.json ./instances/ ; \
				fi


$(LOGS)/heu_bm/tUFLP_nat.csv: experiments/tUFL_hist_sizes_control.py
				mkdir -p $(LOGS)/heu_bm && \
				python -m experiments.tUFL_hist_sizes_control -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.tUFL_hist_sizes_control -n $(TUFL_N) -K $(HEU_BM_K) -p $(TUFL_P) -P {} -l $(INST)/tUFLP_nat_inst.tmp.{} ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				cat $(INST)/tUFLP_nat_inst.tmp* > $(INST)/tUFLP_nat_instances.json && \
				rm $(INST)/tUFLP_nat_inst.tmp* && \
				if [ "$(PREF)" != "./" ]; then \
					mkdir -p ./run_logs/heu_bm ; \
					cp $@ ./run_logs/heu_bm/ ; \
					cp $(INST)/tUFLP_nat_instances.json ./instances/ ; \
				fi


$(LOGS)/heu_bm/tUFLP_rnd.csv: experiments/tUFL_hist_sizes_control.py
				mkdir -p $(LOGS)/heu_bm && \
				python -m experiments.tUFL_hist_sizes_control -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.tUFL_hist_sizes_control -R -n $(TUFL_N) -K $(HEU_BM_K) -p $(TUFL_P) -P {} -l $(INST)/tUFLP_rnd_inst.tmp.{} ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				cat $(INST)/tUFLP_rnd_inst.tmp* > $(INST)/tUFLP_rnd_instances.json && \
				rm $(INST)/tUFLP_rnd_inst.tmp* && \
				if [ "$(PREF)" != "./" ]; then \
					mkdir -p ./run_logs/heu_bm/ ; \
					cp $@ ./run_logs/heu_bm/ ; \
					cp $(INST)/tUFLP_rnd_instances.json ./instances/ ; \
				fi


$(LOGS)/heu_bm/rnd_dia.csv: experiments/rnd_dia_hist_sizes_control.py
				mkdir -p $(LOGS)/heu_bm && \
				rm -rf $(INST)/bm_heu_inst && \
				mkdir -p $(INST)/bm_heu_inst && \
				python -m experiments.rnd_dia_hist_sizes_control -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.rnd_dia_hist_sizes_control -n $(RND_N) -K $(HEU_BM_K) -p $(RND_P) -P {} -l $(INST)/bm_heu_inst ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				tar czf $(INST)/bm_heu_random_diagrams.tar.gz -C $(INST) ./bm_heu_inst --remove-files && \
				if [ "$(PREF)" != "./" ]; then \
					mkdir -p ./run_logs/heu_bm ; \
					cp $@ ./run_logs/heu_bm/ ; \
					cp $(INST)/bm_heu_random_diagrams.tar.gz ./instances/ ; \
				fi


## Random dataset stats
$(FIGS)/orig_lwidth_stats.eps: $(LOGS)/lwidths.csv $(PP)/fig_summary.R
				Rscript $(PP)/fig_summary.R -i $< -o $@

$(LOGS)/lwidths.csv: gen_BDD_pair.py
				mkdir -p $(INST)/orig_stats && \
				parallel mkdir -p $(INST)/orig_stats/{} ::: $(STAT_PS) && \
				parallel -j $(PARFLAG) python -m gen_BDD_pair -v $(RND_N) -K $(STAT_K) -p {1} -R -U $(INST)/orig_stats/{1} '-s $$(( $(STAT_K) * {#} - $(STAT_K)))' --quiet ::: $(STAT_PS) ::: $(shell seq 1 $(PARFLAG)) && \
				parallel find $(INST)/orig_stats/{}/*.bdd ">" $(INST)/orig_stats/{}/inst.list ::: $(STAT_PS) && \
				parallel $(STATS) -s {} $(INST)/orig_stats/{}/inst.list ">" $(LOGS)/lwidths.tmp.{} ::: $(STAT_PS) && \
				cat $(LOGS)/lwidths.tmp.* >> $@ && rm $(LOGS)/lwidths.tmp.* && \
				tar czf $(INST)/orig_stats.tar.gz -C $(INST)/orig_stats . --remove-files && \
				if [ "$(PREF)" != "./" ]; then \
					cp $@ ./run_logs/ ; \
					cp $(INST)/orig_stats.tar.gz ./instances/ ; \
				fi


## Breakdown of runtimes (CPP)

$(FIGS)/tUFLP_runtimes_breakdown.eps: $(LOGS)/tUFLP_runtimes.csv $(PP)/fig_tUFLP_runtimes_breakdown.R
				Rscript $(PP)/fig_tUFLP_runtimes_breakdown.R -i $< -o $@

$(LOGS)/tUFLP_runtimes.csv: experiments/tUFLP_runtimes.py
				python -m experiments.tUFLP_runtimes -H > $@ && \
				python -m experiments.tUFLP_runtimes -K 1000 -n 20 -l $(INST)/tUFLP_steps_breakdown.json >> $@ && \
				if [ "$(PREF)" != "./" ]; then \
					cp $@ ./run_logs/ ; \
					cp $(INST)/tUFLP_steps_breakdown.json ./instances/ ; \
				fi

## clean-up
# ######################################################################
# # auxiliary recipes

save-orig-instances:
				tar czf $(INST)/orig_problem.tar.gz -C $(INST)/orig_problem . --remove-files

prep_dirs:
				cp --preserve=timestamps -r ./run_logs $(LOGS)
				cp --preserve=timestamps -r ./instances $(INST)

clean:
				rm -rf $(INST)/*
				rm -rf $(LOGS)/*
				touch ./clean

fresh: clean
				rm -rf $(FIGS)/*.eps
				mkdir -p ./run_logs
				mkdir -p ./instances
				rm ./clean

check_src:
	egrep -nr --color 'TODO|FIXME|BUG|NOTE'

# install_R_pkg:
# 	Rscript ./aux/R_install_packages.R

cross-checks:
				python -m pytest tests/BDD_test.py
				python -m pytest BDD.py --doctest-modules
				python -m pytest tests/varseq_test.py
				python -m pytest varseq.py --doctest-modules
				python -m pytest tests/BB_search_test.py
				python -m pytest BB_search.py --doctest-modules
				python -m pytest tests/UFLP_test.py
				python -m pytest tests/tUFLP_test.py
