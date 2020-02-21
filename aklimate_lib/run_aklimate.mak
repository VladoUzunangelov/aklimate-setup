# chrisw

THIS_DIR=$(shell pwd)

THIS_DIR_NAME=$(shell pwd | cut.pl -d "/" -f -1)

DATA_DIR=data
LIB_DIR=lib

COMBINED_MATRIX_FILE=$(DATA_DIR)/combined_matrix.tsv
CV_SETS_FILE=$(DATA_DIR)/cv_folds.tsv

TARGETS= \
	feature_set_weights.tsv \
	bal_acc_subtype_50_cutoff.png \
	feature_importance_stats.tsv \
	datatype_stacked_bar_plots.png \
	sickle_plot.png \
	cv_test_sample_predictions_full.tsv \
	datatypes.tsv \
	samples.tsv \
	labels.tsv \
	\

SUMMARY_TARGETS= \
	bacc.tsv \
	bacc_stats.tsv \
	collect_predictions \
	aklimate_sample_predictions.tsv \
	sickle_plot.png \
	datatype_stacked_bar_plots.png \
	feature_importance_stats.tsv \
	bal_acc_subtype_50_cutoff.png \
	feature_set_weights.tsv \

REDUCED_MODELS_CUTOFFS= \
	$(shell find ./models/ -type f -name "*rf_reduced_model_predictions.RData" | cut.pl -d "cutoff_" -f 2 | cut -d "_" -f 1| sort.pl | uniq)

test:

feature_set_weights.tsv:
	rm -f 1.tmp ;
	\
	for fold in $(shell ls -1 ./models | grep "_junkle_final_model.RData" | cut.pl -d "_junkle_final_model.RData" -f 1 ) ; do \
		echo $${fold} ; \
		\
		Rscript get_aklimate_feature_set_weights.R \
			./models/$${fold}_junkle_final_model.RData \
			$${fold}_feature_set_weights.tmp \
		; \
		\
		\
		if [ -f "$${fold}_feature_set_weights.tmp_1" ] ; then \
			echo "multiclass model detected" ; \
			\
			for i in `ls -1 *feature_set_weights.tmp* | cut.pl -d "_" -f -1` ; do \
				echo $${fold}_$${i} ; \
				\
				paste.pl "$${fold}_$${i}" $${fold}_feature_set_weights.tmp_$${i} \
				| tail -n +2 \
				>> 1.tmp ; \
				\
				rm -f $${fold}_feature_set_weights.tmp_$${i} ; \
				\
			done ; \
		\
		\
		else \
			echo "defaultint to handling feature set weights from binary model" ; \
			\
			paste.pl "$${fold}" $${fold}_feature_set_weights.tmp \
			| tail -n +2 \
			>> 1.tmp ; \
		\
		\
		fi ; \
		rm -f $${fold}_feature_set_weights.tmp ; \
		\
	done ;
	\
	cut -f 2 1.tmp \
	| cut.pl -d "_" -f 1--6 \
	> 2.tmp ;
	\
	paste.pl 2.tmp 1.tmp \
	> 3.tmp ;
	\
	join.pl -r -1 1 -2 1 -o "__NONE__" 3.tmp ./p_store_files/pathway_name_mapping.tsv \
	> 4.tmp ;
	\
	mv 4.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp 3.tmp 4.tmp ;
	\

bal_acc_subtype_50_cutoff.png:
	Rscript find_balance_accuracy_subcohort.R ;
	\

feature_importance_stats.tsv:
	rm -f 1.tmp ;
	\
	for file in `find ./models -name "*_aklimate_multiclass_feature_importance.tab" ` ; do \
		echo $${file} ; \
		\
		tail -n +2 $${file} \
		>> 1.tmp ; \
		\
	done ;
	\
	expand.pl 1.tmp \
	> 2.tmp ;
	\
	cat 2.tmp \
	| row_stats.pl -h 0 -k 0 -allstats \
	> 3.tmp ;
	\
	cat 3.tmp \
	| sort.pl -h 1 -k 9 -r \
	> 4.tmp ;
	\
	mv 4.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp 3.tmp 4.tmp ;
	\
	\

datatype_stacked_bar_plots.png:
	rm -f colorfeature_proportion_barplot.png ;
	\
	Rscript feature_type_stacked_barplot.R ;
	\
	mv colorfeature_proportion_barplot.png $@ ;
	\
	rm -f colorfeature_proportion_barplot.png ;
	\

sickle_plot.png:
	Rscript --vanilla ./bacc_sickle_plot.R \
		./models/balanced_accuracy_reduced_feature_sets.tab \
		1.tmp \
	;
	mv 1.tmp $@ ;
	\
	rm -f 1.tmp ;
	\

aklimate_sample_predictions.tsv:
	cat $(shell find ./collect_predictions -type f -name "cv_test_sample_predictions_*" | head -n 1) \
	> 1.tmp ;
	\
	grep -v '^#' 1.tmp \
	| cut -f 1-5 \
	| sed -e 's/	/___/g' \
	> 2.tmp ;
	\
	\
	rm -f comments.tmp ;
	\
	for file in $(shell find ./collect_predictions -type f -name "cv_test_sample_predictions_*") ; do \
		echo $${file} ; \
		\
		grep -v '^#' $${file} \
		> a.tmp ; \
		\
		cut -f 1-5 a.tmp \
		| sed -e 's/	/___/g' \
		> b.tmp ; \
		\
		cut -f 6 a.tmp \
		| paste.pl b.tmp - \
		> c.tmp ; \
		\
		join.pl -1 1 -2 1 -o "__NONE__" 2.tmp c.tmp \
		> d.tmp ; \
		\
		grep '^#' $${file} \
		>> comments.tmp ; \
		\
		mv d.tmp 2.tmp ; \
		\
		rm -f a.tmp b.tmp c.tmp d.tmp ; \
		\
	done ;
	\
	head -n 1 comments.tmp \
	> 3.tmp ;
	\
	cat comments.tmp \
	| grep -vf 3.tmp \
	| cat 3.tmp - \
	> 4.tmp ;
	\
	\
	cut -f 1 2.tmp \
	| sed -e 's/___/	/g' \
	> x.tmp ;
	\
	tail -n +2 x.tmp \
	| sed -e 's/	R/	/' \
		-e 's/	F/	/' \
		-e 's/	TEST_/	/' \
	> y.tmp ;
	\
	head -n 1 x.tmp \
	| cat - y.tmp \
	> z.tmp ;
	\
	cut -f 2- 2.tmp \
	| paste.pl z.tmp - \
	> 5.tmp ;
	\
	\
	cat 4.tmp 5.tmp \
	> 6.tmp ;
	\
	\
	mv 6.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp 3.tmp 4.tmp 5.tmp 6.tmp comments.tmp x.tmp y.tmp z.tmp ;
	\

collect_predictions:
	rm -rf $@_temp ;
	\
	for cutoff in $(REDUCED_MODELS_CUTOFFS) ; do \
		echo $${cutoff} ; \
		\
		make cv_test_sample_predictions_reduced_$${cutoff}.tsv ; \
	done ;
	\
	mkdir -p $@_temp ;
	\
	find . -maxdepth 1 -type f -name "cv_test_sample_predictions_reduced_*.tsv" -exec mv {} ./$@_temp/. \; ;
	\
	mv $@_temp $@ ;
	\

# R5:F5_cutoff_500_rf_reduced_model_predictions.RData
# might need to use jq to write out sample label probabilities.
cv_test_sample_predictions_reduced_%.tsv:
	cut -f 1,2 $(DATA_DIR)/cv_folds.tsv \
	> labeled_samples.tmp ;
	\
	rm -f 1.tmp ;
	\
	for repeat_fold in $(shell ls -1 ./models/*_rf_reduced_model_predictions.RData | grep "_cutoff_$*_rf_" | xargs -n 1 basename | cut -d "_" -f 1 ) ; do \
		echo $${repeat_fold} ; \
		\
		Rscript get_aklimate_prediction_probabilities.R \
			reduced \
			test \
			"./models/$${repeat_fold}_cutoff_$*_rf_reduced_model_predictions.RData" \
		; \
		\
		Rscript get_aklimate_prediction_probabilities.R \
			reduced \
			train \
			"./models/$${repeat_fold}_cutoff_$*_rf_reduced_model.RData" \
			"./models/$${repeat_fold}_junkle_final_model.RData" \
		; \
		\
		\
		cat reduced_train_probs.tsv \
		| sed -e 's/	X/	/g' \
		> train1.tmp ; \
		\
		source ~/softwares/venv/3/bin/activate && \
		python ./tmp_aklimate_prediction_matrix_to_json.py \
		< train1.tmp \
		> train2.tmp \
		; \
		\
		paste.pl train2.tmp "TEST_0" \
		> train3.tmp ; \
		\
		source ~/softwares/venv/3/bin/activate && \
		python ./tmp_aklimate_prediction_matrix_to_json.py \
		< reduced_test_probs.tsv \
		> test1.tmp \
		; \
		\
		paste.pl test1.tmp "TEST_1" \
		> test2.tmp ; \
		\
		cat train3.tmp test2.tmp \
		| paste.pl - "$${repeat_fold}" \
		> a.tmp ; \
		\
		join.pl -1 1 -2 1 -o "__NONE__" labeled_samples.tmp a.tmp \
		> b.tmp ; \
		\
		cat b.tmp \
		| awk 'BEGIN{FS="\t";OFS=FS} {if (NR==1) {$$2="Label"; $$3="predicted_label"; $$4="Test"; $$5="Repeat:Fold"} print $$0}' \
		> c.tmp ; \
		\
		cut.pl -f 1,5,4,2,3 c.tmp \
		| sed -e 's/\:/	/' \
		> d.tmp ; \
		\
		echo $(THIS_DIR_NAME) \
		| cut -d "_" -f 3,4 \
		| paste.pl -d "_" - "p" \
		| tr "_" "|" \
		| sed -e 's/^\([a-zA-Z0-9]\+\)|\(.*\)$$/AKLIMATE_\1_REDUCED_$*|\2/' \
		> x.tmp ; \
		\
		cut -f 6 d.tmp \
		| tail -n +2 \
		| cat x.tmp - \
		> y.tmp ; \
		\
		paste.pl d.tmp y.tmp \
		| cut.pl -f 1--3,-1 \
		> z.tmp ; \
		\
		head -n 1 z.tmp \
		> header.tmp ; \
		\
		tail -n +2 z.tmp \
		>> 1.tmp ; \
		\
		rm -f a.tmp b.tmp c.tmp d.tmp test1.tmp test2.tmp train1.tmp train2.tmp train3.tmp x.tmp y.tmp z.tmp ; \
		\
	done ;
	\
	cat header.tmp 1.tmp \
	> 2.tmp ;
	\
	\
	\
	echo $(THIS_DIR_NAME) \
	| cut -d "_" -f 3,4 \
	| paste.pl -d "_" - "p" \
	| tr "_" "|" \
	| sed -e 's/^\([a-zA-Z0-9]\+\)|\(.*\)$$/AKLIMATE_\1_REDUCED_$*|\2/' \
	> x.tmp ;
	\
	head -n 1 2.tmp \
	| cut.pl -f 1--2 \
	| paste.pl - "`cat x.tmp`" \
	> new_header.tmp ;
	\
	tail -n +2 2.tmp \
	| cat new_header.tmp - \
	> 3.tmp ;
	\
	\
	cat tarball.txt \
	> comment_section.tmp ;
	\
	echo "" \
	>> comment_section.tmp ;
	\
	cat x.tmp \
	| paste.pl - "reduced AKLIMATE model with $* features" \
	>> comment_section.tmp ;
	\
	echo "" \
	>> comment_section.tmp ;
	\
	cat comment_section.tmp \
	| sed -e 's/^/\#/' \
	| cat - 3.tmp \
	> 4.tmp ;
	\
	\
	\
	mv 4.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp labeled_sample.tmp ;
	\
	rm -f 3.tmp 4.tmp comment_section.tmp header.tmp labeled_samples.tmp new_header.tmp x.tmp ;
	\

# Rscript get_aklimate_prediction_probabilities.R full test $${file} ;
cv_test_sample_predictions_full.tsv:
	cut -f 1,2 $(DATA_DIR)/cv_folds.tsv \
	> labeled_samples.tmp ;
	\
	rm -f 1.tmp ;
	\
	for repeat_fold in $(shell ls -1 ./models/*_junkle_final_model_stats_preds.RData | xargs -n 1 basename | cut -d "_" -f 1 | sort | uniq | head -n 1) ; do \
		echo $${repeat_fold} ; \
		\
		Rscript get_aklimate_prediction_probabilities.R \
			full \
			test \
			"./models/$${repeat_fold}_junkle_final_model_stats_preds.RData" \
		; \
		\
		Rscript get_aklimate_prediction_probabilities.R \
			full \
			train \
			"./models/$${repeat_fold}_junkle_final_model.RData" \
		; \
		\
		\
		source ~/softwares/venv/3/bin/activate && \
		python ./tmp_aklimate_prediction_matrix_to_json.py \
			-v \
		< full_train_probs.tsv \
		> train1.tmp \
		; \
		\
		paste.pl train1.tmp "TEST_0" \
		> train2.tmp ; \
		\
		\
		transpose.pl full_test_probs.tsv \
		| tail -n +2 \
		> test1.tmp ; \
		\
		cat test1.tmp \
		| cut -f 1 \
		| tr "\." "-" \
		| paste.pl - test1.tmp \
		| cut -f 1,3 \
		| paste.pl - "TEST_1" \
		> test2.tmp ; \
		\
		cat train2.tmp test2.tmp \
		| join.pl -1 1 -2 1 -o "__NONE__" labeled_samples.tmp - \
		> a.tmp ; \
		\
		cat a.tmp \
		| awk 'BEGIN{FS="\t";OFS=FS} {if (NR==1) {$$2="Label"; $$3="predicted_label"; $$4="Test"} print $$0}' \
		> d.tmp ; \
		\
		paste.pl d.tmp "$${repeat_fold}" \
		> e.tmp ; \
		\
		cut.pl -f 1,5,4,2,3 e.tmp \
		| awk 'BEGIN{FS="\t";OFS=FS} {if (NR==1) {$$2="Repeat:Fold"} print $$0 }' \
		> f.tmp ; \
		\
		tail -n +2 f.tmp \
		>> 1.tmp ; \
		\
		head -n 1 f.tmp \
		> header.tmp ; \
	done ; \
	\
	\
	cat header.tmp 1.tmp \
	| sed -e 's/\:/	/' \
	> 2.tmp ;
	\
	echo $(THIS_DIR_NAME) \
	| cut -d "_" -f 3,4 \
	| paste.pl -d "_" - "p" \
	| tr "_" "|" \
	| sed -e 's/^\([a-zA-Z0-9]\+\)|\(.*\)$$/AKLIMATE_\1_FULL|\2/' \
	> x.tmp ;
	\
	head -n 1 2.tmp \
	| cut.pl -f 1--2 \
	| paste.pl - "`cat x.tmp`" \
	> new_header.tmp ;
	\
	tail -n +2 2.tmp \
	| cat new_header.tmp - \
	> 3.tmp ;
	\
	\
	cat tarball.txt \
	> comment_section.tmp ;
	\
	echo "" \
	>> comment_section.tmp ;
	\
	cat x.tmp \
	| paste.pl - "Full, non-reduced AKLIMATE model" \
	>> comment_section.tmp ;
	\
	echo "" \
	>> comment_section.tmp ;
	\
	cat comment_section.tmp \
	| sed -e 's/^/\#/' \
	| cat - 3.tmp \
	> 4.tmp ;
	\
	\
	mv 4.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp 3.tmp 4.tmp a.tmp comment_section.tmp d.tmp e.tmp f.tmp header.tmp labeled_samples.tmp new_header.tmp test1.tmp test2.tmp train1.tmp train2.tmp x.tmp ;
	\
	rm -f full_train_probs.tsv full_test_probs.tsv ; \
	\

get_bacc_stats:bacc.tsv bacc_stats.tsv

bacc_stats.tsv:bacc.tsv
	cat $< \
	| transpose.pl \
	| row_stats.pl -h 1 -k 0 -allstats \
	| transpose.pl \
	> 1.tmp ;
	\
	mv 1.tmp $@ ;
	\
	rm -f 1.tmp ;
	\

bacc.tsv:
	rm -f 1.tmp ;
	\
	for file in $(shell find . -type f -name "*_junkle_final_model_stats_preds.RData") ; do \
		echo $${file} ; \
		\
		echo "$${file}" \
		> a.tmp ; \
		\
		Rscript ./extract_bacc.R $${file} \
		>> a.tmp ; \
		\
		transpose_fast.pl a.tmp \
		>> 1.tmp ; \
		\
	done ;
	\
	cap.pl "task","bacc" 1.tmp \
	> 2.tmp ;
	\
	mv 2.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp a.tmp ;
	\

setup:datatypes.tsv samples.tsv labels.tsv

datatypes.tsv:
	head -n 1 $(COMBINED_MATRIX_FILE) \
	| cut -f 3- \
	| transpose_fast.pl \
	> 1.tmp ;
	\
	cut.pl -d ":" -f 1-3 1.tmp \
	> 2.tmp ;
	\
	uniqc.sh 2.tmp \
	> 3.tmp ;
	\
	mv 3.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp 3.tmp ;
	\

samples.tsv:
	cut -f 1 $(COMBINED_MATRIX_FILE) \
	| tail -n +2 \
	> 1.tmp ;
	\
	mv 1.tmp $@ ;
	\
	rm -f 1.tmp ;
	\

labels.tsv:
	cut -f 1,2 $(COMBINED_MATRIX_FILE) \
	| tail -n +2 \
	> 1.tmp ;
	\
	mv 1.tmp $@ ;
	\
	rm -f 1.tmp ;
	\

all: $(TARGETS)

summary: $(SUMMARY_TARGETS)

clean_all: clean_targets clean_tmp clean_summary

clean_summary:
	rm -f $(SUMMARY_TARGETS) ;

clean_targets:
	rm -f $(TARGETS) ;

clean_tmp:
	rm -f $(wildcard *.tmp) ;
