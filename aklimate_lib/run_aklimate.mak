# chrisw

THIS_DIR=$(shell pwd)

THIS_DIR_NAME=$(shell pwd | cut.pl -d "/" -f -1)

DATA_DIR=data
LIB_DIR=lib

COMBINED_MATRIX_FILE=$(DATA_DIR)/combined_matrix.tsv
CV_SETS_FILE=$(DATA_DIR)/cv_folds.tsv

TARGETS= \
	cv_test_sample_predictions_full.tsv \
	datatypes.tsv \
	samples.tsv \
	labels.tsv \
	\

SUMMARY_TARGETS= \
	predictions.tar.gz \
	features.tsv.gz \
	bacc.tsv \
	bacc_stats.tsv \
	cv_test_sample_predictions_full.tsv \
	sickle_plot.png \
	datatype_stacked_bar_plots.png \
	feature_importance_stats.tsv \
	bal_acc_subtype_50_cutoff.png \
	feature_set_weights.tsv \
	\

REDUCED_MODELS_CUTOFFS= \
	$(shell find ./models/ -type f -name "*rf_reduced_model_predictions.RData" | cut.pl -d "cutoff_" -f 2 | cut -d "_" -f 1| sort.pl | uniq)

test:

# AKLIMATE pipeline assumes features reverse-ordered by importance and truncates the list at each cutoff.
# The short list of features is then passed to ranger to build RF models.
# The format for this file is described at https://www.synapse.org/#!Synapse:syn8011998/wiki/600601
features.tsv.gz:
	head -n 1 $(COMBINED_MATRIX_FILE) \
	| cut -f 1 \
	> cohort.tmp ;
	\
	rm -f 1.tmp ;
	\
	for fold in $(shell ls -1 ./models | grep "_aklimate_multiclass_feature_importance.tab" | cut.pl -d "_aklimate_multiclass_feature_importance.tab" -f 1 ) ; do \
		echo $${fold} ; \
		\
		tail -n +2 "./models/$${fold}_aklimate_multiclass_feature_importance.tab" \
		| sort.pl -r -h 0 -k 2 \
		> a.tmp ; \
		\
		for cutoff in $(REDUCED_MODELS_CUTOFFS) ; do \
			echo $${fold}_$${cutoff} ; \
			\
			head -n $${cutoff} a.tmp \
			| cut -f 1 \
			| transpose_fast.pl \
			> b.tmp ; \
			\
			cat cohort.tmp b.tmp \
			> c.tmp ; \
			\
			source ~/softwares/venv/3/bin/activate && \
				python ./tsv_to_json_lists.py \
			< c.tmp \
			> d.tmp \
			; \
			\
			\
			echo "$(THIS_DIR_NAME)" \
			| cut.pl -d "_" -f 1,3,4 \
			> aa.tmp ; \
			\
			echo "$${fold}_reduced_$${cutoff}_features" \
			| cat aa.tmp - \
			| transpose.pl \
			| sed -e 's/	/_/g' \
			> bb.tmp ; \
			\
			cat bb.tmp d.tmp \
			> e.tmp ; \
			\
			transpose_fast.pl e.tmp \
			> f.tmp ; \
			\
			cat f.tmp \
			>> 1.tmp ; \
			\
			rm -f aa.tmp bb.tmp b.tmp c.tmp e.tmp f.tmp ; \
			\
		done ; \
		\
		rm -f a.tmp ; \
	\
	done ; \
	\
	cat 1.tmp \
	| cap.pl "Feature_Set_ID","TCGA_Projects","Features" \
	> 2.tmp ;
	\
	gzip -c 2.tmp \
	> 3.tmp ;
	\
	mv 3.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp 3.tmp a.tmp b.tmp c.tmp cohort.tmp d.tmp e.tmp f.tmp aa.tmp bb.tmp ;
	\

# AKLIMATE pipeline assumes features reverse-ordered by importance and truncates the list at each cutoff.
# The short list of features is then passed to ranger to build RF models.
collected_aklimate_features_old.tsv:
	rm -f 1.tmp ;
	\
	for cutoff in $(REDUCED_MODELS_CUTOFFS) ; do \
		echo $${cutoff} ; \
		\
		rm -f a.tmp ; \
		\
		for fold in $(shell ls -1 ./models | grep "_aklimate_multiclass_feature_importance.tab" | cut.pl -d "_aklimate_multiclass_feature_importance.tab" -f 1 ) ; do \
			echo $${cutoff}_$${fold} ; \
			\
			tail -n +2 "./models/$${fold}_aklimate_multiclass_feature_importance.tab" \
			| sort.pl -r -h 0 -k 2 \
			| head -n $${cutoff} \
			>> a.tmp ; \
			\
		done ; \
		\
		cat a.tmp \
		| expand.pl \
		| fill.pl "0" \
		| row_stats.pl -h 0 -k 0 -allstats \
		| sort.pl -h 1 -k -1 -r \
		> $${cutoff}_b.tmp ; \
		\
		head -n 1 $(COMBINED_MATRIX_FILE) \
		| cut -f 1 \
		> $${cutoff}_c.tmp ; \
		\
		cut -f 1 $${cutoff}_b.tmp \
		| transpose.pl \
		| cut -f 2- \
		>> $${cutoff}_c.tmp ; \
		\
		source ~/softwares/venv/3/bin/activate && \
			python ./tsv_to_json_lists.py \
		< $${cutoff}_c.tmp \
		> $${cutoff}_d.tmp \
		; \
		\
		echo "$(THIS_DIR_NAME)_reduced_features_$${cutoff}" \
		| cat - $${cutoff}_d.tmp \
		> $${cutoff}_e.tmp ; \
		\
		transpose.pl $${cutoff}_e.tmp \
		>> 1.tmp ; \
		\
	done ;
	\
	cat 1.tmp \
	| cap.pl "Feature_Set_ID","TCGA_Projects","Features" \
	> 2.tmp ;
	\
	mv 2.tmp $@ ;
	\
	rm -f a.tmp 1.tmp 2.tmp *_a.tmp *_b.tmp *_c.tmp *_d.tmp *_e.tmp ;
	\

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
			echo "defaulting to handling feature set weights from binary model" ; \
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
	| sed -e 's/_MUTA//' \
		-e 's/_CNVR//' \
		-e 's/_METH//' \
		-e 's/_GEXP//' \
		-e 's/_[0-9]\+$$//' \
	| cut.pl -d "_" -f 1--2 \
	> 2.tmp ;
	\
	paste.pl 2.tmp 1.tmp \
	> 3.tmp ;
	\
	join.pl -r -1 1 -2 1 -o "__NONE__" 3.tmp ./p_store_files/pathway_name_mapping.tsv \
	> 4.tmp ;
	\
	cut -f 1,2,5 4.tmp \
	| sed -e 's/	/____/' \
	| expand.pl \
	> 5.tmp ;
	\
	cat 5.tmp \
	| row_stats.pl -h 0 -k 0 -allstats \
	| sort.pl -h 1 -k -1 -r \
	> 6.tmp ;
	\
	mv 6.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp 3.tmp 4.tmp 5.tmp 6.tmp ;
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
	| sort.pl -h 1 -k -1 -r \
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

# collect sample classification predictions
# file format is at https://www.synapse.org/#!Synapse:syn8011998/wiki/600601
predictions.tar.gz:
	\
	rm -rf predictions ;
	\
	mkdir -p predictions ;
	\
	cut -f 1,2 $(DATA_DIR)/cv_folds.tsv \
	| tail -n +2 \
	> labeled_samples.tmp ;
	\
	for fold in $(shell ls -1 ./models | grep "_aklimate_multiclass_feature_importance.tab" | cut.pl -d "_aklimate_multiclass_feature_importance.tab" -f 1 ) ; do \
		echo $${fold} ; \
		\
		cut -f 1 labeled_samples.tmp \
		| cap.pl "Sample_ID" \
		> 1.tmp ; \
		\
		rm -f header.tmp ; \
		\
		for cutoff in $(REDUCED_MODELS_CUTOFFS) ; do \
			echo "$${fold}_$${cutoff}" ; \
			\
			Rscript get_aklimate_prediction_probabilities.R \
				reduced \
				test \
				"./models/$${fold}_cutoff_$${cutoff}_rf_reduced_model_predictions.RData" \
			; \
			\
			Rscript get_aklimate_prediction_probabilities.R \
				reduced \
				train \
				"./models/$${fold}_cutoff_$${cutoff}_rf_reduced_model.RData" \
				"./models/$${fold}_junkle_final_model.RData" \
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
			paste.pl train2.tmp "0" \
			> train3.tmp ; \
			\
			source ~/softwares/venv/3/bin/activate && \
			python ./tmp_aklimate_prediction_matrix_to_json.py \
			< reduced_test_probs.tsv \
			> test1.tmp \
			; \
			\
			paste.pl test1.tmp "1" \
			> test2.tmp ; \
			\
			\
			cat train3.tmp test2.tmp \
			| paste.pl - "$${fold}	$${cutoff}" \
			> a.tmp ; \
			\
			join.pl -1 1 -2 1 -o "__NONE__" labeled_samples.tmp a.tmp \
			> b.tmp ; \
			\
			cat b.tmp \
			| cap.pl "Sample_ID","Label","prediction","Test","Repeat:Fold","cutoff" \
			> c.tmp ; \
			\
			cut.pl -f 1,-2,4,2,3 c.tmp \
			> d.tmp ; \
			\
			cat d.tmp \
			| sed -e 's/\:/	/' \
				-e 's/	R\([0-9]\+\)	/	\1	/' \
				-e 's/	F\([0-9]\+\)	/	\1	/' \
			> e.tmp ; \
			\
			\
			echo "$(THIS_DIR)" \
			| cut.pl -d "/" -f -1 \
			| sed -e 's/^AKLIMATE_TEMPLATE_/AKLIMATE_/' \
			| cut -d "_" -f 1,2,3 \
			> f.tmp ; \
			\
			cut -d "_" -f 1,2 f.tmp \
			> aklimate_cohort.tmp ; \
			\
			cut -d "_" -f 3 f.tmp \
			> model_date.tmp ; \
			\
			cat aklimate_cohort.tmp \
			| cat - model_date.tmp \
			> prediction_colname.tmp ; \
			\
			echo "$${fold}" \
			>> prediction_colname.tmp ; \
			\
			echo "reduced" \
			>> prediction_colname.tmp ; \
			\
			echo "$${cutoff}" \
			>> prediction_colname.tmp ; \
			\
			transpose.pl prediction_colname.tmp \
			| sed -e 's/	/_/g' \
			> prediction_colname_2.tmp ; \
			\
			cat prediction_colname_2.tmp prediction_colname_2.tmp \
			| sed -e '1s/$$/_features/' \
			| tac \
			> prediction_colname_3.tmp ; \
			\
			cat model_date.tmp \
			>> prediction_colname_3.tmp ; \
			\
			echo "p" \
			>> prediction_colname_3.tmp ; \
			\
			transpose.pl prediction_colname_3.tmp \
			| sed -e 's/	/\|/g' \
			> prediction_colname_4.tmp ; \
			\
			\
			head -n 1 e.tmp \
			| transpose.pl \
			| head -n 5 \
			| cat - prediction_colname_4.tmp \
			| transpose.pl \
			> x.tmp ; \
			\
			tail -n +2 e.tmp \
			| cat x.tmp - \
			> y.tmp ; \
			\
			cut.pl -f 1,-1 y.tmp \
			> z.tmp ; \
			\
			join.pl -1 1 -2 1 -o "NA" 1.tmp z.tmp \
			| fill.pl "NA" \
			> zz.tmp ; \
			\
			mv zz.tmp 1.tmp ; \
			\
			\
			echo "#`cat prediction_colname_2.tmp`	AKLIMATE model trained on $${cutoff} features and the $${fold} CV_fold" \
			>> header.tmp ; \
			\
			\
			rm -f reduced_train_probs.tsv reduced_test_probs.tsv ; \
			\
			done ; \
		\
		\
		cut.pl -f 1--2 y.tmp \
		| join.pl -1 1 -2 1 - 1.tmp \
		> 2.tmp ; \
		\
		cat tarball.txt \
		| sed -e 's/_organized//i' \
			-e 's/^/#/' \
		| cat - header.tmp 2.tmp \
		> 3.tmp ; \
		\
		mv 3.tmp predictions/predictions_`cat aklimate_cohort.tmp`_$${fold}.tsv ; \
		\
		rm -f 1.tmp 2.tmp 3.tmp aklimate_cohort.tmp a.tmp b.tmp c.tmp d.tmp e.tmp f.tmp header.tmp model_date.tmp prediction_colname_2.tmp prediction_colname_3.tmp prediction_colname_4.tmp prediction_colname.tmp test1.tmp test2.tmp train1.tmp train2.tmp train3.tmp x.tmp y.tmp z.tmp ; \
		\
		done ; \
	\
	tar -zcvf predictions.tar.gz predictions/ ;
	\
	rm -f labeled_samples.tmp ;
	\
	rm -rf predictions ;
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
	for file in $(shell find ./models -type f -name "*_junkle_final_model_stats_preds.RData") ; do \
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
	rm -rf $(SUMMARY_TARGETS) ;

clean_targets:
	rm -f $(TARGETS) ;

clean_tmp:
	rm -f $(wildcard *.tmp) ;
