# chrisw

THIS_DIR=$(shell pwd)

THIS_DIR_NAME=$(shell pwd | cut.pl -d "/" -f -1)

DATA_DIR=data
LIB_DIR=lib

COMBINED_MATRIX_FILE=$(DATA_DIR)/combined_matrix.tsv
CV_SETS_FILE=$(DATA_DIR)/cv_folds.tsv

TARGETS= \
	datatypes.tsv \
	samples.tsv \
	labels.tsv \
	\

test:cv_test_sample_predictions_full.tsv

cv_test_sample_predictions_full.tsv:
	cut -f 1,2 $(DATA_DIR)/cv_folds.tsv \
	> labeled_samples.tmp ;
	\
	rm -f 1.tmp ;
	\
	for file in $(shell find ./models -type f -name "*_junkle_final_model_stats_preds.RData" ) ; do \
		echo $${file} ; \
		\
		Rscript get_sample_predictions.R $${file} ; \
		\
		\
		transpose.pl preds.tmp \
		| tail -n +2 \
		> a.tmp ; \
		\
		cat a.tmp \
		| cut -f 1 \
		| tr "\." "-" \
		| paste.pl - a.tmp \
		| cut -f 1,3 \
		> b.tmp ; \
		\
		join.pl -1 1 -2 1 -o "TEST_0" labeled_samples.tmp b.tmp \
		> c.tmp ; \
		\
		cat c.tmp \
		| awk 'BEGIN{FS="\t";OFS=FS} {if (NR==1) {$$2="Label"; $$3="predicted_label"; $$4="Test"} else {if ($$3=="TEST_0") {$$3=$$2; $$4="0"} else {$$4="1"} } print $$0}' \
		> d.tmp ; \
		\
		paste.pl d.tmp "from_$${file}" \
		| sed -e 's/\(R[0-9]\+:F[0-9]\+\)_.*$$/	\1/' \
		> e.tmp ; \
		\
		cut.pl -f 1,6,4,2,3 e.tmp \
		| awk 'BEGIN{FS="\t";OFS=FS} {if (NR==1) {$$2="Repeat:Fold"} print $$0 }' \
		> f.tmp ; \
		\
		tail -n +2 f.tmp \
		>> 1.tmp ; \
		\
		head -n 1 f.tmp \
		> header.tmp ; \
		\
	done ;
	\
	cat header.tmp 1.tmp \
	| sed -e 's/\:/	/' \
	> 2.tmp ;
	\
	echo $(THIS_DIR_NAME) \
	| cut -d "_" -f 3,4 \
	| paste.pl -d "_" - "c" \
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
	mv 3.tmp $@ ;
	\
	rm -f labeled_samples.tmp preds.tmp a.tmp b.tmp c.tmp d.tmp e.tmp f.tmp 1.tmp header.tmp 2.tmp x.tmp new_header.tmp 3.tmp ;
	\
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

clean_all: clean_targets clean_tmp

clean_targets:
	rm -f $(TARGETS) ;

clean_tmp:
	rm -f $(wildcard *.tmp) ;

