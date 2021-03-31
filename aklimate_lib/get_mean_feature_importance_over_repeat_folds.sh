\
FILES=`find -L ./models -type f -name "*_aklimate_multiclass_feature_importance.tab"`
\
NUM_FILES=`find -L ./models -type f -name "*_aklimate_multiclass_feature_importance.tab" | wc -l`
\
\rm -f a.tmp
\
for file in $FILES ; do \
	echo ${file} ; \
	\
	cat ${file} \
	| tail -n +2 \
	>> a.tmp ; \
	\
done ;
\
cat a.tmp \
| expand.pl \
| fill.pl "0" \
> b.tmp ;
\
cat b.tmp \
| row_stats.pl -h 0 -k 0 -allstats \
> c.tmp ;
\
cat c.tmp \
| awk -v num_files=$NUM_FILES 'BEGIN{FS="\t";OFS=FS} {if (NR==1) {$11="mean_importance_score"} else {$11=$10/num_files} print $0}' \
> d.tmp ;
\
cat d.tmp \
| sort.pl -h 1 -k 11 -r \
> e.tmp ;
\
cat e.tmp \
| cut -f 1,11 \
> f.tmp ;
\
mv f.tmp mean_feature_importance_over_repeat_folds.tsv
\
rm -f a.tmp b.tmp c.tmp d.tmp e.tmp f.tmp ;
\
