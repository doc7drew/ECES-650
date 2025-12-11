#open Project_1 directory
cd /Project_1

#make visualization file
qiime demux summarize \
	--i-data demux.qza \
	--o-visualization demux.qzv

#Start Dada2
#trims reads to 403 bp
#denoises, puts into ASV sequences and ASV abundance table
qiime dada2 denoise-single \
	--p-trunc-len 403 \
	--i-demultiplexed-seqs demux.qza \
	--o-representative-sequences rep-seqs-dada2.qza \
	--o-table table-dada2.qza \
	--o-denoising-stats stats-dada2.qza

qiime metadata tabulate \
	--m-input-file stats-dada2.qza \
	--o-visualization stats-dada2.qzv

#copy table-dada2.qza to table.qza
mv rep-seqs-dada2.qza rep-seqs.qza
mv table-dada2.qza table.qza

#Feature Table and FeatureData summaries
qiime feature-table summarize \
        --i-table table.qza \
        --o-visualization table.qzv \
        --m-sample-metadata-file database.tsv

qiime feature-table tabulate-seqs \
	--i-data rep-seqs.qza \
	--o-visualization rep-seqs.qzv

#maybe try this next?
#mkdir tmp first in home system
#export MAFFT_TMPDIR=/home/ktf47/Project_1/tmp
#generate phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
	--verbose \
	--p-parttree \
	--i-sequences rep-seqs.qza \
	--o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza

#Explore some results
#computes UniFrac distances, PCoA plots, Shannon diversity, evenness
qiime diversity core-metrics-phylogenetic \
	--i-phylogeny rooted-tree.qza \
	--i-table table.qza \
	--p-sampling-depth 23932 \
	--m-metadata-file database.tsv \
	--output-dir core-metrics-results

qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
	--m-metadata-file database.tsv \
	--o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/evenness_vector.qza \
	--m-metadata-file database.tsv \
	--o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity beta-group-significance \
	--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	--m-metadata-file database.tsv \
	--m-metadata-column Group \
	--o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv \
	--p-pairwise

qiime diversity alpha-rarefaction \
	--i-table table.qza \
	--i-phylogeny rooted-tree.qza \
	--p-max-depth 4000 \
	--m-metadata-file database.tsv \
	--o-visualization alpha-rarefaction.qzv

#conversion for BIOM for MicrobiomeAnalyst
qiime tools export \
	--input-path table.qza \
	--output-path ./

#BIOM to TSV
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv 

#Microbiome Analyst Steps:
#remove CP_Mauge
#make nomauge_table.qza
qiime feature-table filter-samples \
	--i-table table.qza \
	--m-metadata-file database.tsv \
	--p-where "[Sample_name]='CP_Mauge'" \
	--p-exclude-ids \
	--o-filtered-table nomauge_table.qza

#create feature-table.biom folder
qiime tools export \
	--input-path nomauge_table.qza \
	--output-path biomtable

#taxonomic analysis
qiime feature-classifier classify-sklearn \
	--i-classifier gg-13-8-99-515-806-nb-classifier.qza \
	--i-reads rep-seqs.qza \
	--o-classification taxonomy.qza

qiime metadata tabulate \
	--m-input-file taxonomy.qza \
	--o-visualization taxonomy.qzv

#export taxonomy to folder
qiime tools export \
	--input-path taxonomy.qza \
	--output-path biomtable

#combine ASV counts and taxonomy
#make Microbiome Analyst file
biom add-metadata \
	-i biomtable/feature-table.biom \
	-o biomtable/table-with-taxonomy.biom \
	--observation-metadata-fp biomtable/taxonomy.tsv \
	--sc-separated taxonomy

#For MicrobiomeAnalyst Upload, need to create metadata file without Mauge row and without 
#q2:types so it doesn't give errors
#save metadata to CSV
