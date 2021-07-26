# Amy Campbell
# run from local machine using conda environment qiime2-2021.4
# 2021-07-26 
# Making tree from 97% ID clustered sequences

qiime alignment mafft \
  --i-sequences /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/rep-seqs-cr-97.qza \
  --o-alignment /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/aligned-rep-seqs-cr-97.qza

qiime phylogeny fasttree \
  --i-alignment /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/aligned-rep-seqs-cr-97.qza \
  --o-tree /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/cr-97-fasttree.qza

qiime tools export \
  --input-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/cr-97-fasttree.qza \
  --output-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/cr-97-fasttree.newick
