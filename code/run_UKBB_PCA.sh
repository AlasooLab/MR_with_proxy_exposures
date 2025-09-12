#Convert bgen to PLINK
module load any/plink/2-200923
plink2 --bgen PCA_input_original_ukbb.bgen ref-first --sample PCA_input_original_ukbb.sample --make-bed --out UKBB_full_pruned
plink2 --bgen ../ukbb/genotypes/ld_pruned/PCA_input_300K.bgen ref-first --sample ../ukbb/genotypes/ld_pruned/PCA_input_300K.bgen --make-bed --out UKBB_300k_pruned

PCA_input_300K.bgen
#Run FlashPCA
./flashpca_x86-64 --bfile UKBB_full_pruned --outpc UKBB_full_PCA -d20 -n8