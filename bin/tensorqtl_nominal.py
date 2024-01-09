#!/home/workspace/jogrady/env/bin/python3
import pandas as pd
import numpy as np
import tensorqtl
from tensorqtl import genotypeio, cis, post
import qtl.plot
import sys

pr = genotypeio.PlinkReader(sys.argv[1])
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

id = sys.argv[2]

covariates_df = pd.read_csv(sys.argv[3], sep='\t', index_col=0).T
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(sys.argv[4])
group_s = pd.read_csv(sys.argv[5],sep='\t', header=None, index_col=0).squeeze("columns")
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, id, covariates_df=covariates_df, group_s=group_s, output_dir=sys.argv[6])