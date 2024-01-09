import pandas as pd
import qtl.annotation
import sys

annot = qtl.annotation.Annotation(sys.argv[1])
exon_df = pd.DataFrame([[g.chr, e.start_pos, e.end_pos, g.strand, g.id, g.name]
                        for g in annot.genes for e in g.transcripts[0].exons],
                       columns=['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name'])
exon_df.to_csv(sys.argv[2], sep='\t', index=False)