from Bio.Data.IUPACData import ambiguous_dna_values as triv
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import product
import primer3
import pandas as pd
import numpy as np


primers = list(SeqIO.parse('../data/primers_test.fasta', 'fasta'))

def unpack(seq):
    pools = [triv[i] for i in seq]
    return [''.join(pos) for pos in product(*pools)]


variants = {r.id: unpack(r.seq) for r in primers}


def thermo(seq):
    hairpin   = primer3.calc_hairpin(seq)
    homodimer = primer3.calc_homodimer(seq)
    duplex    = primer3.calc_heterodimer(seq, str(Seq(seq).reverse_complement()))

    return {
        'tm':           primer3.calc_tm(seq),
        'dg':           duplex.dg,
        'hairpin_tm':   hairpin.tm if hairpin.structure_found else np.nan,
        'homodimer_tm': homodimer.tm if homodimer.structure_found else np.nan,
    }


res = []

for id, variant in variants.items():
    for seq in variant:
        var_stats = {'primer': id, 'variant': seq}
        var_stats.update(thermo(seq))
        res.append(var_stats)

df = pd.DataFrame(res)


params = ['tm', 'dg', 'hairpin_tm', 'homodimer_tm']
aggs   = ['mean', 'median', lambda x: x.quantile(0.05), lambda x: x.quantile(0.95)]

per_primer = df.groupby('primer')[params].agg(aggs)
overall    = df[params].agg(aggs)

names = {'<lambda_0>': 'p5', '<lambda_1>': 'p95'}
per_primer = per_primer.rename(columns=names, level=1)
overall.index = ['mean', 'median', 'p5', 'p95']


with pd.ExcelWriter('results/task1.xlsx') as writer:
    per_primer.to_excel(writer, sheet_name='stats_per_primer')
    overall.to_excel(writer, sheet_name='stats_overall')