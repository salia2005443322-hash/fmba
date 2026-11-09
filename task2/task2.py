#!/usr/bin/env python3

import pandas as pd
from Bio import Entrez
import taxopy

blast6 = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
          'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

blast_res = pd.read_csv('../data/resp1-1201107674_S85.merged.tsv', sep='\t', header=None, names=blast6)
accessions = blast_res['sseqid'].unique()

handle = Entrez.esummary(db='nuccore', id=','.join(accessions))
acc_to_taxid = {r['AccessionVersion']: int(r['TaxId']) for r in Entrez.read(handle)}
handle.close()

taxdb = taxopy.TaxDb(taxdb_dir='../data/taxdump', keep_files=True)

taxes = []
for taxid in set(acc_to_taxid.values()):
    taxon = taxopy.Taxon(taxid, taxdb)

    if 'Viruses' not in taxon.name_lineage:
        continue

    row = {'taxid': taxid}
    row.update(taxon.rank_name_dictionary)
    taxes.append(row)

taxonomy = pd.DataFrame(taxes)

sp_ranks = ['superkingdom', 'realm', 'kingdom', 'phylum', 'class',
         'order', 'family', 'genus', 'species']
sp_taxonomy = (taxonomy.reindex(columns=sp_ranks).dropna(axis=1, how='all').drop_duplicates(subset='species').reset_index(drop=True))


blast_res['taxid'] = blast_res['sseqid'].map(acc_to_taxid)
virus_blast = blast_res.merge(taxonomy, on='taxid')

virus_blast = (virus_blast.sort_values(['bitscore'], ascending=False).drop_duplicates(['qseqid', 'species']))

stats = virus_blast.groupby('species').agg(
    n_reads=('qseqid', 'nunique'),
    median_identity=('pident', 'median'),
    median_evalue=('evalue', 'median'),
    median_length=('length', 'median')).reset_index()

stats['share_all'] = stats['n_reads'] / blast_res['qseqid'].nunique()
stats['share_virus'] = stats['n_reads'] / virus_blast['qseqid'].nunique()
stats = stats.sort_values('n_reads', ascending=False).reset_index(drop=True)

with pd.ExcelWriter('results/task2.xlsx', engine='openpyxl') as writer:
    sp_taxonomy.to_excel(writer, sheet_name='taxonomy_species', index=False)
    taxonomy.to_excel(writer, sheet_name='taxonomy_subspecies', index=False)
    stats.to_excel(writer, sheet_name='stats', index=False)