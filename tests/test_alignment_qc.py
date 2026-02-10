import textwrap
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

from raccoon.utils import alignment_functions as af


from Bio.Align import MultipleSeqAlignment


def test_find_high_n_sequences(tmp_path):
    r1 = SeqRecord(Seq('ATGATGATGATG'), id='seq1')
    r2 = SeqRecord(Seq('ATGNNNNNNATG'), id='seq2')
    aln = MultipleSeqAlignment([r1, r2])
    path = tmp_path / 'test.fasta'
    AlignIO.write(aln, str(path), 'fasta')

    aln2 = AlignIO.read(str(path), 'fasta')
    flagged = af.find_high_N_sequences(aln2, threshold=0.3)
    assert len(flagged) == 1
    assert flagged[0][0] == 'seq2'


def test_find_clustered_snps(tmp_path):
    # create a small alignment where seq2 has 4 unique snps within a 10bp window
    ref_seq = 'A' * 24
    ref = SeqRecord(Seq(ref_seq), id='ref')
    s1 = SeqRecord(Seq(ref_seq), id='s1')

    # build s2 by mutating a few positions but keeping length the same
    s2_list = list(ref_seq)
    for pos, nt in [(3, 'T'), (5, 'G'), (7, 'A'), (9, 'T')]:
        s2_list[pos] = nt
    s2 = SeqRecord(Seq(''.join(s2_list)), id='s2')

    aln = MultipleSeqAlignment([ref, s1, s2])
    path = tmp_path / 'cluster.fasta'
    AlignIO.write(aln, str(path), 'fasta')
    aln2 = AlignIO.read(str(path), 'fasta')

    sites = af.find_clustered_snps(aln2, window=10, min_snps=3)
    # expect some sites flagged (the mutated positions)
    assert len(sites) > 0
    # ensure 's2' appears in at least one present_in
    present = set()
    for v in sites.values():
        present.update(v['present_in'])
    assert 's2' in present
