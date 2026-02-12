import textwrap
import os

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from raccoon.utils import alignment_functions as af


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


def test_analyze_alignment_unique_snps_and_flags():
    # alignment where s2 has three unique SNPs, one near N and one near a gap
    ref_seq = 'AAAAAAAAAAAA'
    ref = SeqRecord(Seq(ref_seq), id='ref')
    s1 = SeqRecord(Seq(ref_seq), id='s1')
    s2 = SeqRecord(Seq('AATAGNA-CAAA'), id='s2')

    aln = MultipleSeqAlignment([ref, s1, s2])
    unique_mutations, snps_near_n, snps_near_gap, clustered_snps = af.analyze_alignment(
        aln, n_window=2, gap_window=1, snp_window=10, snp_count=3
    )

    assert unique_mutations['s2'] == {2, 4, 8}
    assert snps_near_n['s2'] == {4}
    assert snps_near_gap['s2'] == {8}
    assert clustered_snps['s2'] == {2, 4, 8}


def test_sliding_window():
    assert list(af.sliding_window([1, 2, 3, 4], 2)) == [[1, 2], [2, 3], [3, 4]]
    assert list(af.sliding_window([1, 2], 3)) == [[1, 2]]


def test_find_frame_breaking_indels(tmp_path):
    # reference sequence with a CDS covering the full length
    ref_record = SeqRecord(Seq('ATGAAATTT'), id='ref', name='ref')
    ref_record.annotations['molecule_type'] = 'DNA'
    ref_record.features.append(SeqFeature(FeatureLocation(0, 9), type='CDS'))

    genbank_path = tmp_path / 'ref.gb'
    with open(genbank_path, 'w') as handle:
        from Bio import SeqIO
        SeqIO.write(ref_record, handle, 'genbank')

    # alignment where seq2 has a single gap in the CDS region -> frame break
    ref_aln = SeqRecord(Seq('ATGAAATTT'), id='ref')
    seq2 = SeqRecord(Seq('ATG-AATTT'), id='seq2')
    aln = MultipleSeqAlignment([ref_aln, seq2])

    sites = af.find_frame_breaking_indels(aln, str(genbank_path), reference_id='ref')
    assert sites, "Expected frame-breaking sites to be reported"
    # ensure at least one site includes seq2 in present_in
    assert any('seq2' in v['present_in'] for v in sites.values())
