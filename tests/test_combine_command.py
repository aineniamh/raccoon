import textwrap
from pathlib import Path

from raccoon.commands import combine


class MockArgs:
    def __init__(self, **kwargs):
        self.inputs = kwargs.get("inputs", [])
        self.output = kwargs.get("output", None)
        self.metadata = kwargs.get("metadata", None)
        self.metadata_delimiter = kwargs.get("metadata_delimiter", ",")
        self.metadata_id_field = kwargs.get("metadata_id_field", "id")
        self.metadata_location_field = kwargs.get("metadata_location_field", "location")
        self.metadata_date_field = kwargs.get("metadata_date_field", "date")
        self.header_separator = kwargs.get("header_separator", "|")


def _write_fasta(path, entries):
    lines = []
    for header, seq in entries:
        lines.append(f">{header}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n")


def test_combine_uppercase_and_unwrapped(tmp_path):
    f1 = tmp_path / "a.fasta"
    f2 = tmp_path / "b.fasta"

    _write_fasta(f1, [("seq1", "acgtacgt"), ("seq2", "aaaa")])
    _write_fasta(f2, [("seq3", "tttt")])

    out = tmp_path / "combined.fasta"
    args = MockArgs(inputs=[str(f1), str(f2)], output=str(out))

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1"
    assert lines[1] == "ACGTACGT"
    assert lines[2] == ">seq2"
    assert lines[3] == "AAAA"
    assert lines[4] == ">seq3"
    assert lines[5] == "TTTT"


def test_combine_harmonises_headers_with_metadata(tmp_path):
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            id,location,date
            seq1,UK,2024-01-01
            seq2,US,2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        inputs=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1|UK|2024-01-01"
    assert lines[1] == "ACGT"
    assert lines[2] == ">seq2|US|2024-02-02"
    assert lines[3] == "GGGG"


def test_combine_examples_same_headers(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    examples = repo_root / "examples"
    out = tmp_path / "combined_same_headers.fasta"

    args = MockArgs(
        inputs=[
            str(examples / "combine_a.fasta"),
            str(examples / "combine_b.fasta"),
        ],
        output=str(out),
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">sample1"
    assert lines[1] == "ACGTTGCAACGT"
    assert lines[2] == ">sample2"
    assert lines[3] == "ACGTNNNNACGT"
    assert lines[4] == ">sample3"
    assert lines[5] == "TTGGCCAAGGAA"
    assert lines[6] == ">sample4"
    assert lines[7] == "NNNNACGTACGT"


def test_combine_examples_harmonised_headers(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    examples = repo_root / "examples"
    out = tmp_path / "combined_harmonised_headers.fasta"

    args = MockArgs(
        inputs=[
            str(examples / "combine_a.fasta"),
            str(examples / "combine_b.fasta"),
        ],
        output=str(out),
        metadata=[str(examples / "combine_metadata.csv")],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">sample1|UK|2024-01-01"
    assert lines[1] == "ACGTTGCAACGT"
    assert lines[2] == ">sample2|US|2024-02-02"
    assert lines[3] == "ACGTNNNNACGT"
    assert lines[4] == ">sample3|FR|2024-03-03"
    assert lines[5] == "TTGGCCAAGGAA"
    assert lines[6] == ">sample4|ZA|2024-04-04"
    assert lines[7] == "NNNNACGTACGT"


def test_combine_multiple_metadata_files(tmp_path):
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    meta1 = tmp_path / "meta1.csv"
    meta1.write_text(
        textwrap.dedent(
            """\
            id,location,date
            seq1,UK,2024-01-01
            """
        )
    )
    meta2 = tmp_path / "meta2.csv"
    meta2.write_text(
        textwrap.dedent(
            """\
            id,location,date
            seq2,US,2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        inputs=[str(fasta_path)],
        output=str(out),
        metadata=[str(meta1), str(meta2)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1|UK|2024-01-01"
    assert lines[1] == "ACGT"
    assert lines[2] == ">seq2|US|2024-02-02"
    assert lines[3] == "GGGG"
