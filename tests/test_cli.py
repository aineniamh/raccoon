import pytest

from raccoon.command import build_parser


def test_subcommands_present():
    parser = build_parser()
    # access to subparser choices (stable across argparse implementations)
    subparsers = parser._subparsers._group_actions[0].choices
    assert 'alignment' in subparsers
    assert 'phylo' in subparsers
    assert 'combine' in subparsers


def test_alignment_sets_func_callable():
    parser = build_parser()
    ns = parser.parse_args(['alignment', 'input.fasta'])
    assert callable(ns.func)


def test_phylo_sets_func_callable():
    parser = build_parser()
    ns = parser.parse_args(['phylo', '--phylogeny', 'tree'])
    assert callable(ns.func)


def test_combine_sets_func_callable():
    parser = build_parser()
    ns = parser.parse_args(['combine', 'a.fasta'])
    assert callable(ns.func)
