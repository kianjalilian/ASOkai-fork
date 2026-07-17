#!/usr/bin/env python
"""Tests for antisense oligonucleotide domain classes."""

from ASOkai.Antisense import ASO, Oligonucleotide


def test_aso_is_an_oligonucleotide():
    assert issubclass(ASO, Oligonucleotide)
