#!/usr/bin/env python
"""Tests for shared ASOkai types."""

from typing import get_args

from ASOkai import Types
from ASOkai.Types import Scalar, Strand, TargetRegion


def test_types_namespace_is_available_from_package_root():
    assert Types.Scalar is Scalar


def test_scalar_contains_supported_scalar_values():
    for value in ("text", 1, 1.0, True, None):
        assert isinstance(value, Scalar)


def test_shared_domain_types_are_available_from_types_namespace():
    assert set(get_args(Strand)) == {"+", "-"}
    assert set(get_args(TargetRegion)) == {
        "exonic_only",
        "pre-mrna",
        "transcriptomic",
    }
