#!/usr/bin/env python
"""
Filename: src/ASOkai/utils/type_registrations.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: Registers external types for serialization with the Serializable class.
License: LGPL-3.0-or-later
"""
from Bio.Seq import Seq
from GenomeUtils.Genome import Locus

from .serializer import Serializable


def register_types() -> None:
    """
    Register all external types used in ASOkai for serialization.
    Call this once at application startup.
    """
    # Register Bio.Seq.Seq
    Serializable.register_type(
        Seq,
        'Bio.Seq.Seq',
        serialize=lambda s: {'sequence': str(s)},
        deserialize=lambda d: Seq(d['sequence'])
    )
    
    # Register GenomeUtils.Locus with flatten=True
    # Components (chr, start, end, strand) merge directly into parent dict
    Serializable.register_type(
        Locus,
        'Locus',
        serialize=lambda l: {
            'chr': l.chr,
            'start': l.start,
            'end': l.end,
            'strand': l.strand
        },
        deserialize=lambda d: Locus(
            chr=d['chr'],
            start=d['start'],
            end=d['end'],
            strand=d['strand']
        ),
        flatten=True
    )


# Auto-register when this module is imported
register_types()
