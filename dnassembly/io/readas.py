#! /usr/bin/env python3

from enum import Enum

from ..dna import DNA as _DNA
from ..dna import Plasmid as _Plasmid
from ..dna import Part as _Part

class ReadAs(Enum):
    DNA = _DNA
    Plasmid = _Plasmid
    Part = _Part