#! /usr/bin/env python3

from enum import Enum
from Bio.Restriction import BsaI, BsmBI

from .goldengate import GoldenGate

class ModularCloning(GoldenGate):
    """
    The MoClo system is just Golden Gate assembly with extra rules
    """

    def __init__(self, parts, type):
        super(ModularCloning, self).__init__(parts)
        self.restriction_enzyme_list = type
        self.assembly_type = type

    @property
    def restriction_enzyme_list(self):
        return self._restriction_enzyme_list

    @restriction_enzyme_list.setter
    def restriction_enzyme_list(self, type):
        if isinstance(type, MoCloAssemblyType):
            self._restriction_enzyme_list = [type.value]
        else:
            raise AssemblyTypeException('Assembly type must be set using the MoCloAssemblyType class!')


class MoCloAssemblyType(Enum):
    PART = BsmBI
    CASSETTE = BsaI
    MULTICASSETTE = BsmBI


class AssemblyTypeException(Exception):
    pass