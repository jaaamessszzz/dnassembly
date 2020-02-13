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


# --- Circular double linked lists for Part/Cassette Annotation --- #

class PartNode(object):
    """
    Node in circular double linked list
    """
    def __init__(self, part):
        self.part = part
        self.previous = None
        self.next = None


class CircularOrder(object):
    """
    Base class for circular double linked list
    """
    def __init__(self):
        self.head = None
        self.size = 0
        self.assemblelinkedlist()

    def add(self, part):
        new_part = PartNode(part)

        if self.head is None:
            self.head = new_part
            self.head.previous = self.head
            self.head.next = self.head

        else:
            tail = self.head.previous
            tail.next = new_part
            new_part.previous = tail
            new_part.next = self.head
            self.head.previous = new_part

        self.size += 1

    def assemblelinkedlist(self):
        for part in self.parts:
            self.add(part)


class PartOrder(CircularOrder):
    """
    Part plasmid annotation
    """
    parts = ['1', '2', '3a', '3b', '4a', '4b', '5', '6', '7', '8a', '8b']

    def __init__(self):
        super(PartOrder, self).__init__()


class CassetteOrder(CircularOrder):
    """
    Cassette plasmid annotation
    """
    parts = ['LS-1', '1-2', '2-3', '3-4', '4-5', '5-RE', 'RE-LS']

    def __init__(self):
        super(CassetteOrder, self).__init__()


class MoCloAssemblyType(Enum):
    PART = BsmBI
    CASSETTE = BsaI
    MULTICASSETTE = BsmBI


class AssemblyTypeException(Exception):
    pass
