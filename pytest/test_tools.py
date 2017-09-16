#!/usr/bin/env python
# coding: utf-8

from capytaine.tools import MaxLengthDict

def test_MaxLengthDict():
    dc = MaxLengthDict({'a':1, 'b':5, 'c':3}, max_length=4)
    assert dc.__max_length__ == 4
    assert list(dc.keys()) == ['a', 'b', 'c']
    dc['d'] = 8
    assert list(dc.keys()) == ['a', 'b', 'c', 'd']
    dc['e'] = 2 # drop 'a'
    assert list(dc.keys()) == ['b', 'c', 'd', 'e']
    dc['b'] = 7 # move 'b' to front
    assert list(dc.keys()) == ['c', 'd', 'e', 'b']
    dc['f'] = 4 # drop 'c'
    assert list(dc.keys()) == ['d', 'e', 'b', 'f']
    dc.update({'g':6, 'h':9})
    assert list(dc.keys()) == ['b', 'f', 'g', 'h']

    dc2 = MaxLengthDict({'a':1, 'b':5, 'c':3}, max_length=1)
    assert dc2.__max_length__ == 1
    assert list(dc2.keys()) == ['c']

    dc3 = MaxLengthDict({'a':1, 'b':5, 'c':3}, max_length=0)
    assert dc3 == {}
    dc3['d'] = 8
    assert dc3 == {}
