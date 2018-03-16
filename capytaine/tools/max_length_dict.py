#!/usr/bin/env python
# coding: utf-8


from collections import OrderedDict


class MaxLengthDict(OrderedDict):
    """Dictionary with limited number of entries. When maximum size is
    reached, the oldest entry is dropped at the insertion of a new one.
    """
    def __init__(self, *args, max_length=1, **kwargs):
        assert isinstance(max_length, int)
        assert max_length >= 0
        self.__max_length__ = max_length
        OrderedDict.__init__(self, *args, **kwargs)

    def __setitem__(self, key, val):
        if key in self:
            del self[key]
        OrderedDict.__setitem__(self, key, val)
        if len(self) > self.__max_length__:
            self.popitem(last=False)
