#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" QCDLoop Python Wrapper"""

__authors__ = 'Stefano Carrazza, et at.'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@cern.ch'


import qcdloop as ql

mu2 = 1.7**2
p = []
m = [5.0]

td = ql.QCDLoop()

o = td.integral(mu2,m,p)
print(o)
