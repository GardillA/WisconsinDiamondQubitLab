# -*- coding: utf-8 -*-
"""
Interface for vector signal generators

Created on August 29th, 2022

@author: mccambria
"""

from abc import ABC, abstractmethod
from sig_gen import SigGen


class VectorSigGen(SigGen):
    
    @abstractmethod
    def load_iq(self, c):
        """
        Set up IQ modulation controlled via the external IQ ports
        """
        pass