# -*- coding: utf-8 -*-
"""
Created on Thu May  5 17:12:58 2022

@author: Dean
"""

def ppmv_to_gkg(q):
    """
    Converts specific humidity from units of ppmv to g/kg.
    """
    return 0.622*q/1000