#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:36:25 2020

@author: jim53
"""

def GMread(Ground_Motion):
    import os
    import numpy as np
    THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
    GM1 = os.path.join(THIS_FOLDER, Ground_Motion)
    
    GM=[]
    f=open(GM1)
    
    eq = []
    for line in f.readlines():
        eq.append(float(line))
    
    f.close
    
    return np.array(eq)