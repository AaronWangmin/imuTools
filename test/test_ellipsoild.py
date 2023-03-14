
import pytest

import sys

sys.path.insert(0,"E:\\project_201903\\ecarx-tools-202106\\imuTools")

import numpy as np

from src import unitTools
from src.ellipsoildUtil import Ellipsoild

def test_ellipsoid():
    GRS80 = Ellipsoild()
    
    assert(np.max(np.abs(np.subtract(GRS80.a,6378137))) < 0.01)
    
    assert(np.abs(np.subtract(GRS80.computeB(),6356752.3141)) < 0.01)     
    
    assert(np.max(np.abs(np.subtract(GRS80.computeE(),np.power(0.00669437999013,0.5)))) < 0.0001)      
    
    gravity = GRS80.computeGravity(unitTools.degree2radian(31.16559042328),10)
    gravity_check = np.array([0,0,9.794139999246411]) 
    assert(np.max(np.abs(np.subtract(gravity,gravity_check))) < 0.01)
   