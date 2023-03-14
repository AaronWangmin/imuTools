import sys
sys.path.insert(0,"E:\\project_201903\\ecarx-tools-202106\\imuTools")

import numpy as np

from src import unitTools

def test_angle():
    assert(unitTools.degree2radian(90) - np.pi/2 < 0.001)

    assert(unitTools.radian2degree(np.pi) - 180 < 0.001)