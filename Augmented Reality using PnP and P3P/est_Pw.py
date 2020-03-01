import numpy as np


def est_Pw(s):
    """ 
    Estimate the world coordinates of the April tag corners, assuming the world origin
    is at the center of the tag, and that the xy plane is in the plane of the April 
    tag with the z axis in the tag's facing direction. See world_setup.jpg for details.

    Input:
        s: side length of the April tag
        
    Returns:
        Pw: 4x3 numpy array describing the world coordinates of the April tag corners
            in the order of a, b, c, d for row order. See world_setup.jpg for details.
        
    """
    
    ##### STUDENT CODE START #####


    Pw = np.zeros([4, 3])
    Pw[0,0]= -s/2
    Pw[0, 1] = -s / 2
    Pw[0, 2] = 0
    Pw[1, 0] = s / 2
    Pw[1, 1] = -s / 2
    Pw[1, 2] = 0
    Pw[2, 0] = s / 2
    Pw[2,1]= s/2
    Pw[2, 2] = 0
    Pw[3,0]= -s/2
    Pw[3,1]= s/2
    Pw[3,2]= 0
    ##### STUDENT CODE END #####

    return Pw

