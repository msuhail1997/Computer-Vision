import numpy as np

from line_from_pts import line_from_pts
from line_intersection import line_intersection

def compute_ref_line_segment(ref, ll, lr, ur, ul):
    """ 
    This function finds the end points of the line segment
    where the ref is located on the field. The results will
    be used to plot the virtual line onto the field.

    Input:
        ref: np.array of size (3,) representing a point of ref on the field
        ll:  np.array of size (3,) representing lower left point of rectangle in frame
        lr:  np.array of size (3,) representing lower right point of rectangle in frame
        ur:  np.array of size (3,) representing upper right point of rectangle in frame
        ul:  np.array of size (3,) representing upper left point of rectangle in frame
        
    Returns:
        vanishing_pt: np.array of size (3,) representing scene vanishing point 
        top_pt:       np.array of size (3,) representing top of ref's line segment
        bottom_pt:    np.array of size (3,) representing bottom of ref's line segment

    """

    ##### STUDENT CODE START #####
    
    vanishing_pt = []
    top_pt = []
    bottom_pt = []

    ##### STUDENT CODE END #####
    
    return vanishing_pt, top_pt, bottom_pt