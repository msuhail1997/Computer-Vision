import numpy as np
from est_homography import est_homography
from est_Pw import est_Pw
from numpy.linalg import inv


def PnP(Pc, Pw, K):
    """ 
    Solve Perspective-N-Point problem with collineation assumption, given correspondence and intrinsic

    Input:
        Pc: 4x2 numpy array of pixel coordinate of the April tag corners in (x,y) format
        Pw: 4x3 numpy array of world coordinate of the April tag corners in (x,y,z) format
        K:  3x3 numpy array for camera intrisic matrix (given in run_PnP.py)
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: 3x1 numpy array describing camera translation in the world (t_wc)
        
    """
    
    ##### STUDENT CODE START #####
    Pw = np.delete(Pw, -1, axis=1)
    H=est_homography(Pw,Pc)
    R = np.eye(3)
    t = np.zeros([3])
    K= inv(K)
    h=np.matmul(K,H)
    t=h[:,2]/np.linalg.norm(h[:,0])
    h = np.delete(h, -1, axis=1)
    h1=h[:,[0]]
    h2=h[:,[1]]
    h1=h1.transpose()
    h2=h2.transpose()
    y=np.cross(h1,h2)
    y=y.transpose()
    h=np.append(h, y, axis=1)
    U, S, V = np.linalg.svd(h)
    I=np.identity(3)
    final=np.matmul(U,V)
    I[2,2]=np.linalg.det(final)
    R=np.matmul(U,np.matmul(I,V))
    R=R.transpose()
    t=np.matmul(-R,t)

    ##### STUDENT CODE END #####

    return R, t
