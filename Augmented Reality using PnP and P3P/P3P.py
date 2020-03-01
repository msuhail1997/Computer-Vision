import numpy as np
import math


def P3P(Pc, Pw, K=np.eye(3)):
    """
    Solve Perspective-3-Point problem, given correspondence and intrinsic

    Input:
        Pc: 4x2 numpy array of pixel coordinate of the April tag corners in (x,y) format
        Pw: 4x3 numpy array of world coordinate of the April tag corners in (x,y,z) format
        K:  3x3 numpy array for camera intrisic matrix (given in run_P3P.py)
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: 3x1 numpy array describing camera translation in the world (t_wc)

    """

    ##### STUDENT CODE START #####

    R = np.eye(3)
    t = np.zeros([3])
    points = (0, 1, 2)
    Pc = np.linalg.inv(K) @ np.hstack((Pc[points, :], np.ones((3, 1)))).T
    Pc = Pc.transpose()
    f=823.3
    P1=Pw[0, :2]
    P2 = Pw[1, :2]
    P3 = Pw[2, :2]
    a = math.sqrt(((P2[0] - P3[0]) ** 2) + (P2[1] - P3[1]) ** 2)
    b = math.sqrt(((P1[0] - P3[0]) ** 2) + (P1[1] - P3[1]) ** 2)
    c = math.sqrt(((P1[0] - P2[0]) ** 2) + (P1[1] - P2[1]) ** 2)
    j1 = Pc[0] / np.linalg.norm(Pc[0])
    j2 = Pc[1]/ np.linalg.norm(Pc[1])
    j3 = Pc[2]/ np.linalg.norm(Pc[2])
    alpha = np.dot(j2, j3)
    beta = np.dot(j1, j3)
    gamma = np.dot(j1, j2)
    # define coefficients of the 4th degree polynomial

    A4 = (((a * a - c * c) / (b * b) - 1) ** 2) - 4 * c * c *alpha *alpha / (b * b)
    A3 = 4 * ((a * a - c * c) * (1 - (a * a - c * c) / (b * b) * beta) - (1 - (a * a + c * c) / (b * b)) * alpha * gamma + 2 * c * c * alpha * alpha * beta / (b * b))
    A2 = 2 * ((a * a - c * c) ** 2 / (b ** 4) - 1 + 2 * (beta * beta * (a * a - c * c) ** 2 / (b ** 4)) + 2 * ((b * b - c * c) * alpha * alpha / (b * b)) - 4 * ((a * a + c * c) * alpha * beta * gamma / (b * b)) + 2 * (b * b - a * a) *gamma * gamma / (b * b))
    A1 = 4 * (-((a*a - c*c)/(b*b) * (1 + (a*a - b*b)/(c*c)) * beta) + (2 * a * a * gamma * gamma * beta / (b * b)) - ((1 - (a*a + c*c)/(b*b)) * alpha * gamma))
    A0 = (1 + (a * a - c * c) / (b * b)) ** 2 - (4 * a * a * gamma * gamma) / (b * b)
    # calculate real roots u and v

    roots = np.roots(np.array([A4, A3, A2, A1, A0]))
    uv = []
    for root in roots:
        if np.isreal(root):
            v = root.real
            u = (((-1 + (a * a - c * c) / (b * b)) * v * v) - (2 * (a * a - c * c) / (b * b) * beta * v) + 1 + (
                        a * a - c * c) / (b * b)) / (2 * (gamma - v * alpha))
            uv.append([u, v])
    # check for valid distances

    for x in uv:
        u = x[0]
        v = x[1]
        s1_sq = c * c / (1 + u * u - 2 * u * gamma)
        if s1_sq < 0:
            continue
        s1 = math.sqrt(s1_sq)
        s2 = u*s1
        s3 = v*s1
        if s2 > 0 and s3 > 0:
            break
    # calculate 3D coordinates in Camera frame

    p1 = s1*j1
    p2 = s2*j2
    p3 = s3*j3
    p1w = Pw[0]
    p2w = Pw[1]
    p3w = Pw[2]
    # Calculate R,t using Procrustes

    X = np.array([p1, p2, p3])
    Y = np.array([p1w, p2w, p3w])
    R, t = Procrustes(X, Y)
    ##### STUDENT CODE END #####
    return R, t


def Procrustes(X, Y):
    """
    Solve Procrustes: Y = RX + t

    Input:
        X: Nx3 numpy array of N points in camera coordinate (returned by your P3P)
        Y: Nx3 numpy array of N points in world coordinate
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: 1x3 numpy array describing camera translation in the world (t_wc)

    """

    ##### STUDENT CODE START #####

    X = X.transpose()
    Y = Y.transpose()
    a_cent = np.reshape(np.mean(Y, axis=1), (3, 1))
    b_cent = np.reshape(np.mean(X, axis=1), (3, 1))
    Y = Y - a_cent
    X = X - b_cent
    U, S, V = np.linalg.svd(np.matmul(X, Y.transpose()))
    V = V.transpose()
    U = U.transpose()
    I = np.identity(3)
    final = np.matmul(U, V)
    I[2, 2] = np.linalg.det(final)
    R = np.matmul(V, np.matmul(I, U))
    t = a_cent - np.matmul(R, b_cent)
    ##### STUDENT CODE END #####
    return R, t
