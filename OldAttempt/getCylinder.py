# Based on code created by Author: Kjetil A. Johannessen, Abdullah Abdulhaque
#   Last edit: 09-10-2019


#Description:
#   Generate a mesh triangulation of the cylinder of height h
#   with unit disc as base
#
# Arguments:
#   N       Number of nodes in each circle layer.
#   n       Number of circle layers
#   h       Height of cylinder
#
# Returns:
#   p		Nodal points, (x,y,z)-coordinates for point i given in row i.
#   tri   	Elements. Index to the four corners of element i given in row i.
#   BottomEdgeNodes  	Nodal points on boundary at z = 0.
#   MiddleEdgeNodes  	Nodal points on boundary at 0 < z < 1.
#   TopEdgeNodes  	Nodal points on boundary at z = 1.


import numpy as np
import scipy.spatial as spsa

    
def GetCylinder(N, n, h = 0.1):
    # Controlling the input.
    if N < 4:
        print("Error. N >= 4 reguired for input.")
        return
    if n<2:
        print("Error. n>= 2 required for input.")
        return

    # Defining auxiliary variables.
    M,r,alpha,theta = CircleData(N)

    # Generating the nodal points.
    p = NodalPoints(M,N,alpha,theta,r, n, h)
    
    # Generating the elements.
    mesh = spsa.Delaunay(p)
    tri = mesh.simplices

    # Generating the boundary elements.
    BottomEdgeNodes, MiddleEdgeNodes, TopEdgeNodes = FreeBoundary(N,alpha, n)

    return p,tri,BottomEdgeNodes, MiddleEdgeNodes, TopEdgeNodes

def NodalPoints(M,N,alpha,theta,r,n, h):
    # Auxiliary function for generating nodal points.
    z = np.linspace(0, h, n)
    p = np.zeros([])
    for l in range(n):
        #t_ = (theta[1] - theta[0]) * l%2
        p_ = np.zeros((N,3))
        p_[:,2] = z[l]
        k = 1
        for i in range(1,M+1):
            t = theta[i]
            t_ = np.pi/alpha[i] * l%2
            for j in range(0,alpha[i]):
                p_[k,:] = [np.cos(t+t_)*r[i]+np.random.uniform(0,0.0001),
                           np.sin(t+t_)*r[i]+np.random.uniform(0,0.0001),
                           z[l] + np.random.uniform(0,0.0001)]
                t += 2*np.pi/alpha[i]
                k += 1
        if l == 0:
            p = p_
        else:
            p = np.concatenate((p, p_), axis = 0)
    return p



def FreeBoundary(N,alpha, n):
    # Auxiliary function for generating boundary nodes.
    BottomEdgeNodes = np.arange(0,N)
    MiddleEdgeNodes = np.array([])
    TopEdgeNodes = np.arange(N*(n-1), N*n)
    for j in range(1,n-1):
        E = np.arange(N-alpha[-1]+1,N+1)
        edge = np.zeros((len(E),2),dtype=np.int)
        for i in range(0,len(E)):
            edge[i,:] = [E[i],E[i]+1]
        edge[-1,-1] = N-alpha[-1]+1
        edge -= 1
        edge += N*j
        edge = np.unique(edge)
        if (j == 1):
            MiddleEdgeNodes = edge
        else:
            MiddleEdgeNodes = np.concatenate((MiddleEdgeNodes, edge))
    MiddleEdgeNodes = np.unique(MiddleEdgeNodes)
    return BottomEdgeNodes, MiddleEdgeNodes, TopEdgeNodes


def CircleData(N):
    # Number of outward circles,excluding the origin.
    M = np.int(np.floor(np.sqrt(N/np.pi)))

    # Radius of the different circles.
    r = np.linspace(0,1,M+1)

    # Number of DOF in each circle.
    alpha_temp = np.floor((2*np.pi*M)*r)
    alpha = np.zeros(len(alpha_temp),dtype=np.int)
    for i in range(0,len(alpha_temp)):
        alpha[i] = np.int(alpha_temp[i])

    # Fine-tuning to get the right amount of DOF.
    alpha[0] = 1
    i = 1
    while sum(alpha) > N:
        if alpha[i] > 0:
            alpha[i] -= 1
        i += 1
        if sum(alpha[1:M]) == 0:
            i = M
        elif i > M:
            i = 1
    while sum(alpha) < N:
        alpha[-1] += 1

    # Creating the starting angle.
    theta = np.pi/alpha
    theta[0:len(alpha):2] = 0

    return M,r,alpha,theta
