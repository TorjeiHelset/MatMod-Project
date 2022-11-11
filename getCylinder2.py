# Description:
#   Generate a mesh triangulation of the unit disc.
#
# Arguments:
#   N       Number of nodes in the mesh.
#
# Returns:
#   p		Nodal points, (x,y)-coordinates for point i given in row i.
#   tri   	Elements. Index to the three corners of element i given in row i.
#   edge  	Edge lines. Index list to the two corners of edge line i given in row i.
#
#   Author: Kjetil A. Johannessen, Abdullah Abdulhaque
#   Last edit: 09-10-2019


import numpy as np
import scipy.spatial as spsa

def GetCylinder(N,n,h=0.1):
    # Controlling the input 
    if n < 3:
        print("Error. n >= 3 required for input.")
        return
    dh = h / (n-1) # Height of each layer
    p_, tri_ = GetDisc(N)
    
    p_first = np.zeros((len(p_), 3))
    index = np.arange(len(p_))
    p_first[index,0] = p_[index,0]
    p_first[index,1] = p_[index, 1]
    p_first[index,2] = np.zeros(len(p_))
    
    p = np.zeros((len(p_) * n ,3))
    p[index] = p_first[index]
    tri_first = np.zeros((len(tri_)*3, 4))
    for i in range(len(tri_)):
        v1 = tri_[i,0]
        v2 = tri_[i,1]
        v3 = tri_[i,2]
        v4 = v1 + N
        v5 = v2 + N
        v6 = v3 + N
        tri_first[i*3] = np.array([v1, v2, v3, v4])
        tri_first[i*3+1] = np.array([v2, v4, v5, v6])
        tri_first[i*3+2] = np.array([v2, v3, v4, v6])
    tri = tri_first    
    p_index = index
    index = np.arange(len(p_))
    for i in range(n-2):               
        # Go through each layer
        # For each 2D element in base, add 3 3D elements
        # Add z-component to nodes in each layer
        # Node directly above has index N higher
        p_new = np.zeros((len(p_), 3))
        p_new[index,0] = p_[index,0]
        p_new[index,1] = p_[index, 1]
        p_new[index,2] = np.ones(len(p_)) * dh*(i+1)
        p_index += N
        p[p_index] = p_new[index]

        tri_new = tri_first + N
        tri = np.concatenate((tri, tri_new))

    p_new = np.zeros((len(p_), 3))
    p_new[index,0] = p_[index,0]
    p_new[index,1] = p_[index, 1]
    p_new[index,2] = np.ones(len(p_)) * h
    p_index += N
    p[p_index] = p_new[index]
    
    M,r,alpha,theta = CircleData(N)
    BottomEdgeNodes, MiddleEdgeNodes, TopEdgeNodes = FreeBoundary(N,alpha, n)
    
    return p, tri.astype(int), BottomEdgeNodes, MiddleEdgeNodes, TopEdgeNodes
        
        

def GetDisc(N):
    # Controlling the input.
    if N < 4:
        print("Error. N >= 4 reguired for input.")
        return

    # Defining auxiliary variables.
    M,r,alpha,theta = CircleData(N)

    # Generating the nodal points.
    p = NodalPoints(M,N,alpha,theta,r)

    # Generating the elements.
    mesh = spsa.Delaunay(p)
    tri = mesh.simplices

    # Generating the boundary elements.
    #edge = FreeBoundary(N,alpha)

    return p,tri

def NodalPoints(M,N,alpha,theta,r):
    # Auxiliary function for generating nodal points.
    p = np.zeros((N,2))
    k = 1
    for i in range(1,M+1):
        t = theta[i]
        for j in range(0,alpha[i]):
            p[k,:] = [np.cos(t)*r[i],np.sin(t)*r[i]]
            t += 2*np.pi/alpha[i]
            k += 1

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