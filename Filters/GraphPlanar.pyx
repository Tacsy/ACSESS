# PYREX code that creates an interface between C-language graph theory
# embedding routines and the rest of the python code
#######################################################################


# C structures from graph.h
# We declare only what we will need
cdef extern from "graph.h":
    ctypedef struct BM_Graph:
        pass
    ctypedef BM_Graph *graphP

    graphP gp_New()
    int    gp_InitGraph(graphP theGraph, int N)
    void   gp_Free(graphP *pGraph)
    int    gp_AddEdge(graphP theGraph, int u, int ulink, int v, int vlink)
    int    gp_Embed(graphP theGraph, int embedFlags)


#
# Convert OEChem molecule to graph, check planarity
def IsPlanar(mol):
    cdef graphP thegraph
    cdef int err,N,node1,node2,result

    # Macros from the header file - this works fine
    cdef EMBEDFLAGS_PLANAR,OK,NONPLANAR
    EMBEDFLAGS_PLANAR=1
    OK=0
    NONPLANAR=-3
    
    thegraph=gp_New()
    N=mol.GetNumAtoms()
    #N=mol.GetMaxAtomIdx()
    err=gp_InitGraph(thegraph,N)

    for bond in mol.GetBonds():
        node1=bond.GetBeginAtomIdx()
        node2=bond.GetEndAtomIdx()
        err=gp_AddEdge(thegraph,node1,0,node2,0)

    result=gp_Embed(thegraph,EMBEDFLAGS_PLANAR)

    gp_Free(&thegraph)

    if result==-3:
        return False
    else: return True
