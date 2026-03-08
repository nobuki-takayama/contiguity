def e_i_connected(x,y):
    d=len(x)
    z=vector(x)-vector(y)
    sum=0
    for i in range(d):
        sum += abs(z[i])
    if sum > 1: return False
    return True

def construct_mesh_edges(L):
    if L==[]: return [True,None]
    d=len(L[0])
    edges={}
    for i in range(len(L)):
        x=L[i]
        i_adjacent=[]
        for j in range(i+1,len(L)):
            y=L[j]
            if e_i_connected(x,y): i_adjacent.append(j)
        edges[i]=i_adjacent
    return edges

def check_if_mesh_type(L):
    g=Graph(construct_mesh_edges(L))
    return [g.is_connected(),g]

r"""
p=construct_mesh_edges([[1,1,1],[2,1,1],[1,2,1],[2,2,1]])
print(p)
print(Graph(p).is_connected())
pp=check_if_mesh_type([[1,1,1],[2,1,1],[1,2,1],[2,2,1]])
print(pp)
p2=construct_mesh_edges([[1,1,1],[2,1,1],[1,2,1],[2,2,3],[2,2,1]])
print(p2)
print(Graph(p2).is_connected())
p3=construct_mesh_edges([[1,1,1],[1,2,1],[2,2,2],[2,1,1]])
print(p3)
print(Graph(p3).is_connected())
"""
