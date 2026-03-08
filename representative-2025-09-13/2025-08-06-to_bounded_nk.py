# p : Polyhedron
# a : -a <= x_i <= a
def to_bounded(p, a):
    l = p.inequalities_list()
    #print(l)
    np = len(l[0])
    n = np-1
    for i in range(1, np):
        t = [0]*np
        t[0] = a
        t[i] = 1
        l.append(t)
        #print(i, l)
        tt = [0]*np
        tt[0] = a
        tt[i] = -1
        l.append(tt)
        #print(i, l)
    pp = Polyhedron(ieqs=l)
    return pp

def test_to_bounded():
    l = [[0, -1, 0, 1], [0, 1, 0, 0], [0, 0, 1, -1], [1, 1, 0, -1]]
    p = Polyhedron(ieqs=l)
    pp = to_bounded(p, 5)
    return pp

# return interior integral points of p 
# p : Polyhedron
# a : -a <= x_i <= a
def ub_int_p(p, a):
    bp = to_bounded(p, a)
    l = bp.integral_points()
    interior_pts = []
    for i, pt in enumerate(l):
        if p.interior_contains(pt):
            interior_pts.append(pt)
    return interior_pts

"""
test data
gauss unbounded regions, representative_points
ur_index = [2,3,5,8,10,11,12,13,14,15,16,17,19,22,24,25]
0 (2, 2, 0)
1 (-1, 1, -2)
2 (1, 5/2, 3/2)
3 (-4/3, 4/3, -2/3)
4 (1, 4, 3)
5 (-5/2, 5/2, 1/2)
6 (1, -1, -2)
7 (-2, -2, -4)
8 (-5/2, -1, -2)
9 (-4, -1, -2)
10 (5/2, 1, 3/2)
11 (5/4, 5/4, 7/4)
12 (1, 5/2, 3)
13 (-4/3, 4/3, 5/3)
14 (4/3, -4/3, -2/3)
15 (-1, -5/2, -2)
16 (-5/4, -5/4, -3/4)
17 (-5/2, -1, -1/2)
18 (4, 1, 3)
19 (5/2, 1, 3)
20 (2, 2, 5)
21 (-1, 1, 3)
22 (5/2, -5/2, 1/2)
23 (-1, -4, -2)
24 (4/3, -4/3, 5/3)
25 (-1, -5/2, -1/2)
26 (1, -1, 3)
27 (-2, -2, 1)
"""

##############################################
# generate a hyperplane arrangement from bf 
##############################################

# e.g. p2list(a+2*b+3*c+4, [a,b,c]) ---> [(1,2,3),4]
def p2list(f, vl):
    l = []
    for t in vl:
        l.append(f.coefficient(t))
    return [tuple(l), f.constant_coefficient()]

def test_p2list():
    ring = PolynomialRing(QQ, 'a,b,c')
    a,b,c = ring.gens()
    bf = [(-a + c - 1, 1), (-a + c, 1), (-b + c - 1, 1), (-b + c, 1), (b, 1), (a, 1)]
    return p2list(bf[0][0], [a,b,c])

def bf2ha(bf, vl):
    l = []
    for t in bf:
        c = p2list(t[0], vl)
        l.append(c)
    svl = list(map(str, vl))
    H = HyperplaneArrangements(QQ, tuple(svl))
    ha = H(l)
    return ha

def test_bf2ha():
    ring = PolynomialRing(QQ, 'a,b,c')
    a,b,c = ring.gens()
    bf = [(-a + c - 1, 1), (-a + c, 1), (-b + c - 1, 1), (-b + c, 1), (b, 1), (a, 1)]
    return bf2ha(bf, [a,b,c])

##############################################
# return representative points of each regions of hyperplane arrangements
##############################################

# h : hyperplane arrangement (e.g., the output of bf2ha)
# output: list of elements s.t. [region, representative point or None]
def rep_pts(h):
    br = h.bounded_regions()
    ur = h.unbounded_regions()

    br_rep = []
    for t in br:
        exists = False
        pts = t.integral_points()
        for p in pts:
            if t.interior_contains(p):
                br_rep.append([t, p])
                exists = True
                break
        if exists == False:
            br_rep.append([t,None])

    ur_rep = []
    for t in ur:
        exists = False
        tt = to_bounded(t, 5)
        pts = tt.integral_points()
        for p in pts:
            if t.interior_contains(p):
                ur_rep.append([t, p])
                exists = True
                break
        if exists == False:
            ur_rep.append([t,None])
    return [br_rep, ur_rep]
    
def test_rep_pts():
    ring = PolynomialRing(QQ, 'a,b,c')
    a,b,c = ring.gens()
    bf = [(-a + c - 1, 1), (-a + c, 1), (-b + c - 1, 1), (-b + c, 1), (b, 1), (a, 1)]
    vl = [a,b,c]
    h = bf2ha(bf, vl)
    return rep_pts(h)

def bf2rep(bf, vl):
    return rep_pts(bf2rep(bf, vl))
    

msg="""
type test_rep_pts()
"""
print(msg)