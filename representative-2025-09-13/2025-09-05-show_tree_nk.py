# load("2025-07-23-representative.py")

##############################################
# from file "2025-08-06-to_bounded.py"
##############################################

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

# bounded polyhedron version
def b_int_p(p):
    l = p.integral_points()
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
    return rep_pts(bf2ha(bf, vl))


##############################################
##############################################
    
"""
l はrepresentative の返り値 
e.g. 
l = load("sb-F2-representative.sobj") 

l[0] : 関数情報
l[1] : 隣接関係式の情報 
    l[1][1] : bf 
l[2] 以降 l[1] の結果に対応して，制限した部分に対する計算結果(l[0],l[1],l[2] と同様の構造を持つ) 

["0-dim face", 点] が木の節点 
"""
def show_tree(l, depth):
    n = len(l)
    if l[0]=="0-dim face":
        print("    "*depth, end="")
        print(l)
        return
    print("    "*depth, end="")
    print(l[0][1][0], l[0][1][1])
    print("    "*depth, end="")
    print("bf", l[1][1])  # bf 

    # 代表点の計算
    bf = l[1][1]
    h = bf2ha(bf, bf_vars(bf))
    pts = rep_pts(h)
    print("    "*depth, end="")
    print("pts")
    # print(pts) 
    show_rep_pts(pts, l[0][1][0], l[0][1][1], depth)

    for i in range(2, n):
        print("    "*(depth+1), end="")
        print(l[1][1][i-2])     # bf の因子 l[i] に対応する
        fctr=l[1][1][i-2] 
        # l[1][0] 隣接関係式情報  
        # l[1][0][0][0] パラメータ情報
        # l[1][0][0][1] bf 情報 
        # fctr in l[1][0][k][1] なる k を探し，パラメータ情報を表示 
        for k in range(len(l[1][0])):
                if fctr in l[1][0][k][1]:
                    break
        print("    "*(depth+1), end="")
        print(l[1][0][k][0])     # パラメータ情報
        """
        print("    "*(depth+1), end="")
        print(l[1][0][k][2])     # cont. 
        print("    "*(depth+1), end="")
        print(l[1][0][k][3])     # cont.
        """
        show_tree(l[i], depth+1)

# l は多項式を因数分解した結果のリスト --> 変数リストを返す
# e.g. [(-a + c - 1, 1), (-a + c, 1), (-b + c - 1, 1), (-b + c, 1), (b, 1), (a, 1)]
def bf_vars(l):
    s = set()
    for t in l:
        s = s.union(t[0].variables())
    return list(s)


# pts : rep_pts の返り値
# l, c: 元パラメータと新パラメータの間の変換のためのリスト
# d : 再帰の深さ(表示用) 
# e.g., (a,b,c) = a'[0, 1, 0]+b'[1, 0, 1]+(0, 0, 1) 
# l=[[0,1,0],[1,0,1]], c=[0,0,1]
def show_rep_pts(pts, l, c, d):
    ur = pts[0]
    br = pts[1]
    n = len(c) # パラメータの個数 
    # unbounded region
    for t in ur:
        pd(d); print(t[0].Hrepresentation()) # 領域の定義方程式
        pd(d); print(t[1])                   # 領域の点
        if t[1] == None: 
            continue
        # 元のパラメータへの変換
        param = [0] * n
        for i in range(n):
            param[i] = c[i]
            for j in range(len(t[1])):
                param[i] += t[1][j] * l[j][i]
        pd(d); print(param)
                
    # bounded region
    for t in br:
        pd(d); print(t[0].Hrepresentation()) # 領域の定義方程式
        pd(d); print(t[1])                   # 領域の点
        if t[1] == None: 
            continue
        # 元のパラメータへの変換
        param = [0] * n
        for i in range(n):
            param[i] = c[i]
            for j in range(len(t[1])):
                param[i] += t[1][j] * l[j][i]
        pd(d); print(param)

"""
Memo: 
元のパラメータと新しいパラメータでは変換が行われる
(a,b,c) = a'[0, 1, 0]+b'[1, 0, 1]+(0, 0, 1) 
a = b'
b = a'
c = b'+1
"""

# print_depth(d) 
def pd(d):
    print("    "*d, end="")

