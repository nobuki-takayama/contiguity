#
# Implementation of the algorithm by Micciancio (2007)
#
import itertools

def add_zero_rows(a):
    h=a.nrows()
    n=a.ncols()
    if n>h:
        z=zero_matrix(n-h,n)
        b=block_matrix([[a],[z]])
        return matrix(list(b))
    elif n==h:
        return a
    else:
        return add_zero_rows(a.T).T
"""
M=matrix([[1,2,3],[4,5,6]])
print(add_zero_rows(M))

M2=matrix([[1,2,3,2],[4,5,6,5]])
print(add_zero_rows(M2))
"""

#
# B and H are given in Micciancio. It returns tilde B.
#
#  span_ZZ(BT[0], BT[1], ...) cap { x | H[0]*x=0, H[1]*x=0, ...}
#  BT[i] and H[j] must be the same size m.
#
def lattice_basis(BT,H):
#    print('lattice_basis(\n ',BT,',\n ',H,')') #for debug
    BT=matrix(BT); H=matrix(H)
    m=BT.ncols()
    B=BT.T
    hb=H*B
    c=add_zero_rows(hb)
    s=c.smith_form(integral=True)
    v=s[1]
    u=s[2]
    n=s[0].ncols()
    dd=n
    for i in range(n):
        if s[0][i][i]==0:
            dd=i
            break
    d=n-dd
    tb=B*u
    tb2=list(reversed(list(tb.T)))
    ans=[]
    if d==0:
        return [[0]*m]
    for i in range(d):
        ans.append(tb2[i])
    return ans

"""
L=lattice_basis(matrix([[1,0,0],[0,1,0]]),matrix([1,1,1]))
print(L)
L=lattice_basis(matrix([[1,0,0],[0,1,0]]),matrix([[1,1,1],[1,2,3]]))
print(L)
L=lattice_basis(matrix([[1,0,0,0],[0,0,1,0]]),matrix([1,2,3,4]))
print(L)
L=lattice_basis(matrix([[1,0,0,0],[0,0,1,0],[0,0,0,1]]),matrix([[1,1,1,1],[1,2,3,4]]))
print(L)
L=lattice_basis(matrix([[1,0,0,0],[0,0,1,0],[0,0,0,1]]),matrix([[1,1,1,1]]))
print(L)
"""
#
# BT is a basis of lattice space in QQ^m
# AH  defines an affine space; AH_i=(h,h_m), h in QQ^m, hx+h_m=0
# S is a point in BT cap V(AH)
#
def affine_lattice_basis(BT,AH,S):
    BT=matrix(BT); AH=matrix(AH); S=vector(S)
    m=BT.ncols()
    n=BT.nrows()
    if S[m] != 1:
        print("Input error: The last element of S must be 1. Return 0")
        return 0
    BT2=matrix(list(block_matrix([[BT.T],[zero_matrix(1,n)]]))).T
    unit=zero_matrix(1,m+1); unit[0,m]=1
    B=matrix(list(block_matrix([[BT2],[unit]]))).T
    #return B
    hb=AH*B
    c=add_zero_rows(hb)
    s=c.smith_form(integral=True)
    v=s[1]
    u=s[2]
    n=s[0].ncols()
    dd=n
    for i in range(n):
        if s[0][i][i]==0:
            dd=i
            break
    d=n-dd
    tb=B*u
    tb2=list(reversed(list(tb.T)))
    ans=[]
    if d==0:
        return [vector([0]*m)]
    for i in range(d):
        ans.append(tb2[i])
    tt=last_one_subspace(ans,S)
    return remove_last_coord(tt)

#
#  BT の Z-span で 最後の座標が 1 であるものを求める.
#   S はこの Z-span の元で最後の座標が 1.
#   戻り値 + S が答え.
#
def last_one_subspace(BT,S):
#    print('last_one_subspace(',BT,')') #for debug
    BT=matrix(BT)
    m=BT.ncols()-1
    n=BT.nrows()
    H=matrix(1,m+1); H[0,m]=1
    ans=lattice_basis(BT,H)
    return ans

def remove_last_coord(H):
    n=len(H)
    ans=[]
    for i in range(n):
        v=list(H[i])[0:-1]
        ans.append(v)
    return ans

"""
L=affine_lattice_basis(matrix([[1,0,0],[0,0,1]]),matrix([1,2,3,4]),vector([-4,0,0,1]))
print(L)

# Example 1.
# a-c+1=0 in Z^3 の lattice basis. homogenize して [1,0,-1,1] の直交補空間
# <[1,0,-1,1],[0,0,1,1]>=0   
L=affine_lattice_basis(matrix([[1,0,0],[0,1,0],[0,0,1]]),matrix([1,0,-1,1]),vector([0,0,1,1]))
print(L)
# p1*(0,1,0,0)+p2*(1,0,1,0)+(0,0,1,1), p1, p2 in Z, の最初の３つの成分が a-c+1=0 を満たす.
"""

"""
# Example 2.
# (a,b,c)=p1*(0,1,0)+p2*(1,0,1)+(0,0,1), p1, p2 in Z と b-c+1=0 の intersection.
# まず c->c'+1 と変数変換. b-c+1=b-c' および  (a,b,c')=p1*(0,1,0)+p2*(1,0,1) となる.
L=affine_lattice_basis(matrix([[0,1,0],[1,0,1]]),matrix([0,1,-1,0]),vector([1,1,1,1]))
print(L)
# [(1,1,1,0)] が戻るので,
# (a,b,c')=p1*(1,1,1,0)+(1,1,1,1), p1 in Z, の最初の３つの成分が b-c'=0, a-c+1=0 を満たす lattice points 全体.
"""

# 2025.07.18
#

def mycons(a,b):
    return [a]+b
#  x_1+...+x_n=s, x_i >= 0
def nsum_is_s(n,s):
    if n==0: return []
    if n==1: return [[s]]
    sub=[]
    for i in range(s+1):
        sub.append(nsum_is_s(n-1,s-i))
    ans=[]
    for i in range(s+1):
        for j in sub[i]:
            ans.append(mycons(i,j))
    return ans
# print(nsum_is_s(3,4))

# list_sign_patterns([a,0,c])
# [a,0,c]->[[a,0,c],[a,0,-c],[-a,0,c],[-a,0,-c]]
def list_sign_patterns(v):
    n=len(v)
    if n==0: return []
    non_zero=0;
    for i in range(n):
        if v[i]!=0: non_zero += 1
    if non_zero==0: 
        sign=[[0]*n]
    else:
        sign=itertools.product(range(-1,2,2),repeat=non_zero)
    ans=[]
    for s in sign:
        new_s=v.copy()
        j=0;
        for i in range(n):
            if new_s[i]!=0:
                new_s[i] *= s[j]
                j += 1
        ans.append(new_s)
    return ans
# print(list_sign_patterns([2,0,3]))            

#  |x_1|+...+|x_n|=s 
def abs_nsum_is_s(n,s):
    x_list=nsum_is_s(n,s)
    ans=[]
    for x in x_list:
        ans=ans+list_sign_patterns(x)
    return ans
#abs_nsum_is_s(2,2)
#print(abs_nsum_is_s(3,4))

##Ref: goodnotes/2025-07-18-my-note-lattice.*
# BT: lattice basis of L
# AH=0, INEQ >0 をみたす L の点を一つ見つける. しらみつぶし.
# limit を大きくすると検索範囲が広がる. 
# BT.T*y,  |y|<limit
#
def find_ip_sol1(BT,AH,INEQ,limit=4):
    Debug=False
#    Debug=True
    BT=matrix(BT); AH=matrix(AH); INEQ=matrix(INEQ);
    m=BT.ncols()
    n=BT.nrows()
    BT2=matrix(list(block_matrix([[BT.T],[zero_matrix(1,n)]]))).T

    if len(list(AH))==0: AH=matrix([0]*(m+1))
    if len(list(INEQ))==0: INEQ=matrix([0]*(m+1)); INEQ[0,m]=1

    Zero_cond=[AH*(BT2.T),list(AH)[0][-1]]
    #return AH*(BT2.T)
    Ineq_cond=[INEQ*(BT2.T),matrix(list((INEQ.T))[-1])]
    #return INEQ*(BT2.T)
    #return [Zero_cond,Ineq_cond]
    ans=[]
    for s in range(0,limit+1):
        y_list=abs_nsum_is_s(n,s)
        for y in y_list:
            if Debug: print('y=',y)
            y=matrix(y).T
            if Zero_cond[0]*y+Zero_cond[1] != 0:
                continue
            success=1
            if Debug: print('Ineq_cond[0]=\n',Ineq_cond[0])
            vec0=Ineq_cond[0]*y; vec1=Ineq_cond[1]
            if Debug: print('y=',y,'\nvec0=\n',vec0,'\nvec1=',vec1)
            vec=vec0.T+vec1
            if Debug: print('list(vec[0])=',list(vec[0]))
            for pp in list(vec[0]):
                if pp < 0:
                    success=0; continue
            if success==1:
                return [y.T,BT.T*y]
    return []

# Example: テスト入力
"""
ring=PolynomialRing(QQ,'a1,a2,a3,a4,c1,c2,c3,c4,x1,x2')
a1,a2,a3,a4,c1,c2,c3,c4,x1,x2=ring.gens()
f=find_ip_sol1([[1,2,3],[4,5,6]],[a1,a2,a3,a4],[[c1,c2,c3,c4],[a1,a2,a3,a4]])
print(f)
"""


# Example: テスト入力
# Lattice = y1*[1,0,0]+y2*[0,1,0] =:[x1,x2,x3]
# 1*x1+1*x2+1*x3-1=0
# 1*x1+0*x2+0*x3+0>0, 0*x1+1*x2+1*x3+0>0
# |y1|+|y2| <= 4 
# ans: f[0]=[y1,y2], f[1]=[x1,x2,x3] 
f=find_ip_sol1([[1,0,0],[0,1,0]],[1,1,1,-1],[[1,0,0,0],[0,1,1,0]],limit=4)
ring=PolynomialRing(QQ,'a1,a2,a3,a4,c1,c2,c3,c4,x1,x2')
a1,a2,a3,a4,c1,c2,c3,c4,x1,x2=ring.gens()
#f=find_ip_sol1([[1,0,0],[0,1,0]],[1,1,1,-1],[[a1,a2,a3,a4],[c1,c2,c3,c4]],limit=4)
print(f)

# Example: テスト
# Lattice = y1*[1,2,3]+y2*[1,1,1] =:[x1,x2,x3]
# 2*x1+1*x2+1*x3-11=0
# -1*x1+1*x2+0*x3+0>0, 2*x1-1*x2+0*x3+0>0
# |y1|+|y2| <= 2
# ans: f[0]=[y1,y2], f[1]=[x1,x2,x3] 
f=find_ip_sol1([[1,2,3],[1,1,1]],[2,1,1,-11],[[-1,1,0,0],[2,-1,0,0]],limit=2)
print(f)

# 2025.07.20
# (lattice BT) + S と AH=0 の交わりの表現 lattice + T を戻す
#  T の検索には find_ip_sol1(BT,AH(y+s),limit) を用いる.
def affine_lattice_basis2(BT,S,AH,limit=3,T=0):
    Debug=False
    BT=copy(list(BT)); S=copy(list(S)); AH=copy(list(AH))
    BT=matrix(BT)
    m=BT.ncols()
    n=BT.nrows()
    AH=matrix(AH)
    S=vector(list(S))
    S0=matrix(list(S)+[0])
    AH_new=copy(AH)
    if Debug: print('AH*S0.T=',(AH*S0.T)[0,0])
    AH_new[0,-1]=AH_new[0,-1]+(AH*S0.T)[0,0]
# AH_new が AH(y+s) の y の多項式としての係数.
    #return AH_new
##Ref: 2025-07-18-my-note-lattice-Tokyo.goodnotes 
# 2025.07.21
    if T==0:
        Tp=find_ip_sol1(BT,AH_new,[],limit=limit)
        if len(list(Tp))==0:
            print('Warning! find_ip_sol1 failed with [BT,AH_new,[],limit]=',[BT,AH_new,[],limit])
            print('  Make limit larger if necessary.')
            return [[],[]]  # 2025-09-05
        # みつからない時は [] が戻るので次でエラー.
        Tp=Tp[1]
        Tp=Tp.T
        Tp=Tp[0]
    else:
        Tp=vector(list(T))-S
    #return Tp
    #return [BT,list(AH_new[0]),list(Tp)+[1]]
    L2=affine_lattice_basis(BT,list(AH_new)[0],list(Tp)+[1])
    return [L2,Tp+S]

#f=affine_lattice_basis2([[1,0,1],[0,1,0]],[1,2,3],[2,1,1,-4],T=[0,2,2])
#f=affine_lattice_basis2([[1,0,1],[0,1,0]],[1,2,3],[2,1,1,-4])

f=affine_lattice_basis2([[1,0,0],[0,1,0],[0,0,1]],[1,2,3],[2,1,1,-4],limit=5)
