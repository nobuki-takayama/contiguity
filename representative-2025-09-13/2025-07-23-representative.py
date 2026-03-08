load("2025-07-18-a2s.py");
load("2025-07-07-lattice-base.py");
import traceback
import sys
# 2025-07-20, start

# global. contiguity の動きについての strategy
# 0:  a<->a+1 のみ. 1: a<->a+1, a<->a-1 (平行な bf が現れる)
# 2025-07-27.  一般には Markov basis に対応する動きで.
Rep_strategy=0   

def set_strategy(s):
    global Rep_strategy
    if s<0:
        return Rep_strategy
    Rep_strategy=s
    return Rep_strategy

r"""
Examples:
L=gauss_contiguity_by_GKZ('[a+b,b,b]','[a+1+b+1,b+1,b+1]')

H=affine_lattice_basis([[0,1,0],[1,0,1]],[0,1,-1,0],[1,1,1,1])
  [1,1,1,1] 最後は1で, (lattice basis) \cap AH=0 を満たす点. 戻り値はこの点を起点としたこの intersection の lattice basis.
  戻り値 H=[[1,1,1]]
  remove_last_coord で最後の座標を消去.

               lattice gen in Q^3, 2*x1+x2+x3-11=0   
f=find_ip_sol1([[1,2,3],[1,1,1]],[2,1,1,-11],[],limit=2)
ノート.
Gauss の例: m=3, Z^m の座標は (a,b,c)
Z^m 内の V は lattice base + S(V) と表現されてる.
   たとえば
   base=[[1,0,1],[0,1,1]], S(V)=[2,1,1] とする.
   この時 contiguity を求める式は
   (2,1,1)+(a,b,a)->(2,1,1)+(a+1,b,a+1)
   (2,1,1)+(a,b,b)->(2,1,1)+(a,b+1,b+1)
      (2,1,1) は略してよい.
   0 の位置の変数は独立. 0 でない位置の変数は同じに.

   b関数の一次因子を [a0,a1,...,am] とする.
        L(x)=(a0*x1+a1*x2+...+am)
   L(x)=0 と V の交わりの lattice base と参照点を求める.

"""

##Ref: 2025-07-18-my-note-lattice-Tokyo.goodnotes
# affine_lattice_basis2(BT,S,AH,limit=3,T=0) in 2025-07-07-lattice-base.py
#  (lattice BT)+S と AH=0 の交わりを [L2,T] (lattice L2)+T) と表現.
# T は上記交わりの点で, 予め与えるか 0 なら | | が limit までの範囲で探索.
# Ex. f=affine_lattice_basis2([[1,0,1],[0,1,0]],[1,2,3],[2,1,1,-4],T=[0,2,2])
#  AH=[2,1,1,-4] は 2*x1+x2+x3-4 の意味.
# Ex2. f=affine_lattice_basis2([[1,0,1,2],[0,1,1,0]],[1,0,0,1],[2,1,1,1,-4],limit=5)
# 2025.07.23

def gauss_contiguity(old,new,r):
    return gauss_contiguity_by_GKZ(old,new,ring=r)
def f1_contiguity(old,new,r):
    return f1_contiguity_by_GKZ(old,new,ring=r)
def f2_contiguity(old,new,r):
    return f2_contiguity_by_GKZ(old,new,ring=r)
def f0134_contiguity(old,new,r):
    return f0134_contiguity_by_GKZ(old,new,ring=r)
# 注意: 整数係数多項式のみ.
# Ex. linear_poly_coefficients(a-c,[a,b,c],[c,a,a])
# linear_poly_coefficients は 2025-07-18-a2s.py へ. 2025.07.25
def linear_poly_coefficients_old1(f,vars,new_vars):
    #return f.coefficients(sparse=False)
    m=len(vars)
    clist=[0]*(m+1)
    h=0
    for i in range(m):
        try:
            pos=new_vars.index(vars[i])
        except ValueError:
            pos=-1
        if pos>=0:
            clist[pos]=int(f.coefficient(vars[i]));
            h += clist[pos]*vars[i]
         
    try:
        clist[m]=int(f-h)
    except TypeError:
        print('Error: f=',f,', h=',h,', [vars,new_vars]=',[vars,new_vars])
        try:
            raise Exception
        except:
            traceback.print_exc()
            print('Type in ret to exit')
            input()
            sys.exit(1)

    return clist
def is_zero(vec):
    vec=list(vec)
    n=len(vec)
    for i in range(n):
        if vec[i]!=0: return False
    return True

# (lattice BT)+S で conti_func を用いて 代表元を列挙する. vars は使うパラメータ
#  conti_func の例: gauss_contiguity_by_GKZ(old,new,ring=vars[0].parent())  
#
def representative(BT,S,conti_func,vars):
    # 2025-07-27.  一般には Markov basis に対応する動きで.
    global Rep_strategy  # 0: old <-> up のみ.  1: old <-> up と old <-> down 両方. 

    old_new_list = set_old_new_param_list(BT,S,vars)
    n=len(old_new_list)
    msg=['representative [BT,S,conti_func,vars]',[BT,S,conti_func,vars]]
    print(msg)
    Conti=[msg]
    #todo. BT をマニュアル的に作り変える関数を呼ぶ.
    bf=vars[0]**0
    conti_set=[]
    for i in range(n):
        old=old_new_list[i][0]
        new_up=old_new_list[i][1]
        new_down=old_new_list[i][2]  # new_down is used only for Rep_strategy=1
        #return[[old,new_up],vars[0].parent()]
        print('calling ',[old,new_up])
        conti_up=conti_func(old,new_up,vars[0].parent())
        print('Done\nCalling ',[new_up,old])
        conti_down=conti_func(new_up,old,vars[0].parent())
        print('Done')
        bf_this=lcm(conti_up[1],conti_down[1])
        bfi=list(factor(bf_this))
        bf=lcm(bf,bf_this)
        conti_set.append([[old,new_up],bfi,conti_up[0:2],conti_down[0:2]])
        if Rep_strategy==1:
            [conti_down2,conti_up2,bf2,bfi2]=find_contiguity_and_bf(old,new_down,conti_func,vars)
            conti_set.append([[old,new_down],bfi2,conti_down2[0:2],conti_up2[0:2]])
            bf=lcm(bf,bf2)
    # breadth-first search is done.
    bf=list(factor(bf))
    print(" [latticce BT,origin S,param]=",[BT,S,old]," と bf=",bf,"で代表点 Rep を探す. Conti に代表点 Rep を加える (todo)")
    Rep=0  # todo
    Conti.append([conti_set,bf,Rep])
    # bf の各因子について lattice, origin を求めて representative を求める. Conti に戻り値を append
    for j in range(len(bf)):
        print('j=',j,', [bf[j][0],vars,old]=',[bf[j][0],vars,old])
        AH=linear_poly_coefficients(bf[j][0],vars,old)
        print('[BT,S,AH]=',[BT,S,AH]);
        BT_new,T=affine_lattice_basis2(BT,S,AH)
        print('BT_new=',BT_new,' T=',T)
        if (len(BT_new)==0): continue  # 2025-09-05
        # b=0 との交わりに変化がない場合.
        if (list(BT)==list(BT_new)) and (list(S)==list(T)):
            msg=['intersection with b=0 does not shirink, [BT,S,AH]',[BT,S,AH]]
            print(msg)
            continue   # fixed 2025.09.04
        # 交わりが 0 次元の場合.
        if (len(BT_new)==1) and is_zero(BT_new[0]):
            msg=['0-dim face',T,bf[j]]
            print(msg)
            Conti.append(msg)
            continue
        Conti.append(representative(BT_new,T,conti_func,vars))
    return Conti

# new_up は　new_down でもいい.  2025-07-27
def find_contiguity_and_bf(old,new_up,conti_func,vars):
    print('calling ',[old,new_up])
    conti_up=conti_func(old,new_up,vars[0].parent())
    print('Done\nCalling ',[new_up,old])
    conti_down=conti_func(new_up,old,vars[0].parent())
    print('Done')
    bf_this=lcm(conti_up[1],conti_down[1])
    return [conti_up,conti_down,bf_this,list(factor(bf_this))]


# [old, new_up, new_down]
# Ex. dir=[1,0,0,2]
def set_old_new_param_0(BT,i,vars,old,S):
    n=len(list(BT))
    new_up=0
    new_down=0
    for j in range(n):
        if i==j:
            new_up += (vars[j]+1)*vector(BT[j])
            new_down += (vars[j]-1)*vector(BT[j])
        else:
            new_up += vars[j]*vector(BT[j])
            new_down += vars[j]*vector(BT[j])
    new_up += vector(list(S))
    new_down += vector(list(S))
    return [list(old),list(new_up),list(new_down)]

def set_old_new_param_list(BT,S,vars):
    n=len(list(BT))
    old=0
    for i in range(n):
        old += vars[i]*vector(BT[i])
    old += vector(list(S))
    L=[]
    for i in range(n):
        old_new=set_old_new_param_0(BT,i,vars,old,S)
        L.append(old_new)
    return L

# Example: Gauss    
ring=PolynomialRing(QQ,'x,dx,a,b,c')
x,dx,a,b,c=ring.gens()
dic={'x':x,'dx':dx,'a':a,'b':b,'c':c}

#f=set_old_new_param_list([[1,0,0,2],[0,1,1,1]],[a,b,c])
# gauss_contiguity_by_GKZ([b,a,a],[b,a+1,a+2])  # fails
#f=representative([[1,0,1],[0,1,2]],[2,1,1],gauss_contiguity,vars=[a,b,c]) # 多分isomでない.
#print(linear_poly_coefficients(2*a+3*b+4*c+5,[a,b,c]))    
#f=representative([[1,1,0],[0,1,1]],[2,1,1],gauss_contiguity,vars=[a,b,c]) 

##Ref: 2025-07-24-my-note-lattice.goodnotes
#f=representative([[1,0,0],[0,1,0],[0,0,1]],[0,0,0],gauss_contiguity,vars=[a,b,c]) 
print('Example: Try\n  f=representative([[1,0,0],[0,1,0],[0,0,1]],[0,0,0],gauss_contiguity,vars=[a,b,c])\nfor Gauss hg classification.') 
print('         See also Example: Appell F1 and Example: Appell F2')

r"""
# Example: 1F1
ring=PolynomialRing(QQ,'x,dx,a,c')
x,dx,a,c=ring.gens()
dic={'x':x,'dx':dx,'a':a,'c':c}
f=representative([[1,0],[0,1]],[0,0],c1f1_contiguity,vars=[a,c])
"""

r"""
# Example: Appell F1
ring=PolynomialRing(QQ,'x,dx,y,dy,a,b,bp,c')
x,dx,y,dy,a,b,bp,c=ring.gens()
dic={'x':x,'dx':dx,'y':y,'dy':dy,'a':a,'b':b,'bp':bp,'c':c}
f=representative([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],[0,0,0,0],f1_contiguity,vars=[a,b,bp,c])
"""

r"""
BUG:  2025-07-24 夜.  linear_poly_coefficients()  --> fixed 2025-07-25
 [latticce BT,origin S,param]= [[[1, 0, 0, 0], [0, -1, 1, 0], [0, 1, 0, 1]], (0, 0, 0, 1), [a, -b + bp, b, bp + 1]]  と bf= [(-b + bp, 1), (-b + bp - 1, 1), (b, 1), (-a + bp, 1), (-a + bp + 1, 1), (a, 1)] で代表点 Rep を探す. Conti に代表点 Rep を加える (todo)

Error: f= -b + bp , h= -b , [vars,new_vars]= [[a, b, bp, c], [a, -b + bp, b, bp + 1]]
Traceback (most recent call last):
  File "./2025-07-23-representative.py", line 66, in linear_poly_coefficients
    clist[m]=int(f-h)
"""


r"""
# Example: Appell F2
ring=PolynomialRing(QQ,'x,dx,y,dy,a,b,bp,c,cp')
x,dx,y,dy,a,b,bp,c,cp=ring.gens()
dic={'x':x,'dx':dx,'y':y,'dy':dy,'a':a,'b':b,'bp':bp,'c':c,'cp':cp}
f=representative([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]],[0,0,0,0,0],f2_contiguity,vars=[a,b,bp,c,cp])

#f=representative([[0,0,1,0,0],[0,0,0,0,1]],[0,0,0,0,0],f2_contiguity,vars=[a,b,bp,c,cp])
"""
