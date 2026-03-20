#
# 2025-06-30-a2s.py の機能強化版.
# このプログラムを動かすには asir.py, 2025-04-01-contiguity.rr, 2025-07-17-F1-F2-by-gkz.rr が必要.
#
#Ref: load('asir.py') or from sage.interfaces.asir import load_reduce_Asir
#
# orange3n:  conda  activate sage
#            sage
#            load("2025-07-18-a2s.py")
#            quit
#            conda deactivate
#
#    conda -> mamba in ~/miniforge3/bin
#
def load_contiguity():
    global Asir_running
    try:
        print(Asir_running)
    except NameError:
        load('asir.py')
        Asir_running='asir is already running'
    asir.eval('load(\"2025-04-01-contiguity.rr\")')
    asir.eval('load(\"2025-07-17-F1-F2-by-gkz.rr\")')
    asir.eval('load(\"2025-07-28-f0134-by-gkz.rr\")')

def gauss_contiguity_by_GKZ(old,new,ring=None):
    if ring==None:
        ring=PolynomialRing(QQ,'x,dx,a,b,c')
    x,dx,a,b,c=ring.gens()
    dic={'x':x,'dx':dx,'a':a,'b':b,'c':c}
    oo=str(old)
    nn=str(new)
    L=asir.eval('gauss_contiguity_by_GKZ('+oo+','+nn+')')
    print(L)  # for debug
    L2=L.split(',') # これは文字列成分のリストにしてくれないので以下で修正.
    L3=[]
    info=''
    for i,s in enumerate(L2):
        if (i<2):
            s = s.replace('[','')
            f=sage_eval(s,locals=dic)
            L3.append(f)
        else:
            info = info+s+','
    L3.append(info[:-1])
    L3.append(ring)
    return L3

"""
##Ref: https://doc.sagemath.org/html/en/reference/misc/sage/misc/sage_eval.html 
poly_asir=asir.eval('(x+1)^2')   #type(poly_asir) is a string
x = PolynomialRing(RationalField(),"x").gen()
f = sage_eval(poly_asir,locals={'x':x})
# str(f)
"""

"""
# 開発メモ 1
load('2025-06-30-a2s.py')
load_contiguity()
L=gauss_contiguity_by_GKZ('[a,b,c]','[a-1,b,c]')
R.<x,dx,a,b,c>=QQ[]
f=sage_eval(str(L[0]),locals={'x':x,'dx':dx,'a':a,'b':b,'c':c})
f.parent()  # f の ring str を取り出す.
"""

# pattern match 版. bug いりかもしれない.
def c1f1_contiguity(old,new,ring):
    print('Warning: c1f1_contiguity, pattern match 版. bug いりかもしれない.')
    if ring==None:
        ring=PolynomialRing(QQ,'x,dx,a,c')
    x,dx,a,c=ring.gens()
    dic={'x':x,'dx':dx,'a':a,'c':c}
    oo=str(old)
    nn=str(new)
    L=asir.eval('c1f1_contiguity('+oo+','+nn+')')
    print(L)  # for debug
    L2=L.split(',') # これは文字列成分のリストにしてくれないので以下で修正.
    L3=[]
    info=''
    for i,s in enumerate(L2):
        if (i<2):
            s = s.replace('[','')
            f=sage_eval(s,locals=dic)
            L3.append(f)
        else:
            info = info+s+','
    L3.append(info[:-1])
    L3.append(ring)
    return L3

def f1_contiguity_by_GKZ(old,new,ring=None):
    if ring==None:
        ring=PolynomialRing(QQ,'x,dx,y,dy,a,b,bp,c')
    x,dx,y,dy,a,b,bp,c=ring.gens()
    dic={'x':x,'dx':dx,'y':y,'dy':dy,'a':a,'b':b,'bp':bp,'c':c}
    oo=str(old)
    nn=str(new)
    L=asir.eval('f1_contiguity_by_GKZ('+oo+','+nn+')')
    print(L)  # for debug
    L2=L.split(',') # これは文字列成分のリストにしてくれないので以下で修正.
    L3=[]
    info=''
    for i,s in enumerate(L2):
        if (i<2):
            s = s.replace('[','')
            f=sage_eval(s,locals=dic)
            L3.append(f)
        else:
            info = info+s+','
    L3.append(info[:-1])
    L3.append(ring)
    return L3

def f2_contiguity_by_GKZ(old,new,ring=None):
    if ring==None:
        ring=PolynomialRing(QQ,'x,dx,y,dy,a,b,bp,c,cp')
    x,dx,y,dy,a,b,bp,c,cp=ring.gens()
    dic={'x':x,'dx':dx,'y':y,'dy':dy,'a':a,'b':b,'bp':bp,'c':c,'cp':cp}
    oo=str(old)
    nn=str(new)
    L=asir.eval('f2_contiguity_by_GKZ('+oo+','+nn+')')
    print(L)  # for debug
    L2=L.split(',') # これは文字列成分のリストにしてくれないので以下で修正.
    L3=[]
    info=''
    for i,s in enumerate(L2):
        if (i<2):
            s = s.replace('[','')
            f=sage_eval(s,locals=dic)
            L3.append(f)
        else:
            info = info+s+','
    L3.append(info[:-1])
    L3.append(ring)
    return L3

"""
# 実行例 1
load_contiguity()  
R=PolynomialRing(QQ,'x,dx,a,b,c')
x,dx,a,b,c=R.gens()
L=gauss_contiguity_by_GKZ([a,b,c],[a-1,b,c],ring=R)
Bf=list(L[1].factor())
"""

"""
# 実行例 2
load_contiguity()  
R=PolynomialRing(QQ,'x,dx,y,dy,a,b,bp,c')
x,dx,y,dy,a,b,bp,c=R.gens()
L=f1_contiguity_by_GKZ([a,b,bp,c],[a+1,b,bp,c],ring=R)
Bf=list(L[1].factor())
"""

def linear_poly_coefficients(f,vars,param):
    ring=vars[0].parent()
    ff=str(f)
    vv=str(vars)
    pp=str(param)
    L=asir.eval('linear_poly_coefficients_for_sage('+ff+','+vv+','+pp+')')
    #print(L)
    return sage_eval(L)

# 2025-07-28
def f0134_contiguity_by_GKZ(old,new,ring=None):
    if ring==None:
        ring=PolynomialRing(QQ,'x,dx,y,dy,a,b')
    x,dx,y,dy,a,b=ring.gens()
    dic={'x':x,'dx':dx,'y':y,'dy':dy,'a':a,'b':b}
    oo=str(old)
    nn=str(new)
    L=asir.eval('f0134_contiguity_by_GKZ('+oo+','+nn+')')
    print(L)  # for debug
    L2=L.split(',') # これは文字列成分のリストにしてくれないので以下で修正.
    L3=[]
    info=''
    for i,s in enumerate(L2):
        if (i<2):
            s = s.replace('[','')
            f=sage_eval(s,locals=dic)
            L3.append(f)
        else:
            info = info+s+','
    L3.append(info[:-1])
    L3.append(ring)
    return L3

def xm_nox():
    asir.eval('Xm_noX=1')
def gb_verbose(level):
    asir.eval('dp_gr_print('+str(level)+')')

# create asir server.
load_contiguity()

