load("2025-07-23-representative.py")
load("2025-08-06-to_bounded_nk.py")

"""
以下 2025-08-28-show_tree_nk.py より (中山プログラム)
l はrepresentative の返る値 
e.g. 
l = load("sb-F2-representative.sobj") 

l[0] : 関数情報
l[1] : 隣接関係式の情報 
    l[1][1] : bf 
l[2] 以降 l[1] の結果に対応して，制限した部分に対する計算結果(l[0],l[1],l[2] と同様の構造を持つ) 

["0-dim face", 点] が木の節点 
"""
def show_tree(l, depth,prev=None):
    n = len(l)
    if l[0]=="0-dim face":
        print("/*&C \n");
        print("    "*depth, end="")
        print(l)
        print("*/\n")
#        T=matrix(prev[0][1][0]).transpose();
#        offset=vector(prev[0][1][1]);
#        print('T=',T,', offset=',offset)
#        print('prev[1][1][0] =',prev[1][1][0])
        return
    print("    "*depth, end="")
    print(l[0][1][0])
    print("    "*depth, end="")
    print("bf", l[1][1])  # bf 

    # 代表点の計算
    bf = l[1][1]
    h = bf2ha(bf, bf_vars(bf))
    pts = rep_pts(h)
# By NT
#    print("    "*depth, end="")
#    print("pts", pts); 
    T=matrix(l[0][1][0]).transpose();
    offset=vector(l[0][1][1]);
    print('/*&C \n\ncontiguity---');
    print(l[1][0]);
    print('---contiguity \n*/');
    print('/*&C T.transpose()=',l[0][1][0],' offset=', l[0][1][1]);
    print(' bf = ',l[1][1],' \n*/')
    print('matrix T=',T);
    print('offset =',offset);
    for i in range(len(pts)):
        for j in range(len(pts[i])):
            if pts[i][j][1] != None:
                new_abc=pts[i][j][1];
                print('new_abc=',new_abc)
                print('/*&C ',new_abc,' --> ', T*vector(new_abc)+offset,' \n*/')
## end of by NT

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
        show_tree(l[i], depth+1,prev=l)

# l は多項式を因数分解した結果のリスト --> 変数リストを返す
# e.g. [(-a + c - 1, 1), (-a + c, 1), (-b + c - 1, 1), (-b + c, 1), (b, 1), (a, 1)]
def bf_vars(l):
    s = set()
    for t in l:
        s = s.union(t[0].variables())
    return list(s)

