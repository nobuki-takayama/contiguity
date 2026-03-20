# horn_contiguity.py を用いて representative を計算する機能のテスト

# 1. パッケージのロード
load("horn_contiguity.py")
load("2026-03-21-representative.py")


# 2. 行列 A のセットと変数リングの自動生成
A = [[1,0,0,-1], [0,1,0,1], [0,0,1,1]] # Gauss超幾何の例
R, dic = set_horn_A(A)
Id=get_horn_system(R,dic)  # parameter が方程式にどう入っているか?
print(Id)

# 3. パラメータ変数の取得
a1, a2, a3 = dic['a1'], dic['a2'], dic['a3']

# 4. 同型分類の実行 (3次元の一般空間からスタート)
BT_init = [[1,0,0], [0,1,0], [0,0,1]] # 3x3 単位行列 (基底)
S_init = [0, 0, 0]                    # 原点
vars_list = [a1, a2, a3]

# 自動探索の開始
f = representative(BT_init, S_init, horn_contiguity_by_GKZ, vars=vars_list,dic=dic)

# 結果は f に階層的に保存される
