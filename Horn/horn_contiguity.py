# by gemini, 2026.03.07
from sage.misc.sage_eval import sage_eval
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR

def load_horn_contiguity():
    global Asir_running
    try:
        print(Asir_running)
    except NameError:
        # SageMathの組み込みasirを利用するか、外部のasir.pyを利用する設定
        load('asir.py') 
        Asir_running = 'asir is already running'
    
    # Risa/Asir側の最新ファイルをロード
    asir.eval('load("2026-03-01-horn-by-gkz.rr");')

def set_horn_A(A_matrix):
    """
    行列AをAsirにセットし、対応する変数情報を取得して
    Sage側のPolynomialRingと変換用辞書を自動生成する。
    """
    asir.eval(f'set_A({A_matrix});')
    
    # AsirからHorn変数とパラメータのリスト(文字列)を取得
    x_str = asir.eval('H_horn_xvar').strip('[]')
    a_str = asir.eval('H_a').strip('[]')
    
    x_vars = [v.strip() for v in x_str.split(',') if v.strip()]
    a_vars = [v.strip() for v in a_str.split(',') if v.strip()]
    
    # xに対応する微分作用素 dx を生成 (Asirの poly_dvar 規則に準拠)
    dx_vars = ['d' + v for v in x_vars]
    
    # Sage側のリングを作成
    all_vars = x_vars + dx_vars + a_vars
    var_names = ','.join(all_vars)
    
    R = PolynomialRing(QQ, var_names)
    gens = R.gens()
    dic = dict(zip(all_vars, gens))
    
    print(f"[*] Asir に行列 A をセットしました。")
    print(f"[*] 自動生成されたリング変数: {var_names}")
    
    return R, dic

def horn_contiguity_by_GKZ(old, new, ring, dic):
    """
    Hornパラメータのシフトからcontiguityを計算し、Sageのオブジェクトとして返す。
    old, new: [a1, a2, ...] のようなリスト形式
    """
    oo = str(old).replace("'", "")
    nn = str(new).replace("'", "")
    
    # Asir側で計算し、一時変数 Ans_tmp に格納する
    # ※直接結果文字列を受け取るとカンマ分割が難しいため
    asir.eval(f'Ans_tmp = horn_contiguity_from_a_shift({oo}, {nn});')
    
    # 各要素を個別に文字列として取得
    c0_str = asir.eval('Ans_tmp[0]')
    b0_str = asir.eval('Ans_tmp[1]')
    info_str = asir.eval('Ans_tmp[2]')
    uv_str = asir.eval('Ans_tmp[3]')
    
    # 要素数を確認して Common_factor を取得
    len_ans = int(asir.eval('length(Ans_tmp)'))
    if len_ans >= 5:
        cf_str = asir.eval('Ans_tmp[4]')
    else:
        cf_str = '1'
    
    # Sage の多項式/有理式に変換
    C0 = sage_eval(c0_str, locals=dic)
    B0 = sage_eval(b0_str, locals=dic)
    CF = sage_eval(cf_str, locals=dic) if cf_str != '1' else ring(1)
    
    # 情報リストはそのまま文字列として保持
    return [C0, B0, info_str, uv_str, CF, ring]


def get_horn_system(ring, dic):
    """
    Asir 側から Horn 型超幾何系の定義方程式 (微分作用素) を取得し、
    Sage の Symbolic Expression (SR) として未展開のまま返す。
    """
    print("[*] Asir から Horn system の定義方程式を取得します (未展開)...")
    
    # H_horn_str2 の要素数を取得
    length_str = asir.eval('length(H_horn_str2)')
    try:
        length = int(length_str)
    except ValueError:
        print(f"Error: 予期しない長さのデータです ({length_str})")
        return []
    
    # 現在の辞書 (dic) の変数名を使って、SR用の変数辞書を自動生成
    # 例: 'a1' -> SR.var('a1')
    sr_dic = {var_name: SR.var(var_name) for var_name in dic.keys()}
    
    horn_ops = []
    for i in range(length):
        # Asir から i 番目の作用素の文字列を取得
        op_str = asir.eval(f'H_horn_str2[{i}]')
        
        # Sage の Symbolic Ring (SR) オブジェクトに変換
        try:
            op_sage = sage_eval(op_str, locals=sr_dic)
            horn_ops.append(op_sage)
        except Exception as e:
            print(f"文字列のパースに失敗しました: {op_str}")
            print(f"エラー詳細: {e}")
            horn_ops.append(op_str)
            
    return horn_ops


# 1. モジュールのロード
load_horn_contiguity()

# 2. 行列Aのセットと、Sage側リングの自動構築 (Gauss の例)
A_F1 = [[1,0,0,-1], [0,1,0,1], [0,0,1,1]]
R, dic = set_horn_A(A_F1)

# ここで R は x5, x6, dx5, dx6, a1, a2 などの変数を持つリングになります
# a1, a2 は dic['a1'] などを介してアクセスできます
a1, a2, a3 = dic['a1'], dic['a2'], dic['a3']

# 3. Contiguity の計算 (例: a2 -> a2 + 1)
L = horn_contiguity_by_GKZ([a1, a2, a3], [a1, a2+1, a3], ring=R, dic=dic)

# 結果の表示
print("分子 (Numerator):", L[0])
print("分母 (Denominator):", L[1])
print("共通因子 (Common Factor):", L[4])

# ----------------------
# 1. Appell F1 の行列をセットして、ベースとなるリングを取得
A_F1 = [[1,0,0,-1], [0,1,0,1], [0,0,1,1]]
R_base, dic_base = set_horn_A(A_F1)

# 2. 退化パラメータ用の変数 p, q を追加した「拡張リング」を作成する
base_vars = list(R_base.variable_names())
ext_vars = base_vars + ['p', 'q']  # p, q を追加
R_ext = PolynomialRing(QQ, ext_vars)

# 新しいリングの変数で辞書 (dic) を作り直す
dic_ext = dict(zip(ext_vars, R_ext.gens()))

# パラメータ p, q を取得
p = dic_ext['p']
q = dic_ext['q']

# 3. 退化した状態でのシフトを定義
# 例: a1 = p-1, a2 = p, a3 = q から出発し、p を p+1 にシフトする
A_start = [p - 1, p, q]
A_end   = [p, p + 1, q]

print(f"\n[*] 退化パラメータのシフトを計算します:")
print(f"    Start : {A_start}")
print(f"    End   : {A_end}")

# 4. 拡張したリングと辞書を渡して contiguity を計算
L_deg = horn_contiguity_by_GKZ(A_start, A_end, ring=R_ext, dic=dic_ext)

print("\n=== Asir からの計算結果 ===")
print("分子 (Numerator):")
print(L_deg[0])
print("\n分母 (Denominator):")
print(L_deg[1])
print("\n共通因子 (Common Factor):")
print(L_deg[4])


# 2. Horn system の取得と表示
horn_system = get_horn_system(R_base, dic_base)

print("\n=== Horn System (定義方程式) ===")
for i, op in enumerate(horn_system):
    print(f"P{i+1} = {op}")
