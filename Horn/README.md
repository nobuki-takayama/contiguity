
# Horn型超幾何系 同型分類パッケージ マニュアル

本パッケージは、Risa/Asir の $D$-加群計算機能と SageMath の多面体・整数計画計算機能を連携させ、Horn型超幾何系のパラメータ空間（$\mathbb{Z}^m$ の格子空間）全体における同型分類を自動遂行するための関数群を提供します。

## 1. Risa/Asir コア関数 (GKZ系を経由した Contiguity 導出)

GKZ 超幾何系から制限アルゴリズムを用いて、Horn 系の contiguity relation（隣接関係式）を計算します。

### `set_A(A)`

* **概要:** 行列 $A$ に付随する GKZ 系と、その制限として得られる Horn 型超幾何方程式系を初期化し、必要な大域変数（`H_horn_str`, `H_a` など）を設定します。
* **引数:**
* `A`: 左側が単位行列となっている整数要素の行列（リストのリスト）。



### `horn_contiguity_from_a_shift(A_start, A_end)`

* **概要:** Horn 系のパラメータを `A_start` から `A_end` へシフトさせた際の contiguity relation（微分作用素）と超幾何 $b$-関数を計算します。
* **引数:**
* `A_start`: シフト前のパラメータのリスト（例: `[a, b, c]`）。退化パラメータ（例: `[p-1, p, q]`）も指定可能です。
* `A_end`: シフト後のパラメータのリスト（例: `[a+1, b, c]`）。


* **戻り値:** `[Numerator, Denominator, Info, UV, Common_factor]` の 5 要素リスト。
* `Numerator` / `Denominator`: 正規化された contiguity 作用素の分子と分母。
* `Info`: シフトの方向を示す文字列情報。
* `UV`: GKZ系での計算に用いられた最短シフトベクトル $u, v$ を用いた作用素。
* `Common_factor`: 可換 GCD によって抽出・約分された、余分なパラメータ因子（Torsion成分）。


* **特徴:** 内部で「最短シフト探索（ヒューリスティクス）」を行い、計算時間の爆発と冗長な因子の発生を防ぎます。

---

## 2. SageMath 連携インターフェース (`asir.sage` / `asir.py`)

Risa/Asir の計算エンジンを SageMath から透過的に呼び出し、数式オブジェクトとして安全に操作するためのラッパー関数群です。

### `set_horn_A(A_matrix)`

* **概要:** 行列 $A$ を Risa/Asir に送信して初期化し、SageMath 側の多項式環（`PolynomialRing`）と変数辞書を自動生成します。
* **引数:**
* `A_matrix`: 初期化する行列 $A$（例: `[[1,0,0,-1], [0,1,0,1], [0,0,1,1]]`）。


* **戻り値:** `(ring, dic)` のタプル。
* `ring`: 生成された多項式環（例: `x4, dx4, a1, a2, a3` をジェネレータに持つ）。
* `dic`: 変数名（文字列）からジェネレータへの変換辞書。



### `get_horn_system(ring, dic)`

* **概要:** 現在セットされている Horn 型超幾何方程式の定義方程式を取得します。
* **引数:** `set_horn_A` で取得した `ring` と `dic`。
* **戻り値:** SageMath の Symbolic Ring (`SR`) オブジェクトのリスト。
* **特徴:** 完全に展開された多項式ではなく、因数分解された構造（例: `(dx4*x4 + a2)*(dx4*x4 + a3) - ...`）を保ったまま保持されるため、パラメータの代入や解析が容易です。

### `horn_contiguity_by_GKZ(old, new, ring, dic)`

* **概要:** `horn_contiguity_from_a_shift` を SageMath から安全に呼び出し、結果を Sage の多項式/有理式オブジェクトとして返します。
* **引数:** `old` (シフト前リスト), `new` (シフト後リスト), `ring`, `dic`。
* **戻り値:** `[Numerator, Denominator, Info, UV, Common_factor, ring]` のリスト。各要素は Sage のオブジェクトとしてパース済みです。

---

## 3. 同型分類の自動遂行関数 (`representative.py`)

パラメータ空間の次元を再帰的に落としながら、超幾何系の同型分類（セル分解）を実行します。

### `representative(BT, S, conti_func, vars, dic)`

* **概要:** 与えられた格子空間上で contiguity relation を網羅的に計算し、$b$-関数の一次因子を「壁（超平面）」とみなして空間を分割、代表元（representative）を抽出します。
* **制限事項:** この関数は各セルの点がすべて格子の基底での動きでつながっているという仮定しています(つまりセルの点達の Markov 基底が格子の基底であることを仮定)。 この制限を除いたプログラムは開発中。
* **引数:**
* `BT`: 現在探索している格子の基底ベクトル（リストのリスト）。初期値は標準基底（単位行列）。
* `S`: 現在の格子の起点（シフト量ベクトル）。初期値はゼロベクトル。
* `conti_func`: Asirを呼び出して contiguity を計算する関数（例: `gauss_contiguity`, `f1_contiguity`）。
* `vars`: 探索対象となるパラメータ変数のリスト（例: `[a, b, c]`）。
* `dic`: パラメータ変数とその文字列表現の対応辞書。


* **戻り値:** 再帰的なリスト構造を持つ分類ログ。各次元（セル）における壁の方程式、計算された作用素、到達した次元（`0-dim face` など）の情報が格納されます。
* **アルゴリズムの流れ:**
1. 現在の格子空間（`BT` + `S`）でシフトを発生させ、`conti_func` を呼び出して $b$-関数を計算。
2. 往復の $b$-関数の最小公倍数を因数分解し、同型を壊す「壁（一次因子）」を抽出。
3. 抽出した「壁」と現在の格子の交わり（次元が 1 つ低い部分格子）を計算し、新たな基底 `BT_new` と起点 `T` を取得。
4. 新しい格子空間に対して `representative(BT_new, T, ...)` を再帰的に呼び出す。これを 0 次元（点）に到達するまで繰り返す。



---

## 4. 基本的な実行ワークフロー (SageMath)

```python
# 1. パッケージのロード
load("horn_contiguity.py")
load("2026-03-21-representative.py")

# 2. 行列 A のセットと変数リングの自動生成
A = [[1,0,0,-1], [0,1,0,1], [0,0,1,1]] # Gauss超幾何 2F1 の例
R, dic = set_horn_A(A)
Id=get_horn_system(R,dic)  
print(Id)
# parameter が方程式にどう入っているか表示. 
# 2F1 の通常のパラメータを a, b, c とすると a1=c-1, a2=a, a3=b であることに注意.

# 3. パラメータ変数の取得
a1, a2, a3 = dic['a1'], dic['a2'], dic['a3']


# 4. 同型分類の実行 (3次元の一般空間からスタート)
BT_init = [[1,0,0], [0,1,0], [0,0,1]] # 3x3 単位行列 (基底)
S_init = [0, 0, 0]                    # 原点
vars_list = [a1, a2, a3]

# 自動探索の開始
f = representative(BT_init, S_init, gauss_contiguity, vars=vars_list)

# 結果は f に階層的に保存される

# 5. contiguity を一つ求めたい場合の例.
 horn_contiguity_by_GKZ([a1,a2,a3],[a1+1,a2+1,a3+1],R,dic)
 horn_contiguity_by_GKZ([a1,a1,-1],[a1+1,a1+1,-1],R,dic)  

```

---
## 5. サンプル入力ファイル
* 2026-03-21-test1.rr,  Gauss 超幾何 2F1
* 2026-03-21-test2.rr,  1F1
* 2026-03-21-test3.rr,  Appell function F_1
