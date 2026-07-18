

"""
2026.07.18

adachi_multi_formula_final.py
安達マルチ公式 研究プロジェクト - 統合プログラム

本ファイルには、本研究で得られた主要な2つのアルゴリズムを統合している。

1. 二分法最適化ルーチン (adachi_bsplit)
   従来の素朴な実装（G-Block二重ループ）と比較し、計算速度の向上（約3.2倍）を
   実現した最適化アルゴリズム。(No.41 研究ノート準拠)
   pi = sum_k (-1)^k f(k)/4^k, f(k)=2/(4k+1)+2/(4k+2)+1/(4k+3) を
   A=sum (-1)^k/((4k+1)4^k), B=sum (-1)^k/((4k+2)4^k), C=sum (-1)^k/((4k+3)4^k)
   の3系列に分解し（pi=2A+2B+C）、各系列を二分法で計算する。

2. 16進数BBP桁抽出ルーチン (bbp_digit_extraction)
   安達則男が独自発見した16進4項公式を用い、円周率の任意の桁へ直接アクセス
   する計算手法を実装。(No.42 研究ノート準拠)
   pi = sum_k (1/16^k) * (8/(16k+4) + 4/(16k+6) + 4/(16k+8) - 1/(16k+14))

研究の背景とスタンス：
本研究は、既存のチュドノフスキー・アルゴリズム等の数学的な最高効率手法を
凌駕することを目的とするものではない。あくまで独自のアプローチを通じて、
代数的操作が数値計算の振る舞いにどのような影響を与えるかを「誠実な記録」
として蓄積することを主眼としている。
"""

import sys
import math
from datetime import datetime

sys.set_int_max_str_digits(0)  # Python 3.11以降の桁数制限を解除


# =====================================================================
# 1. 二分法最適化ルーチン（adachi_bsplit）  --  No.41準拠
# =====================================================================

def _bs(a, b, c):
    """sum_{k=a}^{b-1} (-1)^k/((4k+c)*4^k) を二分法で計算し (Q,T) を返す。L=T/Q"""
    if b - a == 1:
        return (4 * a + c), 1
    m = (a + b) // 2
    Qam, Tam = _bs(a, m, c)
    Qmb, Tmb = _bs(m, b, c)
    pow4 = (-4) ** (m - a)  # 公比 r=-1/4（符号付き、交代級数）
    T = Tam * Qmb * pow4 + Tmb * Qam
    Q = Qam * Qmb * pow4
    return Q, T


def adachi_bsplit(digits):
    """安達マルチ公式（基本形）+ 二分法により、円周率をdigits桁計算する。

    戻り値: (pi_str, n_terms) — pi_strは小数点を含む文字列
    """
    n = int(digits / 0.602) + 5  # 1項あたり約0.602桁
    QA, TA = _bs(0, n, 1)
    QB, TB = _bs(0, n, 2)
    QC, TC = _bs(0, n, 3)

    # pi = 2*TA/QA + 2*TB/QB + TC/QC を1本の分数に統合し、最後に1回だけ整数除算
    num = 2 * TA * QB * QC + 2 * TB * QA * QC + TC * QA * QB
    den = QA * QB * QC

    prec = digits + 20
    scale = 10 ** prec
    scaled = (num * scale) // den
    s = str(abs(scaled))
    pi_str = s[0] + '.' + s[1:digits + 15]
    return pi_str, n


# =====================================================================
# 2. 16進数BBP桁抽出ルーチン（bbp_digit_extraction）  --  No.42準拠
# =====================================================================

def _S(j, n):
    """sum_{k=0}^inf 16^(n-k)/(16k+j) の小数部分（BBP桁抽出の核）"""
    s = 0.0
    for k in range(0, n + 1):
        denom = 16 * k + j
        exp = n - k
        term = pow(16, exp, denom) / denom  # 高速べき乗剰余
        s += term
        s -= int(s)
    k = n + 1
    while True:
        denom = 16 * k + j
        term = 16.0 ** (n - k) / denom
        if term < 1e-17:
            break
        s += term
        k += 1
    return s - math.floor(s)


def bbp_digit_extraction(n, num_digits=8):
    """円周率のn桁目（16進、0-indexed）から num_digits 桁分を、
    それ以前の桁を計算せずに直接求める（安達則男 独自発見の16進4項公式を使用）。
    """
    x = 8 * _S(4, n) + 4 * _S(6, n) + 4 * _S(8, n) - _S(14, n)
    x = x - math.floor(x)
    if x < 0:
        x += 1
    digits = []
    for _ in range(num_digits):
        x *= 16
        d = int(x)
        digits.append(d)
        x -= d
    return ''.join('0123456789abcdef'[d] for d in digits)


# =====================================================================
# メイン実行部
# =====================================================================

if __name__ == "__main__":
    print("=== 安達マルチ公式 統合プログラム ===")
    print("1: 二分法によるpi計算（10進、任意桁数）")
    print("2: BBP桁抽出（16進、任意位置へジャンプ）")
    mode = input("モードを選んでください (1 or 2) = ") or "1"

    if mode == "1":
        digits = int(input("何桁計算しますか？ (例: 10000) = ") or "1000")
        start_dt = datetime.now()
        result, n_terms = adachi_bsplit(digits)
        end_dt = datetime.now()
        print()
        print(f"開始時刻　: {start_dt.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"終了時刻　: {end_dt.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"処理時間　: {(end_dt - start_dt).total_seconds():.3f} 秒")
        print(f"計算桁数　: {digits} 桁")
        print(f"使用項数　: {n_terms} 項")
        print()
        print(f"先頭100桁 : {result[:100]}")
        print(f"末尾100桁 : ...{result[-100:]}")
    else:
        n = int(input("何桁目(16進)から求めますか？ (例: 1000000) = ") or "0")
        num_digits = int(input("何桁分表示しますか？ (例: 16) = ") or "16")
        start_dt = datetime.now()
        result = bbp_digit_extraction(n, num_digits)
        end_dt = datetime.now()
        print()
        print(f"開始時刻　: {start_dt.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"終了時刻　: {end_dt.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"処理時間　: {(end_dt - start_dt).total_seconds():.3f} 秒")
        print(f"抽出位置　: 小数点以下 {n} 桁目（16進）から")
        print()
        print(f"抽出結果　: {result}")
