# GRSA
一般化範囲交換を用いた画像処理プログラム 
### 仕様
一般化範囲交換アルゴリズムを用いてエネルギーの最小化を試みます.  
- D<sub>p</sub>(f<sub>p</sub>) = |I<sub>p</sub> - f<sub>p</sub>|
- V<sub>pq</sub>(f<sub>p</sub>, f<sub>q</sub>) = |f<sub>p</sub> - f<sub>q</sub>|
#### 記号とか
- k:ラベルの最大値   
- i = 0: source  
- i = k + 1: sink  
#### 各枝の重み設定
そのうちちゃんと書く
- p<sub>i</sub>とq<sub>j</sub>を結ぶ枝  
    - w(e<sup>pq</sup><sub>ij</sub>) =(| i − j + 1 | − 2 | i − j | + | i − j − 1 |) / 2.0  
- p<sub>i</sub>とp<sub>i+1</sub>を結ぶ枝  
    - w(e<sup>p</sup><sub>i</sub>) = |I<sub>p</sub> - i| - ∑<sub>(p, q) ∈ E</sub> (-(|k + 1 - i| + |i + 1|) / 2.0)  
- p<sub>i+1</sub>とp<sub>i</sub>を結ぶ枝  
    - w(e<sup>p</sup><sub>i-</sub>) = ∞  
---
### 不具合
~~V<sub>pq</sub>(f<sub>p</sub>, f<sub>q</sub>) = (f<sub>p</sub> - f<sub>q</sub>)<sup>2</sup>と設定した際,期待する結果が出力されず~~  
~~V<sub>pq</sub>(f<sub>p</sub>, f<sub>q</sub>)の条件は凸関数より,どこかにバグ🐜~~
### 更新時メモ
- 11/8 基にするプログラムの大きなバグ修正に伴い,masterブランチ初期化.
---
### 参考  
[石川 博 グラフカット 2007 年 3 月 CVIM 研究会チュートリアル／情報処理学会研究報告 2007-CVIM-158-(26) pp. 193-204.](http://www.vision.cs.chubu.ac.jp/CV-R/pdf/HiroshiCVIM2007.pdf)
