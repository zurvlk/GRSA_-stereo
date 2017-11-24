# GRSA for stereo imaging
一般化範囲交換アルゴリズムを用いた画像処理プログラム 
### 仕様
一般化範囲交換アルゴリズムを用いて2つの入力画像から視差を求めます.  
- D<sub>p</sub>(f<sub>p</sub>) = ||left(p) - right(p - f<sub>p</sub>)||<sub>2</sub>
- V<sub>pq</sub>(f<sub>p</sub>, f<sub>q</sub>) = λ*(f<sub>p</sub> - f<sub>q</sub>)<sup>2</sup>

### 更新時メモ
- 11/24 大域解の導出可能に,範囲移動時に最適解を選択できない場合あり.
---
### 参考  
[Kangwei Liu et al.  GRMA: Generalized Range Move Algorithms for the Efficient Optimization of MRFs International Journal of Computer Vision February 2017, Volume 121, Issue 3, pp 365–390](https://link.springer.com/article/10.1007/s11263-016-0944-z)
