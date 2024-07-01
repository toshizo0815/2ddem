import numpy as np
import matplotlib.pyplot as plt
import math


#===================■ フォントやフォントサイズの設定 ■=====================
plt.rcParams["text.usetex"] = True  #Tex文字使用の有無（True: Tex文字使用）
plt.rcParams["font.size"] = 18              # 基本となるフォントの大きさ
#===================■ おまじない ■=========================================
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams["mathtext.cal"] = "serif"      # TeX表記に関するフォント設定
plt.rcParams["mathtext.rm"] = "serif"       # TeX表記に関するフォント設定
plt.rcParams["mathtext.it"] = "serif:italic"# TeX表記に関するフォント設定
plt.rcParams["mathtext.bf"] = "serif:bold"  # TeX表記に関するフォント設定
plt.rcParams["mathtext.fontset"] = "cm"     # TeX表記に関するフォント設定
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
#===================■ 凡例に関する設定（好み）==============================
plt.rcParams["legend.loc"] = "best"
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.facecolor"] = "None"
plt.rcParams['legend.handlelength'] = 0.5
plt.rcParams['legend.numpoints'] = 1
#===================■ グリッドに関する設定 ■===============================
plt.rcParams["grid.color"] = "gray"        # グリッドの色
plt.rcParams["grid.linewidth"] = 1.0        # グリッドの線幅

#---パスの設定---
path_data="/home/terai/Documents/gotolab/dem_2d/"  #プロットしたいデータのディレクトリパス
path_out="/home/terai/Documents/gotolab/dem_2d/"  #グラフを保存するディレクトリパス

#---データの読み込み---
state='6_10'
x, y, yerr=np.loadtxt('xt_err_'+state+'_', unpack=True)
# y1, y2, y3=np.loadtxt('xt_cylinder_all.d', unpack=True)


#---グラフの作成---
fig = plt.figure(figsize=(7.15,5), dpi=100) # figsize=(px, py)：px, pyはグラフの縦と横の大きさ
ax=fig.add_subplot(111) #グラフ描画領域の生成（おまじない）

#---プロット（linestyle：線の種類 / linewidth：線の太さ / marker：マーカーの種類etc）---


ax.plot(x, y, color="red", linestyle='solid', linewidth=2, marker='o', markersize=0)

# ax.plot(x, y1, color="blue", linestyle='solid', linewidth=2, marker='o', markersize=0)
# ax.plot(x, y2, color="blue", linestyle='solid', linewidth=2, marker='o', markersize=0)
# ax.plot(x, y3, color="blue", linestyle='solid', linewidth=2, marker='o', markersize=0)

# ax.set_yscale('log')  # y軸をlogスケールで描く
# ax.set_xscale('log')  # x軸をlogスケールで描く
# ax.set_xlim(0, 200)     #---x軸の描画範囲---
# ax.set_ylim(-1, 1)   #---y軸の描画範囲---
ax.set_xlabel(r"$t/\sqrt{d/g}$")           #---x軸のラベル（TeXを用いる場合はrをつける）---
ax.set_ylabel(r"$\Delta x/(\pi D)$")           #---y軸のラベル（TeXを用いる場合はrをつける）---
plt.errorbar(x, y, yerr, capsize=3, fmt='o', markersize=5, ecolor='red', markeredgecolor = "red", color='w')
plt.grid(which='major',color='lightgray',linestyle=':')
plt.grid(which='minor',color='lightgray',linestyle=':')
ax.legend()
fig.tight_layout()  #---おまじない---
plt.show()          #---グラフを描画---

fig.savefig(path_out+"graph"+state+".png")   #---ファイル名+保存形式（.png, .svg, .pdf）---