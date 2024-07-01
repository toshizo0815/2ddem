import numpy as np
import matplotlib.pyplot as plt
import math


#===================■ フォントやフォントサイズの設定 ■=====================
plt.rcParams["text.usetex"] = True  #Tex文字使用の有無（True: Tex文字使用）
plt.rcParams["font.size"] = 15              # 基本となるフォントの大きさ
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
d = 1
theta = 6
d_str = str(d).zfill(1)
theta_str = str(theta).zfill(2)
dg = math.sqrt(d/(9.8*1000))

state=d_str+'_'+theta_str+'_'    #---粒子径と斜面角度---
x, y0=np.loadtxt('xt_'+state+'000', unpack=True)
x, y1=np.loadtxt('xt_'+state+'016', unpack=True)
x, y2=np.loadtxt('xt_'+state+'032', unpack=True)
x, y3=np.loadtxt('xt_'+state+'048', unpack=True)
x, y4=np.loadtxt('xt_'+state+'064', unpack=True)
x, y5=np.loadtxt('xt_'+state+'080', unpack=True)

xn = x / dg
#---x^1とx^2のライン---
# ypow = 0.003*xn**2    #---10度---
# ylin = 0.004*xn

# ypow = 0.002*xn**2    #---8度---
# ylin = 0.004*xn

ypow = 0.001*xn**2    #---6度---
ylin = 0.01*xn

# ypow = 2*xn**2    #---5度---
# ylin = 0.15*xn

# ypow = 0.001*xn**2    #---4度---
# ylin = 0.001*xn

# ypow = 0.0005*xn**2    #---2度---
# ylin = 0.00003*xn

#---グラフの作成---
fig = plt.figure(figsize=(7.15,5), dpi=100) # figsize=(px, py)：px, pyはグラフの縦と横の大きさ
ax=fig.add_subplot(111) #グラフ描画領域の生成（おまじない）

#---プロット（linestyle：線の種類 / linewidth：線の太さ / marker：マーカーの種類etc）---
ax.plot(xn[0], y0[0], color="#000000", linestyle='solid', linewidth=2, marker='o', markersize=3, label="$\phi=0.00$")    #ラベルを一回だけ出力する
ax.plot(xn[0], y1[0], color="#FF4B00", linestyle='solid', linewidth=2, marker='^', markersize=3, label="$\phi=0.16$")
ax.plot(xn[0], y2[0], color="#005AFF", linestyle='solid', linewidth=2, marker='v', markersize=3, label="$\phi=0.32$")
ax.plot(xn[0], y3[0], color="#03AF7A", linestyle='solid', linewidth=2, marker='p', markersize=3, label="$\phi=0.48$")
ax.plot(xn[0], y4[0], color="#4DC4FF", linestyle='solid', linewidth=2, marker='D', markersize=3, label="$\phi=0.64$")
ax.plot(xn[0], y5[0], color="#F6AA00", linestyle='solid', linewidth=2, marker='x', markersize=3, label="$\phi=0.80$")

for i in range(1, 3900, 15):
    ax.plot(xn[i], y0[i], color="#000000", linestyle='solid', linewidth=2, marker='o', markersize=3)
    ax.plot(xn[i], y1[i], color="#FF4B00", linestyle='solid', linewidth=2, marker='^', markersize=3)
    ax.plot(xn[i], y2[i], color="#005AFF", linestyle='solid', linewidth=2, marker='v', markersize=3)
    ax.plot(xn[i], y3[i], color="#03AF7A", linestyle='solid', linewidth=2, marker='p', markersize=3)
    ax.plot(xn[i], y4[i], color="#4DC4FF", linestyle='solid', linewidth=2, marker='D', markersize=3)
    ax.plot(xn[i], y5[i], color="#F6AA00", linestyle='solid', linewidth=2, marker='x', markersize=3)

ax.plot(xn, ypow, color="lightgray", linestyle='dashed', linewidth=2, marker='o', markersize=0)    
ax.plot(xn, ylin, color="lightgray", linestyle='dashed', linewidth=2, marker='o', markersize=0)


ax.set_yscale('log')  # y軸をlogスケールで描く
ax.set_xscale('log')  # x軸をlogスケールで描く
ax.set_xlim(10, 1000)     #---x軸の描画範囲---
ax.set_ylim(0.01, 1000)   #---y軸の描画範囲---
ax.set_xlabel(r"$t/\sqrt{d/g}$")           #---x軸のラベル（TeXを用いる場合はrをつける）---
ax.set_ylabel(r"$\Delta x/(\pi D)$")           #---y軸のラベル（TeXを用いる場合はrをつける）---
plt.grid(which='major',color='lightgray',linestyle=':')
plt.grid(which='minor',color='lightgray',linestyle=':')
ax.legend()
fig.tight_layout()  #---おまじない---
plt.show()          #---グラフを描画---

fig.savefig(path_out+"graph_"+state+".png")