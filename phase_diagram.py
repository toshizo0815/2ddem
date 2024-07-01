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
d=1
d_str=str(d).zfill(1)


x = np.arange(start=0, stop=14, step=2)
y = np.arange(start=0, stop=0.96, step=0.16)

x1, y1=np.loadtxt(path_data+"phase_data_1_2", unpack=True)
x2, y2=np.loadtxt(path_data+"phase_data_2", unpack=True)

phase = np.genfromtxt('./phase_d='+d_str+'_rigid')

counter1 = 0
counter2 = 0


#---グラフの作成---
fig = plt.figure(figsize=(7.15,5), dpi=100) # figsize=(px, py)：px, pyはグラフの縦と横の大きさ
ax=fig.add_subplot(111) #グラフ描画領域の生成（おまじない）

#---プロット（linestyle：線の種類 / linewidth：線の太さ / marker：マーカーの種類etc）---
for i in range(7):
    for j in range(6):
        if phase[5-j,i] == 0:
            ax.plot(x[i], y[j], color="#FF4B00", marker='o', markersize=18, clip_on=False)
            if counter1 == 0:
                ax.plot(x[i], y[j], color="#FF4B00", marker='o', markersize=18, label='rigid', clip_on=False)
                counter1 = counter1 + 1
        if phase[5-j,i] == 1:
            ax.plot(x[i], y[j], color="#005AFF", marker='^', markersize=15, clip_on=False)
            if counter2 == 0:
                ax.plot(x[i], y[j], color="#005AFF", marker='^', markersize=15, label='non-rigid', clip_on=False)
                counter2 = counter2 + 1
ax.plot(x1, y1*0.8, color="black", linestyle='solid', linewidth=2, marker='o', markersize=0, label='theory')
# ax.plot(x2, y2, color="orange", linestyle='solid', linewidth=2, \
#     marker='o', markersize=0)



ax.set_xlim(0, 12)     #---x軸の描画範囲---
ax.set_ylim(0, 0.80)   #---y軸の描画範囲---
ax.set_xlabel(r"$\theta[^{\circ}]$")           #---x軸のラベル（TeXを用いる場合はrをつける）---
ax.set_ylabel(r"$\phi$")           #---y軸のラベル（TeXを用いる場合はrをつける）---
ax.tick_params(pad=10)
plt.xticks([0, 2, 4, 6, 8, 10, 12])
plt.yticks([0, 0.16, 0.32, 0.48, 0.64, 0.80])
plt.grid(which='major',color='lightgray',linestyle=':')
plt.grid(which='minor',color='lightgray',linestyle=':')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=18)
fig.tight_layout()  #---おまじない---
plt.show()          #---グラフを描画---
fig.savefig(path_out+"phase_diagram_"+d_str+".png")
