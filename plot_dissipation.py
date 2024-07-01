#coding: UTF-8
from curses.textpad import rectangle
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
import cv2
import matplotlib.cm as cm
import tqdm
import math
# matplotlib.use("Agg")

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

path1="/home/terai/Documents/gotolab/dem_2d/"
path2="/home/terai/Documents/gotolab/dem_2d/"

x, y, dis = np.loadtxt('dissipation_per_unittime_2.d', unpack=True)

range = 30
y=y[::range]
x=x[0:range]
dis = dis.reshape(range,range)

# fig, ax = plt.subplots()                        #よくわからんおまじない
# ax.set_aspect('equal')                          #x軸とy軸の比率を1:1にする。これがないと横長で作成される。



fig = plt.figure(figsize=(7.15,5), dpi=100) # figsize=(px, py)：px, pyはグラフの縦と横の大きさ
plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.15)
ax=fig.add_subplot(111) #グラフ描画領域の生成（おまじない）

ax.set_xticks([0, 5, 10, 15, 20, 25, 30])
ax.set_yticks([0, 5, 10, 15, 20, 25, 30])
# ax.set_xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
# ax.set_yticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
cax=ax.pcolormesh(x, y, dis, cmap='Blues')
pp = fig.colorbar(cax, ax=ax, orientation="vertical")
# pp.set_clim(0.00,1.00)
# fig.colorbar(cax)

circle1 = patches.Circle(xy = (14.5,14.5), radius = 15, ec ='#000000', fc = "none", lw = 1)
# circle1 = patches.Circle(xy = (30,30), radius = 30, ec ='#000000', fc = "none", lw = 1)
ax.add_patch(circle1)

# ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
# ax.axes.xaxis.set_visible(False)
# ax.axes.yaxis.set_visible(False)

filename = 'dissipation_t_9.png' 
plt.xlabel(r"$x$")                                     #x軸の名前の設定
plt.ylabel(r"$y$")                                     #y軸の名前の設定
plt.show()
plt.savefig(filename, dpi=300)                           #画像を出力

