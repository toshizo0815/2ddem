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
matplotlib.use("Agg")

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

d = 2
theta = 6
phi = 32
step = 300000

d_str = str(d).zfill(1)
theta_str = str(theta).zfill(2)
phi_str = str(phi).zfill(3)
step_str = str(step).zfill(7)
time_str = str(step/2000).zfill(5)

path1="/home/terai/Documents/gotolab/dem_2d/vector_data_"+d_str+'_'+theta_str+'_'+phi_str+'/'
path2="/home/terai/Documents/gotolab/dem_2d/"

x, y, u, v = np.loadtxt(path1+'v_vector_data_'+step_str, unpack=True)


fig, ax = plt.subplots()                        #よくわからんおまじない
ax.set_aspect('equal')                          #x軸とy軸の比率を1:1にする。これがないと横長で作成される。

# circle3 = patches.Circle(xy = (0,0), radius = 33, ec = "lightgray", fc = "none", lw = 20)
circle1 = patches.Circle(xy = (14.5,14.5), radius = 15, ec ='#000000', fc = "none", lw = 1)
# circle2 = patches.Circle(xy = (30,30), radius = 36, ec = "black", fc = "none", lw = 1)

# ax.add_patch(circle3)
ax.add_patch(circle1)
# ax.add_patch(circle2)

# circle0 = patches.Circle(xy = (x1 - xc1, y1 - yc1), radius = r, ec = "black", fc = "none", lw = 1)
# ax.add_patch(circle0)
ax.quiver(x, y, u, v, width = 0.004, scale = 0.6, color='#005AFF')

ax.set_xlim(-0.5, 29.5)     #---x軸の描画範囲---
ax.set_ylim(-0.5, 29.5)   #---y軸の描画範囲---
ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)

filename = 'vector_'+d_str+'_'+theta_str+'_'+phi_str+'_'+time_str+'.png' 
plt.xlabel(r"$x$")                                     #x軸の名前の設定
plt.ylabel(r"$y$")                                     #y軸の名前の設定
plt.savefig(filename, dpi=300)                           #画像を出力

