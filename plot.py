import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("macVelOut.txt")

xi  = 1
yi  = 1
vxi = 1
vyi = 1

# 5×5サイズのFigureを作成してAxesを追加
fig = plt.figure(figsize = (10, 10))
ax = fig.add_subplot(111)

# 格子点を表示
ax.grid()

# 軸ラベルの設定
ax.set_xlabel("x(m)", fontsize = 16)
ax.set_ylabel("y(m)", fontsize = 16)
plt.xticks( np.arange(0.0, 1.1, 0.1))
plt.yticks( np.arange(0.0, 1.1, 0.1))

# 軸範囲の設定
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)

# x軸とy軸
ax.axhline(0, color = "gray")
ax.axvline(0, color = "gray")

# ベクトルを表示
# quiver(始点x,始点y,成分x,成分y)
ax.quiver(data[:,0], data[:,1], data[:,2], data[:,3], color = "red",
          angles = 'xy', scale_units = 'xy', scale = 8)

plt.savefig("macOutput.png")
plt.show()
