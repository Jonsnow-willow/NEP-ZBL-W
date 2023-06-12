import numpy as np
import matplotlib.pyplot as plt

# 设置图片大小和图片风格
plt.rcParams["figure.figsize"] = (6, 6)
plt.rcParams.update({"font.size": 10, "text.usetex": False})

# 如果是测试数据，设置为 True；否则设置为 False
test = True

# 加载数据
loss = np.loadtxt("loss.out")
if test:
    energy = np.loadtxt("energy_test.out")
    virial = np.loadtxt("virial_test.out")
    force = np.loadtxt("force_test.out")
    loss = np.hstack((loss[:, 2:4], loss[:, 7:10]))
else:
    energy = np.loadtxt("energy_train.out")
    virial = np.loadtxt("virial_train.out")
    force = np.loadtxt("force_train.out")
    loss = loss[:, 2:7]

for i in range(virial.shape[0]):
    if virial[i, 1] == 0:
        virial[i, 0] = 0

# 绘制图形
fig, axes = plt.subplots(2, 2)

# subplot (a)
axes[0, 0].loglog(loss)
axes[0, 0].set_title("(a)")
axes[0, 0].set_xlabel("Generation")
axes[0, 0].set_ylabel("Loss functions")
axes[0, 0].legend(["L1", "L2", "Energy", "Force", "Virial"], loc="lower left")

# subplot (b)
axes[0, 1].plot(energy[:, 1], energy[:, 0], ".", markersize=10, color=[0, 0.45, 0.74])
axes[0, 1].plot(axes[0, 1].get_xlim(), axes[0, 1].get_xlim(), linewidth=2, color=[0.85, 0.33, 0.1])
axes[0, 1].set_title("(b)")
axes[0, 1].set_xlabel("DFT energy (eV/atom)")
axes[0, 1].set_ylabel("NEP energy (eV/atom)")
axes[0, 1].axis([-8.8, -4.2, -8.8, -4.2])

# subplot (c)
axes[1, 0].plot(force[:, 3:6], force[:, 0:3], ".", markersize=10, color=[0, 0.45, 0.74])
axes[1, 0].plot(axes[1, 0].get_xlim(), axes[1, 0].get_xlim(), linewidth=2, color=[0.85, 0.33, 0.1])
axes[1, 0].set_title("(c)")
axes[1, 0].set_xlabel("DFT force (eV/Å)")
axes[1, 0].set_ylabel("NEP force (eV/Å)")
axes[1, 0].axis([-30, 30, -30, 30])

# subplot (d)
axes[1, 1].plot(virial[:, 1], virial[:, 0], ".", markersize=10, color=[0, 0.45, 0.74])
axes[1, 1].plot(axes[1, 1].get_xlim(), axes[1, 1].get_xlim(), linewidth=2, color=[0.85, 0.33, 0.1])
axes[1, 1].set_title("(d)")
axes[1, 1].set_xlabel("DFT virial (eV/atom)")
axes[1, 1].set_ylabel("NEP virial (eV/atom)")
axes[1, 1].axis([-6, 11, -6, 11])

plt.subplots_adjust(hspace=0.4, wspace=0.3)

# 保存图片为 PDF 或其他格式
# plt.savefig("output.pdf", bbox_inches="tight")
plt.show()