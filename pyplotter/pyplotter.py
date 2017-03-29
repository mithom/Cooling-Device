from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from math import sqrt

"""
def multiply_b(a):
    print("Will compute", a, 'times', b)
    c = 0
    for i in range(0, a):
        c = c + emb.b
    return c"""

def test(a,distr,third):
    fig = plt.figure()
    fig2 = plt.figure()
    fig3 = plt.figure()
    ax = fig.gca(projection='3d')
    ax2 = fig2.gca(projection='3d')
    ax3 = fig3.gca(projection='3d')
    # Make data.
    X = np.arange(0, 0.01, 0.01/( sqrt(len(a))))
    Y = np.arange(0, 0.01, 0.01/( sqrt(len(a))))
    X, Y = np.meshgrid(X, Y)
    Z = np.reshape(np.asarray(a),(int(sqrt(len(a))),int(sqrt(len(a)))))
    Z_distr = np.reshape(np.asarray(distr),(int(sqrt(len(distr))),int(sqrt(len(distr)))))
    Z_third = np.reshape(np.asarray(third),(int(sqrt(len(third))),int(sqrt(len(third)))))

    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    surf = ax2.plot_surface(X, Y, Z_distr, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    surf = ax3.plot_surface(X,Y,Z_third, cmap=cm.coolwarm,
                            lineWidth=0, antialiased=False)
    plt.show()
