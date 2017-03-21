from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

"""
def multiply_b(a):
    print("Will compute", a, 'times', b)
    c = 0
    for i in range(0, a):
        c = c + emb.b
    return c"""

def test(a):
    print(a)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make data.
    X = np.arange(0, 0.01, 0.01/( len(a)))
    Y = np.arange(0, 0.01, 0.01/( len(a)))
    X, Y = np.meshgrid(X, Y)
    Z = np.asarray(a)

    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    plt.show()