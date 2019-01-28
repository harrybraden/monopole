from PIL import Image
import data
from skimage import measure
import matplotlib.pyplot as plt
from numpy import sqrt, int, arange, rot90

def draw_slice(layer, ax, z, k):
    #levels = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5]
    contours = ax.contour(layer)#, levels=levels)
    ax.clabel(contours, inline=1, fontsize=9, inline_spacing=2)#, levels=levels)

    #for n, contour in enumerate(contours):
    #    ax.plot(contour[:, 1], contour[:, 0], linewidth=2)

    ax.imshow(layer, interpolation='nearest', alpha=0.5)

    ax.axis('image')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("z=%.02f, k=%.02f" % (z, k), y=0.9)

def draw_contours():
    print "# Trace contours"

    rows = 5
    cols = 5

    plt.rcParams["figure.figsize"] = (cols * 5, rows * 5)

    fig, axes = plt.subplots(ncols=cols, nrows=rows)
    fig.subplots_adjust(hspace=0, wspace=0)

    for c in range(0, cols):
        krange = arange(0.01, 0.99, 0.01)
        k = float(c) / float(cols) * 0.99 + 0.01

        print "# Load data", k 
        volume = data.load_data(k, smoothness=2) 
        volume = rot90(volume, 1, [0, 1])

        lim = len(volume) / 2
        for r in range(0, rows):
            zind = int(float(r) / float(rows) * lim) # lim - int(pow(float(r) / float(rows), 2) * lim)
            zrange = arange(0.025, 2.975, 0.05)
            z = float(zind + 1) / float(len(volume)) * 2.975 - 0.025
            print "LAYER", r, zind, len(volume), len(axes)
            layer = volume[zind + 1]
            draw_slice(layer, axes[r, c], z, k)

    plt.savefig('contours.png', dpi=100, bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    draw_contours()
