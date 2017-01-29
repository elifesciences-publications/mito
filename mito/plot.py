# -*- coding: utf-8 -*-

# Utilities for plotting
import matplotlib.pyplot as plt
import numpy as np

from functools import reduce
from skimage.color import gray2rgb, rgb2hsv, hsv2rgb
from skimage.exposure import rescale_intensity


def colourise(img, hue):
    """Colourise a grayscale images according to a hue.

    Inspired by http://nbviewer.jupyter.org/gist/jeanpat/9665267
    """
    hsv = rgb2hsv(gray2rgb(img))
    hsv[:, :, 0] = hue
    hsv[:, :, 1] = 1  # Full saturation
    return hsv2rgb(hsv)


def rgbstack(img, hues, rescale=False):
    if img.shape[0] != len(hues):
        raise ValueError("Number of images in stack must be "
                         "equal to number of hues")

    out = (colourise(i, h) for i, h in zip(img, hues))
    out = np.array(reduce(lambda x, y: x + y/3, out,
                          np.zeros(img[0].shape + (3, ),
                                   dtype=img.dtype)))

    if rescale:
        out = rescale_intensity(out,
                                in_range=(out.min(),
                                          out.max()-0.15),
                                out_range=(0, 255))

    return np.uint8(out)


def output_r90_img(img, df, cal, output_path, hues=[0.6, 0, 0.4, 1]):
    colour_img = rgbstack(img, hues, rescale=True)

    fig, ax0 = plt.subplots()
    ax0.imshow(colour_img)
    for n, g in df.groupby("Infected"):
        ax0.scatter(g["xc"].tolist(), g["yc"].tolist(), s=50,
                    c="r" if n else "b",
                    marker="*" if n else "o")
        for r in g.itertuples():
            circ = plt.Circle((r.xc, r.yc), r._5/cal,
                              facecolor="none",
                              edgecolor="red")
            ax0.add_artist(circ)

    ax0.axis("off")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def vertical_line(x, **kwargs):
    plt.axvline(0, **kwargs)


# Polar histogram
def polar_hist(a, bins=10, hist_kws=None, density=False, axlabel=None,
               color=None, label=None, ax=None):
    if ax is None:
        ax = plt.gca()

    # Intelligently label the support axis (from seaborn distplot)
    label_ax = bool(axlabel)
    if axlabel is None and hasattr(a, "name"):
        axlabel = a.name
        if axlabel is not None:
            label_ax = True

    a = np.rad2deg(np.asarray(a).squeeze())

    # Handle dictionary defaults
    if hist_kws is None:
        hist_kws = dict()

    # Get the color from the current color cycle
    if color is None:
        line, = ax.plot(a.mean(), 0)
        color = line.get_color()
        line.remove()

    # Plug the label into the right kwarg dictionary
    if label is not None:
        if hist_kws:
            hist_kws["label"] = label

    h, b = np.histogram(a, bins=bins, density=density)
    centre = (b[:-1] + b[1:]) / 2
    width = np.pi*2/bins

    hist_kws.setdefault("alpha", 0.4)
    hist_color = hist_kws.pop("color", color)
    ax.bar(np.deg2rad(centre), h, width=width, bottom=0.0,
           color=hist_color, **hist_kws)
    if hist_color != color:
        hist_kws["color"] = hist_color

    if label_ax:
        ax.set_xlabel(axlabel)

    return ax
