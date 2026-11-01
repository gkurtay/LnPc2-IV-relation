"""Common figure style of the article (palette, panel background, marker shading, image cards).

polish(fig) is applied to every figure just before it is saved: plot panels get a light background
with a white grid, isolated markers a soft shadow, and structure images a rounded card."""
import numpy as np
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch
from matplotlib.collections import PathCollection
from matplotlib.image import AxesImage

BG = '#f6f7fb'          # panel background
CARD_EDGE = '#dfe3ec'
MC = {'Lu': '#2f5fd0', 'Y': '#e3792f', 'La': '#1a9e77', 'Eu': '#5b5f7a'}
SHADOW = [pe.SimplePatchShadow(offset=(0.8, -0.8), shadow_rgbFace='#1b2a4a', alpha=0.18), pe.Normal()]
STATBOX = dict(boxstyle='round,pad=0.35,rounding_size=0.6', fc='white', ec='#d5d9e4', lw=0.6)


def chip(ax, text, color, x=0.03, y=0.97, size=7.0):
    """Coloured title chip in the upper left corner of a panel."""
    ax.text(x, y, text, transform=ax.transAxes, ha='left', va='top', fontsize=size, color='white', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.35,rounding_size=0.5', fc=color, ec='none'), zorder=10)


def _card(ax):
    for im in ax.get_images():
        a = np.asarray(im.get_array())
        if a.ndim == 3 and a.shape[2] >= 3:
            fl = a.dtype != np.uint8
            b = (np.clip(a, 0, 1) * 255).astype(np.uint8) if fl else a.copy()
            w = b[..., :3].astype(int).sum(2) > 735
            if b.shape[2] == 4:
                w |= b[..., 3] < 20; b[w, 3] = 255
            b[w, 0], b[w, 1], b[w, 2] = 246, 247, 251
            im.set_data(b)
    ax.add_patch(FancyBboxPatch((-0.02, -0.02), 1.04, 1.04, boxstyle='round,pad=0,rounding_size=0.04',
                                transform=ax.transAxes, fc=BG, ec=CARD_EDGE, lw=0.7, zorder=-5, clip_on=False))


def polish(fig, cards=True):
    axs = list(fig.axes)
    for ax in fig.axes:
        axs += [c for c in getattr(ax, 'child_axes', []) if c not in axs]
    for ax in axs:
        if getattr(ax, '_polished', False):
            continue
        ax._polished = True
        if not ax.axison:
            if cards and ax.get_images() and not getattr(ax, 'nocard', False):
                _card(ax)
            continue
        if getattr(ax, 'keep_style', False):
            continue
        ax.set_facecolor('white' if getattr(ax, 'white', False) else BG)
        from matplotlib.collections import PolyCollection
        from matplotlib.patches import Polygon, Rectangle
        for c in ax.collections:
            if isinstance(c, PolyCollection) and not isinstance(c, PathCollection):
                fc = c.get_facecolor()
                if len(fc) and abs(fc[0][0] - 0.945) < 0.01 and abs(fc[0][2] - 0.925) < 0.01:   # WASH band -> tinted band
                    c.set_facecolor(MC['Lu']); c.set_alpha(0.12)
        for pch in ax.patches:
            fc = pch.get_facecolor()
            if abs(fc[0] - 0.945) < 0.01 and abs(fc[2] - 0.925) < 0.01:
                pch.set_facecolor(MC['Lu']); pch.set_alpha(0.12)
        for s in ('left', 'bottom'):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0)
        ax.grid(True, color='white', lw=0.9)
        ax.set_axisbelow(True)
        for c in ax.collections:
            if isinstance(c, PathCollection) and len(c.get_offsets()) and not c.get_path_effects():
                c.set_path_effects(SHADOW)
        for l in ax.lines:
            if l.get_marker() not in (None, 'None', '', ' ') and l.get_linestyle() in ('None', '', ' ') and not l.get_path_effects():
                l.set_path_effects(SHADOW)
    return fig
