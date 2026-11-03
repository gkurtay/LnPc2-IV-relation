"""Common figure style of the article (palette, panel background, marker shading, image cards).

polish(fig) is applied to every figure just before it is saved: plot panels get a light background
with a white grid, isolated markers a soft shadow, and structure images a rounded card."""
import numpy as np
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch
from matplotlib.collections import PathCollection
from matplotlib.image import AxesImage

BG = 'white'            # panel background (v30: white)
GRID = '#e8eaef'; SPINE = '#b9bdc7'
CARD_EDGE = '#dde1e8'
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
            b[w, 0], b[w, 1], b[w, 2] = 255, 255, 255
            im.set_data(b)
    # v31: no card frame around structure images (white background only)


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
                for t in ax.texts:   # v30: panel letter outside the card, above its upper left corner
                    if t.get_fontweight() in ('bold', 700) and abs(t.get_fontsize() - 9) < 0.1 and len(t.get_text()) == 1:
                        t.set_position((-0.02, 1.035)); t.set_va('bottom'); t.set_ha('left')
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
            ax.spines[s].set_visible(True); ax.spines[s].set_color(SPINE); ax.spines[s].set_linewidth(0.6)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=2.2, width=0.6, color=SPINE)
        ax.grid(True, color=GRID, lw=0.55)
        ax.set_axisbelow(True)
        for c in ax.collections:   # v30: no marker shadows
            if c.get_path_effects() and any(isinstance(e, pe.SimplePatchShadow) for e in c.get_path_effects()):
                c.set_path_effects([])
    return fig
