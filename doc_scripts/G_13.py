import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from colorspacious import cspace_converter

#.. Guides: Color Map Utilities - Miscellaneous
#   Named Colors 

# script editted from:
# https://matplotlib.org/3.1.0/gallery/color/named_colors.html

def plot_colortable(colors, title, emptycols=0):

    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12
    topmargin = 40

    # Sort colors by L* and name.
    by_lab = sorted((tuple(cspace_converter("sRGB1", "CAM02-UCS") (mcolors.to_rgb(color)) ),
                        name)
                    for name, color in colors.items())
    names = [name for lab, name in by_lab]
    namStr = lambda c,n : '{:.1f} {}'.format(c,n)
    nStrg = [namStr(val[0][0],val[1]) for val in by_lab]

    n = len(names)
    ncols = 4 - emptycols
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + margin + topmargin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-topmargin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24, loc="left", pad=10)

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        swatch_end_x = cell_width * col + swatch_width
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, nStrg[i], fontsize=12,
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, swatch_start_x, swatch_end_x,
                  color=colors[name], linewidth=18)

    return fig

plot_colortable(mcolors.BASE_COLORS, "Base Colors sorted by L*", emptycols=1)
plot_colortable(mcolors.TABLEAU_COLORS, "Tableau Palette sorted by L*", emptycols=2)
plot_colortable(mcolors.CSS4_COLORS, "CSS Colors sorted by L*")

plt.show()