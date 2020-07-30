import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from collections import OrderedDict
import s3dlib.cmap_utilities as cmu
from cmap_xtra import Lab_cmap_gradient

cmaps = OrderedDict()

cmaps['Binary'] = [
            'binaryDefault', 'binaryRdBu', 'binaryBuYw', 'binaryBuGnLw' ]

cmu.binary_cmap( name='binaryDefault' )
cmu.binary_cmap( [1,0,0], [0,0,1], 'binaryRdBu'  )
cmu.binary_cmap( 'blue', 'yellow', 'binaryBuYw'  )
cmu.binary_cmap( [0,0,.5], [0,.5,0], 'binaryBuGnLw'  )

cmaps['RGB Gradient'] = [
            'rgbDefault', 'mpDefault',
            'RdToGn', 'YWToBu', 'cardboard', 'cardboardMrrd',
            'bnShade', 'WtToBk', 'BkToWtMrrd', 'WtToBkMrrd' ]

cmu.rgb_cmap_gradient( name='rgbDefault')
cmu.rgb_cmap_gradient( 'black', [0.122, 0.467, 0.706], 'mpDefault' )

cmu.rgb_cmap_gradient( 'red', 'green', 'RdToGn' )
cmu.rgb_cmap_gradient( 'yellow', 'blue', 'YWToBu' )
cmu.rgb_cmap_gradient( [0.25,0.15,0], [1,.9,.75], 'cardboard' )
cmu.rgb_cmap_gradient( [0.25,0.15,0], [1,.9,.75], 'cardboardMrrd', mirrored=True )

cmu.rgb_cmap_gradient( [0.25,0.15,0], 'white', 'bnShade' )
cmu.rgb_cmap_gradient( 'white', 'black', 'WtToBk' )
cmu.rgb_cmap_gradient( 'black', 'white', 'BkToWtMrrd', mirrored=True )
cmu.rgb_cmap_gradient( 'white', 'black', 'WtToBkMrrd', mirrored=True )

cmaps['Hue Gradient'] = [
            'hsvDefault',
            'hsvRdToGn', 'hsvRdFromGn', 'hsvYwToBu','hsvYwFromBu',
            'RdToRd', 'RdFromRd', 'RdToRdMrrd', 'RdFromRdMrrd',
            'BlToBl', 'BlFromBl', 'BlToBlMrrd', 'BlFromBlMrrd' ]

cmu.hsv_cmap_gradient( name='hsvDefault')

cmu.hsv_cmap_gradient( [0,1,1], [0.333,1,1], 'hsvRdToGn' )
cmu.hsv_cmap_gradient( [1,1,1], [0.333,1,1], 'hsvRdFromGn' )
cmu.hsv_cmap_gradient( [0.166,1,1], [0.666,1,1], 'hsvYwToBu' )
cmu.hsv_cmap_gradient( [1.166,1,1], [0.666,1,1], 'hsvYwFromBu' )

cmu.hsv_cmap_gradient( [0,1,1], [1,1,1], 'RdToRd' )
cmu.hsv_cmap_gradient( [1,1,1], [0,1,1], 'RdFromRd' )
cmu.hsv_cmap_gradient( [0,1,1], [1,1,1],'RdToRdMrrd', mirrored=True )
cmu.hsv_cmap_gradient( [1,1,1], [0,1,1],'RdFromRdMrrd', mirrored=True )

cmu.hsv_cmap_gradient( [0.666,1,1], [1.666,1,1], 'BlToBl' )
cmu.hsv_cmap_gradient( [1.666,1,1], [0.666,1,1], 'BlFromBl' )
cmu.hsv_cmap_gradient( [0.666,1,1], [1.666,1,1], 'BlToBlMrrd', mirrored=True )
cmu.hsv_cmap_gradient( [1.666,1,1], [0.666,1,1], 'BlFromBlMrrd', mirrored=True )

cmaps['Saturation and Value Gradient'] = [
            'RdToBk', 'BkFromRd', 'RdToBkMrrd', 'BkFromRdMrrd',
            'RdToWt', 'WtFromRd', 'RdToWtMrrd', 'WtFromRdMrrd' ]

cmu.hsv_cmap_gradient( [0,1,1], [0,1,0], 'RdToBk' )
cmu.hsv_cmap_gradient( [0,1,0], [0,1,1], 'BkFromRd' )
cmu.hsv_cmap_gradient( [0,1,1], [0,1,0], 'RdToBkMrrd', mirrored=True )
cmu.hsv_cmap_gradient( [0,1,0], [0,1,1], 'BkFromRdMrrd', mirrored=True )

cmu.hsv_cmap_gradient( [0,1,1], [0,0,1], 'RdToWt' )
cmu.hsv_cmap_gradient( [0,0,1], [0,1,1], 'WtFromRd' )
cmu.hsv_cmap_gradient( [0,1,1], [0,0,1], 'RdToWtMrrd', mirrored=True )
cmu.hsv_cmap_gradient( [0,0,1], [0,1,1], 'WtFromRdMrrd', mirrored=True )

cmaps['Mirrored and Reversed'] = [
            'viridis', 'plasma', 'magma', 'cividis',
            'viridis_r', 'plasma_r', 'magma_r', 'cividis_r',
            'viridis_m', 'plasma_m', 'magma_m', 'cividis_m',
            'viridis_mr', 'plasma_mr', 'magma_mr', 'cividis_mr' ]
            
cmu.reversed_cmap('viridis')
cmu.reversed_cmap('plasma')
cmu.reversed_cmap('magma')
cmu.reversed_cmap('cividis')

cmu.mirrored_cmap('viridis')
cmu.mirrored_cmap('plasma')
cmu.mirrored_cmap('magma')
cmu.mirrored_cmap('cividis')

cmu.mirrored_cmap('viridis', rev=True)
cmu.mirrored_cmap('plasma', rev=True)
cmu.mirrored_cmap('magma', rev=True)
cmu.mirrored_cmap('cividis', rev=True)

cmaps['Cyclic and Named Color'] = [
            'red3', 'blue2R', 'goldteal', 'plgrnOrchR',
            'pgrdkcy', 'hpnktrq', 'pgrdkcy_RGB', 'hpnktrq_RGB' ]

cmu.hsv_cmap_gradient( [0,1,1], [3,1,1], 'red3' )
cmu.hsv_cmap_gradient( [2.666,1,1], [0.666,1,1], 'blue2R' )
cmu.hsv_cmap_gradient( 'gold','teal', 'goldteal' )
cmu.hsv_cmap_gradient( '+palegreen','orchid', 'plgrnOrchR' )

cmu.hsv_cmap_gradient( 'palegoldenrod','darkcyan', 'pgrdkcy' )
cmu.hsv_cmap_gradient( 'hotpink','turquoise', 'hpnktrq' )

cmu.rgb_cmap_gradient( 'palegoldenrod','darkcyan', 'pgrdkcy_RGB' )
cmu.rgb_cmap_gradient( 'hotpink','turquoise', 'hpnktrq_RGB' )


cmaps['Hue HSV'] = [
            'hue_05', 'hsvDefault', 'HUEdefault','hue_80',
            'hueYwToYw8', 'hueYwToBu8', 'hueYwFromBu8', 'hueYwFromBu5' ]

cmu.hue_cmap(0.5,name='hue_05' )
cmu.hue_cmap(    name='HUEdefault' )
cmu.hue_cmap(8.0,name='hue_80')
cmu.hue_cmap(0.5,name='hue_05' )
cmu.hue_cmap(8.0, 0.166, name='hueYwToYw8')
cmu.hue_cmap(8.0,'y', 'b',    'hueYwToBu8')
cmu.hue_cmap(8.0,'+y','b',    'hueYwFromBu8')
cmu.hue_cmap(0.5,'+y','b',    'hueYwFromBu5')


cmaps['Smoothed HSV'] = [
            'RG', 'RG_s', 'GB','GB_s',
            'BR', 'BR_s']

cmu.hsv_cmap_gradient( [0,1,1],      [0.333,1,.65], 'RG' )
cmu.hsv_cmap_gradient( [0,1,1],      [0.333,1,.65], 'RG_s', smooth=1.6 )
cmu.hsv_cmap_gradient( [.333,1,.34], [0.666,1,1],   'GB' )
cmu.hsv_cmap_gradient( [.333,1,.34], [0.666,1,1],   'GB_s', smooth=1.4 )
cmu.hsv_cmap_gradient( [0.666,1,1],  [1,1,.52],     'BR' )
cmu.hsv_cmap_gradient( [0.666,1,1],  [1,1,.52],     'BR_s', smooth=1.6 )


cmaps['Lab Gradient'] = [
            'magma', 'viridis', 'lab_indg','lab_marn',
            'lab_dgrn', 'LabDefault', 'RdToWt' ]

Lab_cmap_gradient('indigo',   'aqua',        'lab_indg')
Lab_cmap_gradient('maroon',   'yellow',      'lab_marn')
Lab_cmap_gradient('darkgreen','lemonchiffon','lab_dgrn')
Lab_cmap_gradient( name='LabDefault')



cmaps['Stitch'] = [
            'fWt_2', 'inferno_3', 'RdCy_4',
            'Diverging*', 'Cyclic*', 'stchG' ]

cmA = cmu.hsv_cmap_gradient('firebrick','wheat',smooth=1.5)
cmB = cmu.hsv_cmap_gradient('wheat','teal',smooth=1.5)
cmu.stitch_cmap(cmA,cmB,name='fWt_2')

cmu.mirrored_cmap('inferno')
testMap3 = cmu.stitch_cmap('inferno','inferno_m','inferno_r',bndry=[0.2,0.6],name='inferno_3' )

cmu.reversed_cmap('RdToBk')
cmu.hsv_cmap_gradient( [0.5,1,1], [0.5,1,0], 'CyToBk' )
cmu.reversed_cmap('CyToBk')
cmu.stitch_cmap('RdToBk_r', 'RdToBk','CyToBk_r', 'CyToBk',bndry=[0.3,0.5,0.8], name='RdCy_4' )

cmC = Lab_cmap_gradient('olive','paleturquoise')
cmD = Lab_cmap_gradient('paleturquoise','mediumslateblue')
cmu.stitch_cmap(cmC,cmD,name='Diverging*')

cmE,cmF = Lab_cmap_gradient('black','mediumvioletred'), Lab_cmap_gradient('mediumvioletred','white'),
cmG,cmH = Lab_cmap_gradient('white','green'), Lab_cmap_gradient('green','black')
cmu.stitch_cmap(cmE,cmF,cmG,cmH,name='Cyclic*')

cmu.hsv_cmap_gradient( [0,1,.8],        [0,0,1],       'redG')
cmu.hsv_cmap_gradient( 'saddlebrown',   'yellow',      'yelG')
cmu.hsv_cmap_gradient( 'darkgreen',     'lime',        'grnG')
cmu.hsv_cmap_gradient( 'darkslategray', 'cyan',        'cynG')
cmu.hsv_cmap_gradient( 'midnightblue',  'deepskyblue', 'bluG')
cmu.hsv_cmap_gradient( [0.833,1,0.1],   'magenta',     'mgnG')
cmu.stitch_cmap( 'mgnG','bluG','cynG', 'grnG','yelG', 'redG', name='stchG' )


'''
Following code is copied directly from the Matplotlib tutorial:
"Choosing Colormaps in Matplotlib"
https://matplotlib.org/tutorials/colors/colormaps.html#sphx-glr-tutorials-colors-colormaps-py
'''
nrows = max(len(cmap_list) for cmap_category, cmap_list in cmaps.items())
nrows = 18    # set to get same size as matplotlib examples ;-)
gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))

def plot_color_gradients(cmap_category, cmap_list, nrows):
    fig, axes = plt.subplots(nrows=nrows)
    fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
    axes[0].set_title(cmap_category + ' colormaps', fontsize=14)

    for ax, name in zip(axes, cmap_list):
        ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
        pos = list(ax.get_position().bounds)
        x_text = pos[0] - 0.01
        y_text = pos[1] + pos[3]/2.
        fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axes:
        ax.set_axis_off()

#i = 1 
for cmap_category, cmap_list in cmaps.items():
    plot_color_gradients(cmap_category, cmap_list, nrows)
    #filename = 'tutorial_images/cmapExamples_' + str(i) + '.png'
    #i += 1
    #print('cat files:',filename)
    #plt.savefig(filename)

plt.show()