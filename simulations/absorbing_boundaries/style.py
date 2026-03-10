import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import matplotlib 
from matplotlib.colors import ListedColormap
import numpy as np
import os
import matplotlib as mpl

scale = 1
colors = ['tab:blue','tab:orange','tab:green','tab:red']


# PRE figsize
mm_to_inch = 0.039370079
pw = 8.3
double_column_width = 7.48#pw-2
single_column_width = 3.543
figsize = { 'inch': {'column_width': single_column_width*scale, #(PRE), #89, #(nature) # 86, #(PRE)
                   'double_column_width': double_column_width*scale, #(PRE) # 183, # (nature) # 178, #(PRE)
                   'units': 'inches',
                   'golden_ratio': 0.618}, }



def load_preset(scale = 12/8,font_path='.'):
    
    '''
    font_manager._get_fontconfig_fonts.cache_clear()
    
    current_dir_path = '../'
    
    font_dirs = [font_path, ]
    font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
    for font_file in font_files:
        font_manager.fontManager.addfont(font_file)
    '''

    mpl.rcParams.update(mpl.rcParamsDefault)


    SMALL_SIZE = 28
    MEDIUM_SIZE = 28
    BIGGER_SIZE = 32
    MS = 10
    
    
    plt.rcParams['font.family']='serif'
    cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
    plt.rcParams['font.serif']=cmfont.get_name()
    plt.rcParams['mathtext.fontset']='cm'
    plt.rcParams['axes.unicode_minus']=False
    
    #plt.rcParams['text.usetex'] = True       
    #plt.rcParams['font.family'] = 'Helvetica Neue'
    #plt.rcParams['font.style'] = 'normal'
    plt.rcParams['font.size'] = 8*scale
    #plt.rcParams['font.weight'] = 'normal' 
    #plt.rcParams['axes.unicode_minus'] = False

    plt.rcParams['axes.linewidth'] = 0.75*scale
    plt.rcParams['axes.labelweight'] = 'normal'#'light'
    plt.rcParams['axes.titleweight'] = 'normal'#'light'
    plt.rcParams['axes.labelsize'] = 8*scale
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [figsize['inch']['column_width'], figsize['inch']['column_width']*figsize['inch']['golden_ratio']]
        
    plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams['figure.constrained_layout.h_pad'] = 0.01*scale
    plt.rcParams['figure.constrained_layout.w_pad'] = 0.01*scale
    plt.rcParams['image.aspect'] = 'auto'
    plt.rcParams['image.interpolation'] = 'None'
    plt.rcParams['image.origin'] = 'lower'
    plt.rcParams['lines.dash_capstyle'] = 'round'
    plt.rcParams['lines.solid_capstyle'] = 'round'
    plt.rcParams['lines.dash_joinstyle'] = 'bevel'
    plt.rcParams['lines.linewidth'] = 1.0*scale
    plt.rcParams['hatch.linewidth'] = 1.0*scale


    
    plt.rcParams['xtick.major.width'] = 0.75*scale
    plt.rcParams['xtick.major.pad'] = 2.5*scale
    plt.rcParams['xtick.major.size'] =  3.5*scale
    plt.rcParams['xtick.labelsize'] =  8*scale

    plt.rcParams['ytick.major.width'] = 0.75*scale
    plt.rcParams['ytick.major.pad'] = 2.5*scale
    plt.rcParams['ytick.major.size'] =  3.5*scale
    plt.rcParams['ytick.labelsize'] =  8*scale
    
    #plt.rcParams['savefig.transparent'] = True
    print('matplotlib preset loaded')

def cmap_nicify(cmap, transparent=False, idx_white = 0, size_white = None):
    """
    Make the bottom of the colormap white
    """
    register = False
    if type(cmap) == str:
        cmap = matplotlib.cm.get_cmap(cmap)
        register = True
    if size_white is None:
        size_white = cmap.N//5
        
    index_white = np.arange(2*size_white-1) - size_white + idx_white + 1
    curve = np.sin(np.linspace(-np.pi/2, np.pi/2, 2*size_white-1))**2
    clip = (index_white>=0)*(index_white<=(cmap.N-1))
    my_cmap_rgba = cmap(np.arange(cmap.N))
    # Set alpha
    my_cmap_rgba[:,-1][index_white[clip][0]:index_white[clip][-1]+1] = curve[clip]
    my_cmap_rgb = my_cmap_rgba.copy()
    
    if not transparent:
        # Transform alpha to color
        my_cmap_rgb[:,0] = (1-my_cmap_rgba[:,-1])*1 + my_cmap_rgba[:,-1]*my_cmap_rgba[:,0]
        my_cmap_rgb[:,1] = (1-my_cmap_rgba[:,-1])*1 + my_cmap_rgba[:,-1]*my_cmap_rgba[:,1]
        my_cmap_rgb[:,2] = (1-my_cmap_rgba[:,-1])*1 + my_cmap_rgba[:,-1]*my_cmap_rgba[:,2]

        my_cmap_rgb[:,-1] *= 0
        my_cmap_rgb[:,-1] += 1
    if register and not idx_white:
        matplotlib.cm.register_cmap(name=cmap.name + '_w', cmap=ListedColormap(my_cmap_rgb))
    else:
        return ListedColormap(my_cmap_rgb)
    
def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    """
    
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red','green','blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

