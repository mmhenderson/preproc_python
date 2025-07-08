import numpy as np
import scipy.io as spio

def load_mat_behav_data(filename, varname='TheData'):
    
    """
    A way to load .mat files that contain lists of data structures
    into python as lists of dicts. This specifically works for the format of
    data used in shapedim, might need to be modified for other experiments!
    """
    
    TheData = spio.loadmat(filename, squeeze_me=True, struct_as_record=False)[varname]
    try:
        TheData = [_todict(d) for d in TheData]
    except:
        # the above will throw error if only one run
        TheData = [_todict(TheData)]
        
    return TheData

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    d = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            d[strg] = _todict(elem)
        elif isinstance(elem, np.ndarray):
            d[strg] = _tolist(elem)
        else:
            d[strg] = elem
    return d

def _tolist(ndarray):
    '''
    A recursive function which constructs lists from cellarrays
    (which are loaded as numpy ndarrays), recursing into the elements
    if they contain matobjects.
    '''
    elem_list = []
    for sub_elem in ndarray:
        if isinstance(sub_elem, spio.matlab.mio5_params.mat_struct):
            elem_list.append(_todict(sub_elem))
        elif isinstance(sub_elem, np.ndarray):
            elem_list.append(_tolist(sub_elem))
        else:
            elem_list.append(sub_elem)
    return elem_list



def read_xfm_mat(file):

    with open(file, 'r') as f:
        lines = f.readlines()
        mat = []
        for line in lines:
            parts = line.split(' ')
            parts = [p for p in parts if (p!='') and (p!='\n')]
            nums = [float(p) for p in parts]
            mat += [nums]
    mat = np.array(mat)
    return mat


def write_xfm_mat(mat, file):

    with open(file, 'w') as f:
        for ii in range(mat.shape[0]):
            line = ''
            for jj in range(mat.shape[1]):
                line += str(mat[ii][jj]) + '  '
            line += '\n'
            f.write(line)