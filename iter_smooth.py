import smooth as sm
from numpy import arange

#try 3 or 5
def iter_smooth(array,loops,window_len):
    for l in range(loops):
        array = sm.smooth(array, window_len,'flat')
    return array
# array=density,
