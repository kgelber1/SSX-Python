from __future__ import division, print_function, absolute_import
import anim_bfield_merging_nBT as an
import anim_bfield_merging as a
import numpy as np


def main():
    """Just a place to specifiy variables"""
    day ='073019'

    first_shot = 12

    # last_shot = 44
    last_shot = 43
    bad_shots = [27]

    all_shots = np.arange(first_shot,last_shot+1)
    shots = [shot for shot in all_shots if shot not in bad_shots]
    sample_Freq = 5# sampling frequency - turn up for faster animations
    t0 = 20
    tf = 60
    for shot in shots:
        # will save each file in the Analyzed folder.
        print("shot", shot)
        try:
            an.run(day, shot, t0, tf, sample_Freq, show = False)
        except:
            a.run(day, shot, t0, tf, sample_Freq, show = False)


if __name__ == '__main__':
    main()
