import anim_bfield_merging_nBT as an
import numpy as np


def main():
    """Just a place to specifiy variables"""
    day ='072419'

    first_shot = 35
    # last_shot = 44
    last_shot = 43
    bad_shots = [13, 14, 18]

    all_shots = np.arange(first_shot,last_shot+1)
    shots = [shot for shot in all_shots if shot not in bad_shots]
    sample_Freq = 5# sampling frequency - turn up for faster animations
    t0 = 25
    tf = 65
    for shot in shots:
        # will save each file in the Analyzed folder.
        print("shot", shot)
        an.run(day, shot, t0, tf, sample_Freq, show = False)

if __name__ == '__main__':
    main()
