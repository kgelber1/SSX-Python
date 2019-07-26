#checker
import pandas as pd
import numpy as np
import os

"""
 This is a checker for the calibration coefficient for the ocean Optics
 system. It will compare the accuracy of the old coefficents to the new ones
"""

def checker(fileLoc):
    df = pd.read_excel(open(fileLoc,'rb'))

    waves = df.iloc[:,0]
    pixnums = df.iloc[:,1]

    runLen = df.shape[0]
    numBetter = 0
    for i in range (0, runLen):
        true = df.iloc[i,0]
        pix = df.iloc[i,1]

        orig =  df.iloc[0,2] + pix* df.iloc[1,2] + (pix**2)*df.iloc[2,2] + (pix**3)* df.iloc[3,2]
        new =  df.iloc[0,3] + pix* df.iloc[1,3] + (pix**2)*df.iloc[2,3] + (pix**3)* df.iloc[3,3]

        peorg = abs(orig - true)/true
        penew = abs(new - true)/true

        if(abs(peorg)> abs(penew)):
            numBetter +=1

    print("Passed " ,numBetter, "/", runLen)


def main():
    fileLoc = os.getcwd() + "Ca;ibration coefficents.xlsx"

    checker(fileLoc)

if __name__ == '__main__':
    main()
