from __future__ import division, print_function, absolute_import
import numpy as np
import pandas as pd


"""
    This code will calulate the calibration coefficents for the Ocean View Optics System.
    Parameters:
        fileloc - the location of an excel sheet containg the information
                  I did hard code this in.
    Outputs:
        print out the outputs and also download another excel sheet containing the
        calibration coefficents for convient copying and pasting

    -- KG 2018
"""


def calculate(fileLoc):

    df = pd.read_excel(open(fileLoc,'rb'))

    #To get a better result, this cacluates several sets of calibration coefficents and then averages them
    i = 0
    while(i < df.shape[0]-4):
        x = df.iloc[i:i+4:,0:4]
        x.iloc[:,0] = 1
        y = df.iloc[i:i+4,0]
    #     print(x)

        ans = np.linalg.solve(x, y)
        if(i==0):
            cOne= ans[0]
            cTwo= ans[1]
            cThree= ans[2]
            cFour= ans[3]
        else:
            cOne= np.append(cOne, ans[0])
            cTwo= np.append(cTwo, ans[1])
            cThree= np.append(cThree, ans[2])
            cFour= np.append(cFour, ans[3])
        i+=1
    out = pd.DataFrame({'Data': [np.average(cOne), np.average(cTwo), np.average(cThree), np.average(cFour)]})
    writer = pd.ExcelWriter('NewCoefficents.xlsx', engine='xlsxwriter')
    out.to_excel(writer, sheet_name='Sheet1')
    writer.save()

    print("The Intercept is " , np.average(cOne))
    print("C1 is " , np.average(cTwo))
    print("C2 is " , np.average(cThree))
    print("C3 is " , np.average(cFour))
    #need to solve the equation lamda = I + c1P + c2P**2 + c3p**3

def main():
    file = r'C:\Users\Katie\Documents\Book2.xlsx'
    calculate(file)


if __name__ == '__main__':
    main()
