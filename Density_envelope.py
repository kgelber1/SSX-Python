# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 16:32:42 2017

@author: Manjit
"""

import matplotlib.pyplot as plt
import ssx_basic as ssxd
#import basic_ids_1 as idsd
import numpy as np

def plot_envs(day, shots, scope_used = '1'):

	path = 'Data\\2019\\'+day+'\\'+day+'-Analysed\\'
	setting = '-envelope'

	for shot in shots:
		print('On Shot',shot)
		fig=plt.figure(num=1)#,figsize=(8.7,8.5),facecolor='w',edgecolor='k')#, dpi=600)
		plt.clf()
		fig.subplots_adjust(top=0.96, bottom=0.1, left = 0.14, right=0.96, hspace=0)
		ax=plt.subplot(1,1,1)
		#plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)

		density = np.zeros([20002])
		dens = ssxd.interferometer(day+'r'+str(shot), [1, 1], scope = scope_used, showPlot=False)
		plt.scatter(dens.signal1, dens.signal2, cmap='Blues')#,lw= 1)

		plt.xlabel(r'Signal 1',fontsize=12, weight='bold')
		plt.ylabel(r'Signal 2',fontsize=12, weight='bold')
		'''
		tp = 60
		d=idsd.ids(day+'r'+str(shot))
		d.processIDS(times=[-2,125])
		d.plotRawIDStp(tp)
		numChannels = d.data.shape[0]
		plt.plot(numChannels, d.data[:,tp], 'bo')
		'''

		plt.title(day+'r'+str(shot)+ ': envelope', fontsize=12, weight='bold')
	#	ax.get_yaxis().set_label_coords(-0.11,0.6) # for aligning the y-labels in one line
		plt.setp(ax.spines.values(), linewidth=1)#changing the axis linewidth
		ax.tick_params(axis='both', direction='in', length=7, width =1, labelsize = 12)
		ax.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
		if max(dens.signal1)>max(dens.signal2):
		    a = max(dens.signal1)+0.2
		    plt.xlim(0, a)
		    plt.ylim(0, a)
		if max(dens.signal2)>max(dens.signal1):
		    a = max(dens.signal2)+0.2
		    plt.xlim(0, a)
		    plt.ylim(0, a)

		ax.axis('square')
		plt.show()
		# fName = path+day+'r'+str(shot)+setting+'.png'
		# fig.savefig(fName,dpi=150,facecolor='w',edgecolor='k')

def main():
	day = '072419'#'052417'#'081717'
	shot_range = [1,5]

	plot_envs(day,np.arange(shot_range[0],shot_range[-1]+1))


if __name__ == "__main__":
	main()
