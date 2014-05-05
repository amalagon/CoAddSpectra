#chitanda.py
#co-adds the spectra made by hyouka.py
import os
import csv
from scipy import signal
from matplotlib.pylab import *
import pandas as pd
#from collections import defaultdict
#from collections import OrderedDict
import glob
#d=defaultdict(list)

#rfile=open('residual.txt','r')
#residual = [float(rows.split()[0]) for rows in rfile.readlines()]

freqsep=34013.6
coaddfreqlist=[33.4*np.power(10,9)+freqsep*j for j in range(int(np.floor((34.53-33.4)*np.power(10,9)/freqsep)))]
coaddpowlist= [[] for _ in range(len(coaddfreqlist))]
coadderrlist= [[] for _ in range(len(coaddfreqlist))]
boollist=[False]*len(coaddfreqlist)
coaddnumblist= [[] for _ in range(len(coaddfreqlist))]
Goodfilesdict=dict()

with open('goodspectralist.txt','rb') as goodfiles:
	Goodfilesdict={rows.split('_')[0].split('/')[1] : rows.split(' ')[1].split('\r')[0] for rows in goodfiles.readlines()}

datadir='E:/Documents/Documents/PythonScripts'
pgmdir='E:/Documents/Documents/PythonScripts'
os.chdir(datadir)
#notfound=open(pgmdir+'/'+'notfoundfiles.txt','w')
daylist= os.walk('.').next()[1]
#d=defaultdict(list)
#e=defaultdict(list)
toohighflucts=['12-17-19-33-40','12-20-16-14-54', '01-31-01-58-28']
donotgolist=['12-18-14-39-57_5e4','12-20-16-14-54_5e5','Dec04-06-49-51_5e3','Dec04-06-49-51_5e4','Dec04-07-37-58_5e3','Dec04-07-37-58_5e4','Dec04-08-09-08_5e3']
for day in daylist:
	os.chdir(datadir+'/'+day)
	print datadir+'/'+day
	timestampdirlist=os.walk('.').next()[1]
	for timestampdir in timestampdirlist:
		if  'Nov22' != day and 'bins2048' not in timestampdir and 'Nov19' != day and timestampdir not in donotgolist:
			#os.chdir(datadir+'/'+day+'/'+timestampdir)
			#print datadir+'/'+day+'/'+timestampdir
			timestamp=timestampdir.split('_',1)[0]
			N=float(timestampdir.split('_',1)[1])
			try:
				os.chdir(pgmdir+'/'+day+'/'+timestampdir+'/'+timestamp+'avgsubfiles')
				stopindex = int(Goodfilesdict[timestamp])
				files=glob.glob('*'+str(stopindex)+'.txt')
				data=np.genfromtxt(files[0])
				Xbband=data[93:170,0]
				Xtruefreq=data[93:170,1]
				flucts=data[93:170,2]-signal.wiener(data[93:170,2],mysize=15)
				#flucts=data[93:170,2]
				#plot(data[93:170,2],'bo')
				#plot(signal.wiener(data[93:170,2],mysize=15),'ro')
				#wi=signal.wiener(flucts)
				if timestamp in toohighflucts:
					print timestamp
				else:
					fit=data[93:170,3]
					error=data[93:170,4]
					for j in range(len(Xtruefreq)):
						id = int(np.floor((Xtruefreq[j]-33.4*np.power(10,9))/freqsep))
						boollist[id]=True
						coaddpowlist[id].append(flucts[j])
						coadderrlist[id].append(error[j])
						coaddnumblist[id].append(stopindex*N)
						#d[Xtruefreq[j]].append(flucts[j])
						#e[Xtruefreq[j]].append(error[j])
			except KeyError:
				print "timestamp %s not found in spreadsheet" % timestamp
				#print>>notfound, timestamp
			except ValueError:
				print "probably no Cavity Freq entered in spreadsheet for timestamp %s" % timestamp

coaddweightlist=[[] for _ in range(len(coaddfreqlist))]
for j in range(len(coaddfreqlist)):
	if boollist[j]==True:
		norm=np.sum(np.square(np.reciprocal(coadderrlist[j])))
		for val in coadderrlist[j]:
			coaddweightlist[j].append(np.square(np.reciprocal(val))/norm)
powlist=[]
errlist=[]
freqlist=[]
numblist=[]
for j in range(len(coaddfreqlist)):
	if boollist[j]==True:
		freqlist.append(coaddfreqlist[j])
		powlist.append(np.average(coaddpowlist[j],weights=coaddweightlist[j]))
		errlist.append(np.sqrt(np.sum(np.square(np.multiply(coaddweightlist[j],coadderrlist[j])))))
		numblist.append(np.sum(coaddnumblist[j]))

#wavgfile = open (pgmdir+'/wavgcoaddfile.txt','w')
for j in range(len(freqlist)):
	errorbar(freqlist[j]/np.power(10,9),powlist[j],yerr=errlist[j],fmt='bo')
	#print>>wavgfile, freqlist[j], powlist[j], errlist[j], numblist[j]
#wavgfile.close()
xlabel('Frequency (GHz)')
ylabel('Fluctations of Noise Power (mW/Hz)')
title('Co-added Spectra')
grid(True)
# show()
show()
os.chdir('E:/Documents/Documents/PythonScripts')