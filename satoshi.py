#chitanda.py
#co-adds the spectra made by hyouka.py
#divides out axion shape to give you the co=added spectra in units of axion excess power
import os
import csv
from scipy import signal
from matplotlib.pylab import *
import pandas as pd
#from collections import defaultdict
#from collections import OrderedDict
import glob
BWdict=dict()
Cavitydict=dict()
Qdict=dict()
#d=defaultdict(list)
with open('pythonreally.csv','rb') as infile:
	reader=csv.reader(infile)
	for rows in reader:
		Cavitydict[rows[1].strip()]=rows[9].strip()
		#LOdict[rows[1].strip()]=rows[10].strip()
		BWdict[rows[1].strip()]=rows[11].strip()
		Qdict[rows[1].strip()]=rows[12].strip()
#rfile=open('residual.txt','r')
#residual = [float(rows.split()[0]) for rows in rfile.readlines()]

freqsep=34013.6
coaddfreqlist=[33.4*np.power(10,9)+freqsep*j for j in range(int(np.floor((34.53-33.4)*np.power(10,9)/freqsep)))]
coaddpowlist= [[] for _ in range(len(coaddfreqlist))]
coadderrlist= [[] for _ in range(len(coaddfreqlist))]
boollist=[False]*len(coaddfreqlist)
coaddnumblist= [[] for _ in range(len(coaddfreqlist))]
Goodfilesdict=dict()

def lorentzianshape(whereyouare,CavityFreq,BW):
    return 1/(1+4*(whereyouare-CavityFreq)**2/BW**2)

def correction(whereyouare, Q):
	#multiply this by 8.3*10^21 eV^4 to get the correct prefactor for the signal
	#whereyouare is in units of GHz
	bare =Q/np.power(10,4)*(33.85)/whereyouare
	return bare

def pluslor(whereyouare, CavityFreq, BW, Q):
	#this is the total correction including the lorentzian shape of axion response
	return lorentzianshape(whereyouare, CavityFreq, BW)*correction(whereyouare, Q)

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
special=[]
donotgolist=['12-18-14-39-57_5e4','12-20-16-14-54_5e5','Dec04-06-49-51_5e3','Dec04-06-49-51_5e4','Dec04-07-37-58_5e3','Dec04-07-37-58_5e4','Dec04-08-09-08_5e3']
for day in daylist:
	os.chdir(datadir+'/'+day)
	print datadir+'/'+day
	timestampdirlist=os.walk('.').next()[1]
	for timestampdir in timestampdirlist:
		if 'Nov20' != day and 'Nov22' != day and 'bins2048' not in timestampdir and 'Nov19' != day and timestampdir not in donotgolist:
			#os.chdir(datadir+'/'+day+'/'+timestampdir)
			#print datadir+'/'+day+'/'+timestampdir
			timestamp=timestampdir.split('_',1)[0]
			N=float(timestampdir.split('_',1)[1])
			#try:
				#Q=float(Qdict[timestamp])
			#except ValueError:
			#	print 'no Q for timestamp %s' % timestamp
			try:
				CavityFreq=float(Cavitydict[timestamp])
			except ValueError:
				print 'no Cavity for timestamp %s' % timestamp
			try:
				BW=.004
				os.chdir(pgmdir+'/'+day+'/'+timestampdir+'/'+timestamp+'avgsubfiles')
				stopindex = int(Goodfilesdict[timestamp])
				files=glob.glob('*'+str(stopindex)+'.txt')
				data=np.genfromtxt(files[0])
				Xbband=data[93:170,0]
				Xtruefreq=data[93:170,1]
				#flucts=data[93:170,2]-signal.wiener(data[93:170,2],mysize=15)
				flucts=data[93:170,2]
				#wi=signal.wiener(flucts)
				if timestamp in toohighflucts:
					print timestamp
				else:
					Pa=4.52/np.power(10,9)/np.power(10,4) #for g_11=7.27 in mW/Hz
					fit=data[93:170,3]
					error=data[93:170,4]
					#plot(pluslor(Xtruefreq,CavityFreq,BW,Q),'ro')
					for j in range(len(Xtruefreq)):
						#plot(pluslor(Xtruefreq[j],CavityFreq,BW,Q),'ro')
						#if N*stopindex > 0.2*np.power(10,9):
						#	print timestamp+'has large number'
						#	special.append(Xtruefreq[j])
						#if Xtruefreq[j]>33.35*np.power(10,9):
						id = int(np.floor((Xtruefreq[j]-33.4*np.power(10,9))/freqsep))
						boollist[id]=True
						#coaddpowlist[id].append(flucts[j]*freqsep/np.power(10,3)/pluslor(Xtruefreq[j]/np.power(10,9),CavityFreq,BW,Q))
						#coadderrlist[id].append(error[j]*freqsep/np.power(10,3)/pluslor(Xtruefreq[j]/np.power(10,9),CavityFreq,BW,Q))
						coaddpowlist[id].append(flucts[j]/lorentzianshape(Xtruefreq[j]/np.power(10,9),CavityFreq,BW))
						coadderrlist[id].append(error[j]/lorentzianshape(Xtruefreq[j]/np.power(10,9),CavityFreq,BW))
						
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
		#powlist.append(np.average(coaddpowlist[j],weights=coaddweightlist[j]))
		powlist.append(np.average(coaddpowlist[j],weights=coaddnumblist[j]/np.sum(coaddnumblist[j])))
		#errlist.append(np.sqrt(np.sum(np.square(np.multiply(coaddweightlist[j],coadderrlist[j])))))
		errlist.append(np.sqrt(np.sum(np.square(coadderrlist[j]))))

		numblist.append(np.sum(coaddnumblist[j]))
#wavgfile = open (pgmdir+'/wavgcoaddfileexcspwr.txt','w')
for j in range(len(freqlist)):
	#if freqlist[j] in special:
		#pass
		#errorbar(freqlist[j]/np.power(10,9),powlist[j],yerr=errlist[j],fmt='r-')
		#plot(freqlist[j]/np.power(10,9),powlist[j],'r-')
	#else:
		#pass
		#errorbar(freqlist[j]/np.power(10,9),powlist[j],yerr=errlist[j],fmt='b-')
	errorbar(freqlist[j]/np.power(10,9),powlist[j],yerr=errlist[j],fmt='b-')
	#print>>wavgfile, freqlist[j], powlist[j], errlist[j], numblist[j]
#wavgfile.close()
xlabel('Frequency (GHz)')
ylabel('Power Fluctuations ($P_{\text{noise}}$) (mW/Hz)')
#text(0.5,0.5,'$g = 8.2')
#xlim([33.88,34.53])
title('Co-added Spectra')
grid(True)
#savefig(pgmdir+'/'+'whatthe.png',bbox_inches='tight')
show()
os.chdir('E:/Documents/Documents/PythonScripts')