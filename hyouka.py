#hyouka.py
#orders the fft files, reads in all data within Separate directory
#averages up to the number in goodspectralist.txt,
#subtracts the mean, cuts out dc noise and wings, and saves plot and text file 
# columns are truefreq, meansub, and movingavg
#reads in the inputspreadsheet with the LO and CavityFreq info and saves the true freq of the bins
#if it can't find the data timestamp in the spreadsheet, saves that timestamp to notfoundfile.txt 
from matplotlib.pylab import *
import glob
import os
import csv
import re
import pandas as pd
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def stderror(N, input):
	return input/np.sqrt(N)

def holt_winters_second_order_ewma( x, span, beta ):
    N = x.size
    alpha = 2.0 / ( 1 + span )
    s = np.zeros(( N, ))
    b = np.zeros(( N, ))
    s[0] = x[0]
    for i in range( 1, N ):
        s[i] = alpha * x[i] + ( 1 - alpha )*( s[i-1] + b[i-1] )
        b[i] = beta * ( s[i] - s[i-1] ) + ( 1 - beta ) * b[i-1]
    return s
	
ticklabel_format(style='sci', axis='x', scilimits=(0,0))

LOdict=dict()
BWdict=dict()
Cavitydict=dict()
Goodfilesdict=dict()

with open('pythonreally.csv','rb') as infile:
	reader=csv.reader(infile)
	for rows in reader:
		Cavitydict[rows[1].strip()]=rows[9].strip()
		LOdict[rows[1].strip()]=rows[10].strip()
		BWdict[rows[1].strip()]=rows[11].strip()
		
with open('goodspectralist.txt','rb') as goodfiles:
	Goodfilesdict={rows.split('_')[0].split('/')[1] : rows.split(' ')[1].split('\r')[0] for rows in goodfiles.readlines()}

with open("Processed235spectraDec04-11-33-24.txt","rb") as infile:
	mylist= infile.readlines()
bband=[rows.split()[0] for rows in mylist[1:]]
power=[rows.split()[1] for rows in mylist[1:]]

with open("E:/Documents/Documents/Mathematica/WeightedMean/mathematicameansubdataavg5.txt","rb") as mathavg5file:
	mathavg5list=mathavg5file.readlines()
	freq = [float(rows.split()[0]) for rows in mathavg5list]
	poweravg5 = [float(rows.split()[1]) for rows in mathavg5list]
	
datadir='E:/Documents/Documents/Separate'
pgmdir='E:/Documents/Documents/PythonScripts'
os.chdir(datadir)
#notfound=open(pgmdir+'/'+'notfoundfiles.txt','w')
daylist= os.walk('.').next()[1]
def f(day):
	os.chdir(datadir+'/'+day)
	print datadir+'/'+day
	timestampdirlist=os.walk('.').next()[1]
	for timestampdir in timestampdirlist:
		if 'heezle' not in timestampdir:
			os.chdir(datadir+'/'+day+'/'+timestampdir)
			print datadir+'/'+day+'/'+timestampdir
			files = sorted(glob.glob('*'),key=numericalSort)
			timestamp=timestampdir.split('_',1)[0]
			N=float(timestampdir.split('_',1)[1])
			try:
				LOFreq=float(LOdict[timestamp])
				CavityFreq=float(Cavitydict[timestamp])
				stopindex = int(Goodfilesdict[timestamp])
				#stopindex=235
				print "timestamp %s has LO %e , Cavity %e , stopindex %d " % (timestamp, LOFreq, CavityFreq, stopindex)
				pltdir=pgmdir+'/'+day+'/'+timestampdir+'/'+timestamp+'avgsubplots'
				filedir=pgmdir+'/'+day+'/'+timestampdir+'/'+timestamp+'avgsubfiles'
				if not os.path.exists(pltdir):
					os.makedirs(pltdir)
				if not os.path.exists(filedir):
					os.makedirs(filedir)
				data0=np.genfromtxt(files[0])
				y0=data0[:,1]
				x0=data0[:,0]
				summand=np.zeros(186)
				before=np.zeros(186)
				variance=np.zeros(186)
				smooth=np.zeros(186)
				i=0
				while i < stopindex:
					index = files[i].split(timestamp)[0].split('FFTS')[1]
					data=np.genfromtxt(files[i])
					#summand += data[:,1]
					i+=1
				#print i
				#summand/= stopindex
					cut=7
					#span=2.001
					#beta=0.99
					X,Y=x0,data[:,1] 
					Xcutleft,Ycutleft=X[30:137],Y[30:137]
					Xcutright,Ycutright=X[157:264],Y[157:264]
					Xcutcombined,Ycutcombined=np.concatenate((Xcutleft,Xcutright)),np.concatenate((Ycutleft,Ycutright))				
					Ymovingavgleft=pd.rolling_mean(Ycutleft,5, center=True)
					Ymovingavgright=pd.rolling_mean(Ycutright,5, center=True)
					#Ymovingavgleft=holt_winters_second_order_ewma(Ycutleft, span, beta)
					#Ymovingavgright=holt_winters_second_order_ewma(Ycutright, span, beta)
					# fwdl = pd.ewma( Ycutleft, span=2.1 ) # take EWMA in fwd direction
					# bwdl = pd.ewma( Ycutleft[::-1], span=2.1 ) # take EWMA in bwd direction
					# cl = np.vstack(( fwdl, bwdl[::-1] )) # lump fwd and bwd together
					# Ymovingavgleft = np.mean( cl, axis=0 ) # average
					# fwdr = pd.ewma( Ycutright, span=2.1 ) # take EWMA in fwd direction
					# bwdr = pd.ewma( Ycutright[::-1], span=2.1 ) # take EWMA in bwd direction
					# cr = np.vstack(( fwdr, bwdr[::-1] )) # lump fwd and bwd together
					# Ymovingavgright = np.mean( cr, axis=0 ) # average
					Ymeansubleft=Ycutleft-Ymovingavgleft
					Ymeansubright=Ycutright-Ymovingavgright
					stoppoint=len(Ymeansubleft)-cut
					Ybeforemeansub=np.concatenate((Ycutleft[cut:stoppoint],Ycutright[cut:stoppoint]))
					Ymeansub=np.concatenate((Ymeansubleft[cut:stoppoint],Ymeansubright[cut:stoppoint]))
					Ymovingavg=np.concatenate((Ymovingavgleft[cut:stoppoint],Ymovingavgright[cut:stoppoint]))
					Yerr=[stderror(N, val) for val in Ybeforemeansub]
					summand+=Ymeansub
					before+=Ybeforemeansub
					variance=np.square(Yerr)
					smooth+=Ymovingavg
				summand/=stopindex
				before/=stopindex
				variance/=stopindex
				smooth/=stopindex
				Xbband=np.concatenate((Xcutleft[cut:stoppoint],Xcutright[cut:stoppoint]))
				centerfreq=(LOFreq+4.09)*np.power(10,9)
				Xtruefreq=Xbband+centerfreq*np.ones(len(Xbband))
				#errorbar(Xbband/np.power(10,6),summand,yerr=np.sqrt(variance),fmt='o')
				#plot(bband,power,'b',label="readTMCavity.C")
				plot(Xbband[93:]/np.power(10,6),summand[93:])
				#plot(freq[118:]/np.power(10,6), poweravg5[118:],'g',label="mathematica script window_size=5")
				xlabel('Freq Offset (MHz) from %s GHz' % str(LOFreq+4.09))
				ylabel('Mean Subtracted Power (mW/Hz)')
				title(timestamp)
				#legend(loc=8)
				grid(True)
				savefig(pltdir+'/'+timestamp+'avgsubplt_'+str(stopindex)+'.png',bbox_inches='tight')
				close()
				errorbar(Xbband[93:]/np.power(10,6),summand[93:],yerr=np.sqrt(variance)[93:],fmt='o')
				xlabel('Freq Offset (MHz) from %s GHz' % str(LOFreq+4.09))
				ylabel('Mean Subtracted Power (mW/Hz)')
				grid(True)
				title(timestamp)
				savefig(pltdir+'/'+timestamp+'avgsubpltposfreq_'+str(stopindex)+'.png',bbox_inches='tight')
				close()
				plot(Xbband[93:]/np.power(10,6),before[93:],'bs')
				plot(Xbband[93:]/np.power(10,6), smooth[93:],'r^')
				#legend(loc=10)
				xlabel('Freq Offset (MHz) from %s GHz' % str(LOFreq+4.09))
				ylabel('Input Power (mW/Hz)')
				grid(True)
				title(timestamp)
				savefig(pltdir+'/'+timestamp+'inputpltposfreq_'+str(stopindex)+'.png',bbox_inches='tight')
				close()
				outfile=open(filedir+'/'+timestamp+'avgsubfile_'+str(stopindex)+'.txt','w')
				for j in range(len(Xtruefreq)):
					print>>outfile, Xbband[j], Xtruefreq[j], summand[j], smooth[j], np.sqrt(variance)[j]
				outfile.close()
			except KeyError:
				print "timestamp %s not found in spreadsheet" % timestamp
				#print>>notfound, timestamp
			except ValueError:
				print "probably no Cavity Freq entered in spreadsheet for timestamp %s" % timestamp

#notfound.close()
#reset current dir to file's working dir
os.chdir('E:/Documents/Documents/PythonScripts')
map(f,daylist)
os.chdir('E:/Documents/Documents/PythonScripts')

