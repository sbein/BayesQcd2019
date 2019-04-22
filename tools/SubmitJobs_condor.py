import os, sys
from glob import glob
from random import shuffle 
from time import sleep

test = False

import argparse
parser = argparse.ArgumentParser()
defaultInfile_ = "../SampleProduction/delphes/delphes_gjet4.root"
parser.add_argument("-analyzer", "--analyzer", type=str,default='tools/ResponseMaker.py',help="analyzer")
parser.add_argument("-v", "--verbosity", type=int, default=1,help="analyzer script to batch")
parser.add_argument("-printevery", "--printevery", type=int, default=100000,help="short run")
parser.add_argument("-fin", "--fnamekeyword", type=str,default=defaultInfile_,help="file")
parser.add_argument("-bootstrap", "--bootstrap", type=str, default='0',help="boot strapping (0,1of5,2of5,3of5,...)")
parser.add_argument("-quickrun", "--quickrun", type=bool, default=False,help="short run")
args = parser.parse_args()
fnamekeyword = args.fnamekeyword
inputFiles = glob(fnamekeyword)
bootstrap = args.bootstrap
verbosity = args.verbosity
printevery = args.printevery
quickrun = args.quickrun
analyzer = args.analyzer.replace('python/','').replace('tools/','')
    

try: 
	moreargs = ' '.join(sys.argv)
	moreargs = moreargs.split('--fnamekeyword')[-1]
	moreargs = ' '.join(moreargs.split()[1:])
except: 
	moreargs = ''
	
moreargs = moreargs.strip()
print 'moreargs', moreargs

cwd = os.getcwd()
shuffle(inputFiles)


filesperjob = 1

def main():
    ijob = 1
    files = ''
    for ifname, fname in enumerate(inputFiles):
        files += fname+','
        print fname
        if (ifname+1)%filesperjob==filesperjob-1:
            print '==='*3
            jobname = analyzer.replace('.py','')+'-'+fname[fname.rfind('/')+1:].replace('.root','_'+str(ijob))
            print 'moreargs.split()', moreargs.split()
            if len(moreargs.split())>0: 
            	print 'trying to beef up jobname', jobname
            	jobname = jobname+moreargs.replace(' ','-')
            	print 'tried to beef up jobname', jobname
            print 'jobname', jobname
            fjob = open('jobs/'+jobname+'.sh','w')
            files = files[:-1]
            fjob.write(jobscript.replace('CWD',cwd).replace('FNAMEKEYWORD',fname).replace('ANALYZER',analyzer).replace('MOREARGS',moreargs).replace('JOBNAME',jobname))
            fjob.close()
            os.chdir('jobs')
            command = 'condor_qsub -cwd '+jobname+'.sh &'
            print command
            if not test: os.system(command)
            os.chdir('..')
            files = ''
            ijob+=1
            #if test: break
            sleep(0.05)
        
jobscript = '''#!/bin/zsh
source /etc/profile.d/modules.sh
source /afs/desy.de/user/b/beinsam/.bash_profile
module use -a /afs/desy.de/group/cms/modulefiles/
module load cmssw
export THISDIR=$PWD
echo "$QUEUE $JOB $HOST"
source /afs/desy.de/user/b/beinsam/.bash_profile
cd CWD
cmsenv
cd $THISDIR
export timestamp=$(date +%Y%m%d_%H%M%S%N)
mkdir $timestamp
cd $timestamp
cp -r CWD/tools .
cp -r CWD/usefulthings .
cp -r CWD/src .
python tools/ANALYZER --fnamekeyword FNAMEKEYWORD MOREARGS > CWD/jobs/JOBNAME.out > CWD/jobs/JOBNAME.out 2> CWD/jobs/JOBNAME.err
mv *.root CWD/output/smallchunks
cd ../
rm -rf $timestamp
'''

main()
