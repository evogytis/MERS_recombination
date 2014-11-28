import glob
import os
import re
path='./'

for filename in glob.glob(os.path.join(path,'*.arp')):
    linecounter=0
    fastafile=open(filename.replace('.arp','.fasta'),'w')
    locfile=open(filename.replace('.arp','_locs.txt'),'w')
    for line in open('%s'%(filename),'r'):
        if linecounter<18:
            pass
        elif linecounter==19:
            sites=map(float,line.strip('\n')[1:].split(','))
            print>>locfile,'%s %.3f L'%(len(sites),30000)
            print>>locfile,'\n'.join(['%.3f'%(x) for x in sites])
            print>>fastafile,'109 %s 1'%(len(sites))
        elif 'Sample' in line or len(line)==1 or '}' in line or 'Structure' in line or 'Group' in line:
            pass
        else:
            cerberus=re.search('([0-9\_]+)\t1\t([ACTG]+)',line)
            if cerberus is not None:
                #pass
                #print cerberus.group(1),cerberus.group(2)
                print>>fastafile,'>%s\n%s'%(cerberus.group(1),cerberus.group(2))
            #print line[:100]
        linecounter+=1
    fastafile.close()
    locfile.close()

print 'Done!'
