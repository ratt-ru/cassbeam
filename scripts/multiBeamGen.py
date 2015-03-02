#!/usr/bin/env python

"""
Batch script to generate a multi-frequency FITS beam from a template cassBeam input file
"""

import sys,os,os.path
from subprocess import call
import numpy as n
import pyrap.tables as pt

#HARDCODE cassbeam and cass2fits.py locations
cbBin='/usr/local/bin/cassbeam'
c2fScript = os.path.dirname(__file__)+'/cass2fits.py'


if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] -t <CASSBEAM INPUT TEMPLATE>')
    o.set_description(__doc__)
    o.add_option('-c','--clobber',dest='clobber',action='store_true',
        help='Clobber/overwrite output FITS files if they exist')
    o.add_option('-f','--freq',dest='freq',default='1000,1200',
        help='Frequency range to generate beams at, can set to an MS and freqs \
        are pulled from the table or an a start,stop notation (MHz). Step sizes \
        are picked based on maintatining pixel resolution. default: 1000,1200')
    o.add_option('-o','--output',dest='output',default='cassBeam',
        help='Output FITS filename prefix')
    o.add_option('-p','--pixels',dest='pixels',default=128, type='int',
        help='Number of pixels in the output FITS file, default: 128')
    o.add_option('-t','--template',dest='template',default='vla-template.in',
        help='cassBeam input file to be used as the template, default %default')
    o.add_option('--pixelsperbeam',dest='ppb',default=32, type='int',
        help='Pixels across the primary lobe for the lowest frequency, default: 32')
    opts, args = o.parse_args(sys.argv[1:])
    
    if not opts.template:
      o.error("template must be specified with -t");
      

    #Sort out the frequencies, either they are input or from an MS
    if os.path.exists(opts.freq):
        print 'Getteing frequencies from %s'%opts.freq
        ms=pt.table(opts.freq+'/SPECTRAL_WINDOW',readonly=True)
        freqs=ms.getcol('CHAN_FREQ')[0]
        startFreq=freqs[0]
        stopFreq=freqs[-1]
    else:
        startFreq,stopFreq=map(lambda x:float(x)*1e+6,opts.freq.split(','))
    #compute the step freqs
    print 'Start Freq: %.4f GHz \t Stop Freq: %.4f GHz'%(startFreq/1.e9,stopFreq/1.e9)
    ppb0=opts.ppb
    freqs=[]
    b=0
    ppb=[]
    cFreq=startFreq #current freq, increment until cFreq>stopFreq
    while cFreq<stopFreq:
        #cFreq=startFreq+(startFreq*((ppb0/float(ppb0-b))-1.))
        cFreq=startFreq*(1.+((ppb0/float(ppb0-b))-1.))
        b+=1
        freqs.append(cFreq)
        ppb.append(ppb0-b)
    freqs=n.array(freqs)
    ppb=n.array(ppb)

    #convert freqs to GHz
    freqs/=1.e9
    print 'Generating slice at freqs (GHz):',freqs

    gridsize=2*(opts.pixels)

    #read template file
    fh=open(opts.template)
    templateLines=fh.readlines()
    fh.close()

    inputFiles=[]
    print 'Generating cassBeam input files'
    for f,pix in zip(freqs,ppb):
        outStr=opts.output+'-%.2fMHz'%(f*1.e3)
        #outStr=opts.output+'-%.2fMHz'+str(f*1.e3)+'MHz'
        outLines=[]
        for l in templateLines:
            if l.startswith('freq'): outLines.append('freq=%f\n'%f)
            elif l.startswith('gridsize'): outLines.append('gridsize=%i\n'%gridsize)
            elif l.startswith('out'): outLines.append('out=%s\n'%outStr)
            elif l.startswith('pixelsperbeam'): outLines.append('pixelsperbeam=%i\n'%pix)
            else: outLines.append(l)
        cassBeamFile=outStr
        inputFiles.append(cassBeamFile)
        fh=open(cassBeamFile+'.in','w')
        fh.writelines(outLines)
        fh.close()

    #run cassbeam
    for fid,fn in enumerate(inputFiles):
        print 'Running cassbeam with %s.in (%i of %i)'%(fn,fid+1,len(inputFiles))
        print [cbBin,fn+'.in']
        rcode=call([cbBin,fn+'.in'])

    #get the beampixelscale from the lowest frequency output param file
    fh=open(opts.output+'-%.2fMHz.params'%(freqs[0]*1.e3))
    #fh=open(opts.output+'-'+str(freqs[0]*1.e3)+'MHz.params')
    paramLines=fh.readlines()
    for l in paramLines:
        if l.startswith('beampixelscale'): bps=float(l.split('=')[-1])
    fh.close()

    #run cass2fits.py
    print 'Running %s'%c2fScript
    c2fArgs=' '
    if opts.clobber: c2fArgs+='-c '
    c2fArgs+='-p %f '%bps
    c2fArgs+='-o '+opts.output+'_ '
    for fn in inputFiles: c2fArgs+=fn+'.jones.dat '
    os.system(c2fScript+' '+c2fArgs)

