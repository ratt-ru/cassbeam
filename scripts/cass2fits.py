#!/usr/bin/env python
"""
Convert cassbeam output Jones dat files into FITS beams for meqtrees

Assumes cassbeam output files contain freq information in the filename of the form <prefix>-<freq>MHz.jones.dat
"""

import sys,os
import datetime
import pyfits as pf
import numpy as n

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] <CASSBEAM DAT FILES>')
    o.set_description(__doc__)
    o.add_option('-c','--clobber',dest='clobber',action='store_true',
        help='Clobber/overwrite output FITS files if they exist')
    o.add_option('-o','--output',dest='output',default='cassBeam_',
        help='Output FITS filename prefix')
    o.add_option('-p','--pixel',dest='pixel',type='float',default=0.01,
        help='Pixel scale factor in radians, this is the beampixelscale in the cassbeam output parameter file, default: 0.01')
    opts, args = o.parse_args(sys.argv[1:])

    freqs=[]
    beamCube=[]
    for fid,fn in enumerate(args):
        print 'Reading %s (%i of %i)'%(fn,fid+1,len(args))
        freq=float(fn.split('/')[-1].split('-')[1].split('MHz')[0]) #parse filename for frequency
        freqs.append(freq)
        data=n.fromfile(fn,dtype=float,sep=' ')
        data=data.reshape((data.shape[0]/8,8))
        beamCube.append(data)
    freqs=n.array(freqs)*1e6
    beamCube=n.array(beamCube)

    dim=n.sqrt(beamCube.shape[1])
    beamCube=beamCube.reshape((beamCube.shape[0],dim,dim,8))

    #normalize cube so that peak value is 1.0, is this legit or required?
    beamCube/=n.max(n.abs(beamCube))
    
    #generate a fits beam header
    templateCube=n.zeros((1,beamCube.shape[0],dim,dim),dtype=float)
    hdu=pf.PrimaryHDU(templateCube)

    ctime=datetime.datetime.today()
    hdu.header.update('DATE','%s'%ctime)
    hdu.header.update('DATE-OBS','%s'%ctime)
    hdu.header.update('ORIGIN', 'GFOSTER')
    hdu.header.update('TELESCOP', 'VLA')
    hdu.header.update('OBJECT', 'beam')
    hdu.header.update('EQUINOX', 2000.0)

    print 'Using pixel scale factor: %f radians'%opts.pixel
    # note: defining M as the fastest moving axis (FITS uses
    # FORTRAN-style indexing) produces an image that when
    # viewed with kview / ds9 etc looks correct on the sky 
    # with M increasing to left and L increasing toward top
    # of displayed image
    if beamCube.shape[1]%2==0: crpixVal=int(beamCube.shape[1]/2)
    else: crpixVal=int(((beamCube.shape[1]-1)/2)+1)

    hdu.header.update('CTYPE1', 'M')
    hdu.header.update('CDELT1', (-1.0) * opts.pixel, 'in radians')
    hdu.header.update('CRPIX1', crpixVal, 'reference pixel (one relative)')
    hdu.header.update('CRVAL1', 0.0, 'M = 0 at beam peak')
    hdu.header.update('CTYPE2', 'L')
    hdu.header.update('CDELT2', opts.pixel, 'in radians')
    hdu.header.update('CRPIX2', crpixVal, 'reference pixel (one relative)')
    hdu.header.update('CRVAL2', 0.0, 'L = 0 at beam peak')

    #determine frequency step by assuming equal frequency steps and taking the difference of the first two frequencies
    sortFreqs=n.sort(freqs)
    if freqs.shape[0]>1: freqStep=(sortFreqs[-1]-sortFreqs[0])/freqs.shape[0];
    else: freqStep=1.
    hdu.header.update('CTYPE3', 'FREQ')
    hdu.header.update('CDELT3', freqStep, 'frequency step in Hz')
    hdu.header.update('CRPIX3', 1, 'reference frequency postion')
    hdu.header.update('CRVAL3', sortFreqs[0], 'reference frequency')
    hdu.header.update('CTYPE4', 'STOKES')
    hdu.header.update('CDELT4', 1) 
    hdu.header.update('CRPIX4', 1)
    hdu.header.update('CRVAL4', -5)
    
    # in case the frequency range is irregularly sampled, add keywords giving the actual frequency values
    for i,fq in enumerate(sortFreqs):
      hdu.header.update('GFREQ%d'%(i+1),fq);

    # create initial HDUList
    hdulist = pf.HDUList([hdu])

    #for XX,XY,YX,YY,real,imag: write a fits file
    for pid,pol in enumerate(['xx','xy','yx','yy']):
        for cid,cmplx in enumerate(['re','im']):
            #write data to FITS data
            hdu.data=beamCube[:,:,:,2*pid+cid]
            ofn=opts.output+pol+'_'+cmplx+'.fits'
            if opts.clobber or not os.path.isfile(ofn): hdulist.writeto(ofn,clobber=opts.clobber)
            else: print 'File: %s exists, and clobber parameter not set, skipping'%ofn

