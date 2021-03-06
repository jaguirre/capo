#! /usr/bin/env python

import aipy as a, sys, optparse, pylab as pb, numpy as np, string

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
o.add_option('--clean', dest='clean', type='float',
    help='Deconvolve delay-domain data by the "beam response" that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
o.add_option('-o', '--outfile', dest='outfile', type='string',
    help='Filename of matrix.npy file.')

opts,args = o.parse_args(sys.argv[1:])

matrix = np.zeros((33,33))
dout = {}
for filename in args:
    print "Reading", filename
    uv = a.miriad.UV(filename)
    uv.select('auto', -1, -1, include=False)
    for p, d, in uv.all():
        (uvw, t, (i, j)) = p
        bl = '%d,%d' % (i,j)
        dflags = np.logical_not(d.mask).astype(np.float)
        dgain = np.sqrt(np.average(dflags**2))
        dker = np.fft.ifft(dflags)
        d = d.filled(0)
        d = np.fft.ifft(d)
        if not opts.clean is None and not np.all(d == 0):
            d, info = a.deconv.clean(d, dker, tol=opts.clean)
            d += info['res'] / dgain
        d = np.ma.array(d)
        d = np.ma.absolute(d.filled(0))
        d = np.ma.concatenate([d[d.shape[0]/2:],
                           d[:d.shape[0]/2]], axis=0)
        d.shape = (1,) + d.shape
        if not dout.has_key(bl): dout[bl] = []
        dout[bl].append(d)
    del(filename)


bls = dout.keys()
for cnt, bl in enumerate(bls):
    blstring = string.split(bl, ',')
    l = float(blstring[0])
    m = float(blstring[1])
    #print np.sum(np.array(dout[bl])[:,:,((uv['nchan']/2) - 10):((uv['nchan']/2) + 10)])
    temp = np.array(dout[bl])
    temp.shape = (np.array(dout[bl]).shape[0],np.array(dout[bl]).shape[2])
    #if cnt == 1:
    #    print bl
    #    pb.imshow((temp[:,uv['nchan']/2 - 10 : uv['nchan']/2 + 10]))
    #    pb.show()
    matrix[l,m] = np.sum(np.array(dout[bl])[:,:,uv['nchan']/2 - 5 : uv['nchan']/2 + 5])
    
    #print np.array(dout[bl]).shape

np.save(opts.outfile, matrix)

#pb.imshow(matrix, aspect='auto', interpolation='nearest')#, vmax=8, vmin=5)
#pb.colorbar()
#pb.show()
