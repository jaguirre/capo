# Calculations relevant for the HERA data volume

import numpy as np
import pylab as pl
import scipy as sci
from scipy import signal

def data_volume(t_int=10,f_day=0.5,n_days=180,n_chan=3072,n_ant=331,n_pol=2,bytes_per_complex=8):
    # Constants
    sec_per_day = 24*60*60
    bits_per_byte = 8
    TB_per_byte = 1.e-12
    Mbits_per_bit = 1e-6
    
    # Intermediate quantities
    n_samp_per_day = sec_per_day/t_int*f_day
    n_corr_in = n_ant * n_pol
    # Count auto-correlations as baselines
    n_bl = n_corr_in * (n_corr_in + 1)/2.
    bytes_per_day = n_chan * n_bl * n_samp_per_day * bytes_per_complex

    # Total daily data volume
    TB_per_day = bytes_per_day * TB_per_byte
    # Total season data volume
    TB_per_season = TB_per_day * n_days
    # This is the rate at which data spews from the correlator
    Mbits_per_sec = bytes_per_day * bits_per_byte * Mbits_per_bit / (sec_per_day * f_day) 
    
    return {'TB_per_day':TB_per_day,'Mbits_per_sec':Mbits_per_sec}
    
    # It is useful to calculate data rates, because this tells you how fast you
    # have to write tape or copy data to keep up
                
# the magic hex numbers
def hexno(mid_num):
    num = 0
    for i in range(mid_num-1,mid_num/2,-1):
#        print i
        num += i
    num *= 2
    num += mid_num
    return num

def makehex(mid_num,spacing):
    pos = {'x':[],'y':[]}
    dy = np.cos(30.*np.pi/180.)
    # Walk along the center
    for x in rng(mid_num):
        pos['x'].append(x*spacing)
        pos['y'].append(0)
    # Walk up and down 
    for row,num_row in enumerate(range(mid_num-1,mid_num/2,-1)):
#        print row
#        print num_row
        for x in rng(num_row):
            pos['x'].append(x*spacing)
            pos['y'].append((row+1)*dy*spacing)
            pos['x'].append(x*spacing)
            pos['y'].append(-(row+1)*dy*spacing)
    #pos = np.array(pos)
    return pos

def rng(num):
    return np.arange(num)-np.mean(np.arange(num))

def circle(x0,y0,r):
    an = np.linspace(0,2*np.pi,100)
    coord ={'x':np.array(r*np.cos(an)+x0),'y':np.array(r*np.sin(an)+y0)}
    return coord

def bin_center(bin_edges):
    # Why the fuck doesn't hist do this?
    nb = len(bin_edges)
    bc = np.zeros(nb-1)
    for i in range(1,nb):
        bc[i-1] = (bin_edges[i]+bin_edges[i-1])/2
    return bc

def xy2pix(x,y,xrange,yrange,nx,ny):
    dx = float(np.diff(xrange)[0])/(nx-1)
    dy = float(np.diff(yrange)[0])/(ny-1)
    print dx, dy
    # This centers the x,y values in the resulting 2-D image
    i = np.array((x-np.average(x))/dx+nx/2,dtype='int64')
    j = np.array((y-np.average(y))/dy+ny/2,dtype='int64')
    return (i,j)

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)    
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)               
    #print size, sizey    
    x, y = sci.mgrid[-size:size+1, -sizey:sizey+1]
    g = sci.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc)    

def make_hera_array(spac=14.6,diam=14.,mn=21):
    """ Define the HERA hexagonal array.  Default parameters are 14.6 meter spacing, 14 meter dishes, and 21 antennas on the center line.  """
    
    # the HERA-331 array
    spac = 14.6
    diam = 14.
    # Mapping from the number of antennas along the center line to the total 
    # number
    # 7 -> 37, 21 -> 331
    mn = 21
# Calculate positions
pos = makehex(mn,spac)
n_ant = len(pos['x'])
# Calculate baselines

# ---
# This chunk was implemented for speed, but we don't really need the speed, and brute force makes the results easier to classify

bl = {'x':np.zeros([n_ant,n_ant]),'y':np.zeros([n_ant,n_ant]),'len':np.zeros([n_ant,n_ant])}
# Avoid looping over antennas
for coord in ['x','y']:
    bl[coord] = np.outer(pos[coord],np.ones(n_ant))-np.outer(np.ones(n_ant),pos[coord])

# But now, how to classify by redundancy type efficiently?

#bl['len'] = np.sqrt(np.power(bl['x'],2)+np.power(bl['y'],2))

# Just for giggles, let's plot up the distribution of lengths
#h,b = np.histogram(bl['len'],bins=10000,range=[0.1,bl['len'].max()])
#bc = bin_center(b)
#pl.bar(bc,h)
#pl.show()
# ---

# I guess there's nothing better than brute force for identifying
# redundant baseline types ...
nbl = n_ant*(n_ant-1)/2
baselines = {'x':np.zeros(nbl),'y':np.zeros(nbl),'len':np.zeros(nbl),'ang':np.zeros(nbl),'type':np.zeros(nbl)}

ibl = 0
for i in range(n_ant):
    for j in range(i+1,n_ant): # print i,j
        
        # Take the convention that the ordering of antennas is such that
        # x_bl = x_i - x_j where x_i > x_j

        # You get a different answer if you replace y with x.  I had originally
        # chosen x so that all baselines would have the property that the sense
        # of the source delay is the same for all of them (which depends on the
        # E-W component)
        if pos['y'][i] >= pos['y'][j]:
            ii = i
            jj = j
        else:
            ii = j
            jj = i
        baselines['x'][ibl] = pos['x'][ii] - pos['x'][jj]
        baselines['y'][ibl] = pos['y'][ii] - pos['y'][jj]
        baselines['len'][ibl] = np.sqrt(np.power(baselines['x'][ibl],2) + np.power(baselines['y'][ibl],2))
        baselines['ang'][ibl] = np.arctan2(baselines['y'][ibl],baselines['x'][ibl])
        ibl+=1

if False:
    pl.plot(baselines['x'],baselines['y'],'bo')
    pl.axis('equal')
    pl.show()

# Now, let's identify the unique ones
bltype = 1
tol = 0.001 # Tolerance for matching baseline lengths, in meters.
while baselines['type'].min() == 0:
    # Find the baselines which haven't been typed yet
    which_bl = (np.where(baselines['type'] == 0))[0]
    srt = np.argsort(baselines['len'][which_bl])
    which_bl = which_bl[srt[0]]
    which_match = np.where((baselines['x'] > baselines['x'][which_bl]-tol) *
                           (baselines['x'] < baselines['x'][which_bl]+tol) *
                           (baselines['y'] > baselines['y'][which_bl]-tol) *
                           (baselines['y'] < baselines['y'][which_bl]+tol))[0]
    baselines['type'][which_match] = bltype
    bltype+=1

ntypes = baselines['type'].max()
    
baseline_types = {'ID':np.zeros(ntypes),'len':np.zeros(ntypes),'ang':np.zeros(ntypes),'x':np.zeros(ntypes),'y':np.zeros(ntypes),'mult':np.zeros(ntypes)}
    
for blt in range(1,int(baselines['type'].max())+1):
    indx = blt-1
    examples = np.where(baselines['type'] == blt)[0] # God, numpy where sucks
    example = examples[0]

    baseline_types['ID'][indx] = baselines['type'][example]
    baseline_types['len'][indx] = baselines['len'][example]
    baseline_types['ang'][indx] = np.degrees(baselines['ang'][example])
    baseline_types['x'][indx] = baselines['x'][example]
    baseline_types['y'][indx] = baselines['y'][example]
    baseline_types['mult'][indx] = len(examples)

    print '%3d %4d %.3f %.3f'%(baseline_types['ID'][indx], baseline_types['mult'][indx], baseline_types['len'][indx], baseline_types['ang'][indx])

print
print 'Number of distinct types:       %3d'%(ntypes)
print 'Number of distinct lengths:     %3d'%len(set(baseline_types['len']))
print 'Number of distinct angles:      %3d'%len(set(baseline_types['ang']))
print 'Number of distinct E-W lengths: %3d'%len(set(baseline_types['x']))

#pl.bar(baseline_types['ID'],baseline_types['mult'])

# Now need to calclate the fringe rate for each baseline.  For this, I
# need the fully 3-D baseline vector.  I think we have been defining
# the fringe rate as 

npix_uv = 1024
uvplane = np.zeros([npix_uv,npix_uv])
i,j = xy2pix(bl['x'],bl['y'],[bl['x'].min()*1.1,bl['x'].max()*1.1],[bl['y'].min()*1.1,bl['y'].max()*1.1],npix_uv,npix_uv)

i = i.flatten()
j = j.flatten()

for indx in range(len(i)):
    # Stupid MATLAB image convention requires flipped i,j to get it to
    # plot which imshow as x,y
    uvplane[j[indx],i[indx]] += 1

if False:
    # Calculation is slow
    uvplane_conv = blur_image(uvplane,10)
    
if False:
    pl.imshow(uvplane_conv,aspect='auto',interpolation='nearest')
    pl.axis('equal')
    pl.show()
    
if False:
    # plot up the HERA array
    for i in range(hexno(mn)):
        dish = circle(pos['x'][i],pos['y'][i],diam/2.)
        pl.plot(dish['x'],dish['y'],'b')

    pl.axis('equal')
    pl.xlabel('East-West Antenna Position [m]')
    pl.ylabel('North-South Antenna Position [m]')
    pl.show()
    #pl.savefig('/home/jaguirre/PAPER/AAS2014/hera37antenna_positions.png')
    pl.close()


