import numpy as np
import pylab as pl

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

def make_rect_grid(nx,ny,x_spac,y_spac):
    x = np.arange(nx)
    mn = np.mean(x)
    x = x - mn # What the FUCK?
    x *= x_spac
    print x
    y = np.arange(ny)
    mn = np.mean(y)
    y = y - mn
    y *= y_spac
    print y
    pos = {'x':[],'y':[]}
    for i in range(ny):
        for ix in x:
            pos['x'].append(ix)
        for j in range(nx):
            pos['y'].append(y[i])
    return pos

# the HERA-37 array
spac = 14.4
diam = 14.
# 7 -> 37, 21 -> 331
mn = 7
pos = makehex(mn,spac)
#for i in range(hexno(mn)):
#    print '%7.3f %7.3f'%(pos['x'][i],pos['y'][i])

pos_psa32 = make_rect_grid(8,4,30,4)
pos_test = make_rect_grid(4,4,3.5,3.5)

for i in range(hexno(mn)):
    dish = circle(pos['x'][i],pos['y'][i],diam/2.)
    pl.plot(dish['x'],dish['y'],'b')

pl.plot(pos_psa32['x'],pos_psa32['y'],'ro')
pl.plot(pos_test['x'],pos_test['y'],'go')

#pl.xlim(np.array([-1.2,1.2])*np.array(pos['x']).max())
#pl.xlim(np.array([-1.2,1.2])*np.array(pos['y']).max())
pl.axis('equal')
pl.xlabel('East-West Antenna Position [m]')
pl.ylabel('North-South Antenna Position [m]')
pl.savefig('/home/jaguirre/PAPER/AAS2014/hera37antenna_positions.png')
pl.close()


