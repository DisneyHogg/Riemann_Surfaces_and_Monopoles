# Note this code came from a that initiall written by Paul Sutcliffe in March 2023. 
# Linden Disney-Hogg subsequently edited the work flow allowing for multiprocessing
# and spatial extents (x_vol) that are cuboidal. 

#############################################################################
# Practical imports
import numpy as np
from scipy import integrate
from numpy import linalg as LA
from scipy.linalg import null_space
import sys
from scipy import special
from scipy import pi

# Imports for multiprocessing and visualisation
from multiprocessing import Pool
from tqdm import *

# for reading command line arguments 
# see https://opensourceoptions.com/blog/how-to-pass-arguments-to-a-python-script-from-the-command-line/
import sys
import getopt

# This supposedly can be done a lot better with argparse. 
def get_args(argv):
    arg_nahmdata = ""
    arg_output = ""
    arg_nproc = "1"
    arg_xmax = "3"
    arg_ymax = "3"
    arg_zmax = "3"
    arg_spacing = "0.5"
    arg_help = "{0} -T <nahmdata> -p <nproc> -o <output> -x <xmax> -y <ymax> -z <zmax> -s <spacing>".format(argv[0])

    try:
        opts, args = getopt.getopt(argv[1:], "hT:p:o:x:y:z:s:", ["help", "nahmdata=", 
        "nproc=", "output=", "xmax=", "ymax=", "zmax=", "spacing="])
    except:
        print(arg_help)
        sys.exit(2)

    if not len(opts)==7:
        print(arg_help)
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-T", "--nahmdata"):
            arg_nahmdata = arg
        elif opt in ("-p", "--nproc"):
            arg_nproc = arg
        elif opt in ("-o", "--output"):
            arg_output = arg
        elif opt in ("-x", "--xmax"):
            arg_xmax = arg
        elif opt in ("-y", "--ymax"):
            arg_ymax = arg
        elif opt in ("-z", "--zmax"):
            arg_zmax = arg
        elif opt in ("-s", "--spacing"):
            arg_spacing = arg

    arg_output += 'xyz' + arg_xmax + arg_ymax + arg_zmax + 'h' + arg_spacing + 'args' + ''.join(args)

    return arg_nahmdata, np.int0(arg_nproc), arg_output, np.float128(arg_xmax), np.float128(arg_ymax), np.float128(arg_zmax), np.float128(arg_spacing), args

if __name__ == "__main__":
    nahmdata_file, nproc, name, xmax, ymax, zmax, h, args = get_args(sys.argv)
    print(nahmdata_file, nproc, name, xmax, ymax, zmax, h, args)

# Get the nahmdata function from the nahm data file
import os
import importlib
# Split of the path and filename from nahmdata_file
head, tail = os.path.split(nahmdata_file)
# add the path
sys.path.append(os.path.abspath(head))
# strip off the .py
tail = tail[:-3]

nahmdata_module = importlib.import_module(tail)
create_nahmdata = getattr(nahmdata_module, "create_nahmdata")
nahmdata = create_nahmdata(args)
# print(nahmdata(1.))

#############################################################################
####### These parameters may be changed by the user. ########################
# We keep these parameters hard coded, as the average user will not need these
dt = 0.01 # This is the lattice spacing in the t interval [0,2] (default is 0.01)
epsilon = 0.0001 # This is the distance to the pole for finding the residues
#############################################################################

# Have the creation of the x_vol in one cell to make it clear it is separate and an input of the function. 
from math import ceil
px = ceil(2*xmax/h)
py = ceil(2*ymax/h)
pz = ceil(2*zmax/h)
X, Y, Z = px*h/2, py*h/2, pz*h/2
x_vol = np.array([[[((2*ii/px-1)*X, (2*jj/py-1)*Y, (2*kk/pz-1)*Z) for kk in range(pz+1)]
                   for jj in range(py+1)] for ii in range(px+1)], dtype="f, f, f")

#############################################################################

n = nahmdata(1)[0].shape[0]
tau1 = np.array([[0j,1+0j],[1+0j,0j]])
tau2 = np.array([[0j,-1j],[1j,0j]])
tau3 = np.array([[1+0j,0j],[0j,-1+0j]])
idn = np.identity(n)
d = 2*n

def inprod(v1, v2, flag):
    if flag==1:
        values = [(teval[j]-1.0) * sum(np.conj(v1[i][j])*v2[i][j] 
                                       for i in range(np.shape(v1)[0]))
                  for j in range(np.shape(v1)[1])]
    else:
        values = [sum(np.conj(v1[i][j])*v2[i][j] for i in range(np.shape(v1)[0]))
                  for j in range(np.shape(v1)[1])]
    return integrate.simps(values, teval)

def ivshoot(leftright):
    t=epsilon
    if(leftright==1):
        t=2.0-epsilon
    T1,T2,T3=nahmdata(t)
    startv=1j*np.zeros((2*n,2*n))
    pole=-leftright*(-1j*np.kron(T1,tau1)-1j*np.kron(T2,tau2)-1j*np.kron(T3,tau3))*epsilon
    vals,vecs=LA.eig(pole)
    dimc=0
    for i in range(2*n):
        lam=np.real(vals[i])
        if(lam>0):
            startv[:,dimc]=vecs[:,i]
            dimc=dimc+1
    # Note that dimc should be n+1
    return dimc,startv

def higgs_at_ijk(ijk):
    return ijk, higgs_at_x(x_vol[ijk])

def higgs_at_x(xs):
    x1, x2, x3 = xs
    weylcon = np.kron(idn, x1*tau1 + x2*tau2 + x3*tau3)    
    def deriv(t, y):
        T1, T2, T3 = nahmdata(t)
        weyl = weylcon -1j*(np.kron(T1,tau1) + np.kron(T2,tau2) + np.kron(T3,tau3))
        return weyl@y    
    vleft = 1j*np.zeros((dimc1, 2*n, (np.shape(teval1)[0])))
    vright = 1j*np.zeros((dimc2, 2*n, (np.shape(teval2)[0])))
    vjoin = 1j*np.zeros((2*n, dimc1+dimc2))
    v = 1j*np.zeros((2, 2*n, np.shape(teval)[0]))
    # Shoot from the left
    for j in range(dimc1):
        y0 = lhs[:,j]
        soln = integrate.solve_ivp(deriv, t_span1, y0, t_eval=teval1)
        tlist = soln.t
        vv = soln.y
        vleft[j] = 1.0*vv        
    # Shoot from the right
    for j in range(dimc2):
        y0 = rhs[:,j]
        soln = integrate.solve_ivp(deriv, t_span2, y0, t_eval=teval2)
        tlist = soln.t
        vv = soln.y
        vright[j] = 1.0*vv
    # Combine left and right solutions
    for j in range(dimc1):
        vjoin[:,j] = vleft[j,:,-1]
    for j in range(dimc2):
        vjoin[:,j+dimc1] = vright[j,:,-1]
    # Make the solutions match in the middle
    ns = null_space(vjoin)
    if (np.shape(ns)[1]!=2):
        print("!!!!! ERROR : The kernel is not 2-dimensional")
        sys.exit(" ABORT ")
    ker1 = 1.0*ns[:,0]
    ker2 = 1.0*ns[:,1]
    for j in range(np.shape(teval1)[0]):
        v[0,:,j+1] += sum(ker1[jj]*vleft[jj,:,j] for jj in range(dimc1))
        v[1,:,j+1] += sum(ker2[jj]*vleft[jj,:,j] for jj in range(dimc1))      
            
    for j in range(1,np.shape(teval2)[0]):
        v[0,:,-(j+1)] -= sum(ker1[jj+dimc1]*vright[jj,:,j-1] for jj in range(dimc2))
        v[1,:,-(j+1)] -= sum(ker2[jj+dimc1]*vright[jj,:,j-1] for jj in range(dimc2))      
            
    # Calculate an orthonormal basis
    ip11 = inprod(v[0], v[0], 0)
    v1 = v[0]/np.sqrt(ip11)    
    ip12 = inprod(v1, v[1], 0)
    v2t = v[1] - v1*ip12
    ip22 = inprod(v2t, v2t, 0)
    v2 = v2t/np.sqrt(ip22)
    h11 = inprod(v1, v1, 1)
    h12 = inprod(v1, v2, 1)
    h22 = -h11
    h21 = np.conj(h12)
    h2 = 0.5*(h11**2+h22**2+2*h12*h21)
    return np.real(h2)

# Define the points inside the interval [0,2]
t_span1 = (dt, 1+dt)
teval1 = np.arange(dt, 1+dt, dt)
t_span2 = (2-dt, 1-dt)
teval2 = np.arange(2-dt, 1-dt, -dt)
teval = np.arange(0, 2+dt, dt)
dimc1, lhs = 1*ivshoot(-1)
dimc2, rhs = 1*ivshoot(1)

hi = np.zeros((px+1, py+1, pz+1))
N = (px+1)*(py+1)*(pz+1)
cs = max([px+1, py+1, pz+1])

if __name__ == '__main__':
    with Pool(processes=nproc) as pool:
        with tqdm(total=N) as pbar:
            # Here we are using something sneaky about numpy, namely that it has a function that 
            # given the shape of an array, returns an iterable over that array giving the indices
            # of the corresponding elements. This means we can loops over x_vol without having
            # to do anything smart about the iterable, and get an array hi that corresponds nicely
            # to x_vol.
            for iis, his in pool.imap_unordered(higgs_at_ijk, np.ndindex(x_vol.shape), chunksize=cs):
                pbar.update()
                hi[iis] = his

#############################################################################

with open('higgs_array_'+name, 'wb') as f:
    np.save(f, hi)

with open('x_array_'+name, 'wb') as f:
    np.save(f, x_vol)

