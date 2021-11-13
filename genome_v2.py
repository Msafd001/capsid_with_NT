param = {'r_ver': 0.3,'r_wall':0.2,'r_terminal':0.05,'r_gen':0.05}
from hoomd import md
import pandas as pd
import hoomd
import itertools
import numpy
import math
import matplotlib.cm
from matplotlib.colors import Normalize
import numpy as np
from random import seed
from random import random
from hoomd import _hoomd
from hoomd.md import _md
from hoomd.md import force
d1 = param['r_gen']*2
N_b = 368
ddd = str(param['r_gen']).replace('.','_')
def distance(x1,y1,z1,x2,y2,z2):
    return math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    
def radius_gyration(snap):
    
    size = len(snap.particles.position)
    
    x_cm = 0
    y_cm = 0
    z_cm = 0
    
    for i in range(size):
        x_cm = x_cm + snap.particles.position[i][0]
        y_cm = y_cm + snap.particles.position[i][1]
        z_cm = z_cm + snap.particles.position[i][2]
        
    x_cm = x_cm / size
    y_cm = y_cm / size
    z_cm = z_cm / size
    
    g = 0
    for i in range(size):
        g = g + (snap.particles.position[i][0] - x_cm) ** 2
        g = g + (snap.particles.position[i][1] - y_cm) ** 2
        g = g + (snap.particles.position[i][2] - z_cm) ** 2
    
    g = g /size
    return g
    print("radius of gyration = ",g)

def end_to_end_length(snap):
    
    size = len(snap.particles.position)
    
    x1 = snap.particles.position[0][0]
    y1 = snap.particles.position[0][1]
    z1 = snap.particles.position[0][2]
    x2 = snap.particles.position[size-1][0]
    y2 = snap.particles.position[size-1][1]
    z2 = snap.particles.position[size-1][2]
    return(distance(x1,y1,z1,x2,y2,z2))
    print("end to end length = ",distance(x1,y1,z1,x2,y2,z2))
    
    
# initialize position of particles
energy_t = []
ene_to = []
rad_g = []
for j in range(1):
    
    d = d1
    N = N_b
    #     step = math.floor(number_of_bids**(1./3)+0.1)  
    #     print("each size : ", step)
    nn = 0
    #     r = ((step*step)*d ) / (2*math.pi)
    r = math.sqrt((N*d*d ) / (4*math.pi))
    N_z = math.floor(math.pi*r/(d))
    for i in range(N_z) :
        r2 = r*math.sin(((i+1)/N_z)*math.pi)
        N_r = math.floor((2*math.pi*r2)/d)
        nn += N_r
    hoomd.device.GPU()
    sim = hoomd.Simulation(device=hoomd.device.GPU(),seed=1)
    # hoomd.context.initialize("--mode=cpu")
    number_of_bids = nn
    bond_size = number_of_bids - 1
    BOX = hoomd.Box( Lx=50, Ly=50, Lz=50)

    snap = hoomd.Snapshot()
    snap.particles.N = number_of_bids
    snap.configuration.box = BOX#[10, 10, 10, 0, 0, 0]
    # for i in range(number_of_bids):
    #     snap.particles.position[i] = random()*10-5
    snap.particles.types = ['P']
    # snap.particles.typeid[:] = [0]

    # initial_condistion = hoomd.data.make_snapshot(N = number_of_bids,
    #         box=hoomd.data.boxdim(Lx=50,Ly=50,Lz=50),
    #         particle_types=['P'],
    #         bond_types=['genome'])
    
    d = d1
    step = math.floor(number_of_bids**(1./3)+0.1)  
    
    d = d1
    N = N_b
    
    nn = 0
    r = math.sqrt((N*d*d ) / (4*math.pi))
    N_z = math.floor(math.pi*r/(d))
    for i in range(N_z) :
        r2 = r*math.sin(((i+1)/N_z)*math.pi)
        N_r = math.floor((2*math.pi*r2)/d)
        nn += N_r

    N_N = 0
    for i in range(N_z) :
        r2 = r*math.sin(((i+1)/N_z)*math.pi)
        N_r = math.floor((2*math.pi*r2)/(d*1))
        nn += N_r
        for j in range(N_r):
            snap.particles.position[N_N][0] = (r2 * math.cos((2*j*math.pi)/N_r)) * 1#((-1)**(j%2)) * (- step / 2 * d + k * d) + d
            snap.particles.position[N_N][1] = (r2 * math.sin((2*j*math.pi)/N_r)) * 1#((-1)**(i%2)) * (- step / 2 * d + j * d) + d
            snap.particles.position[N_N][2] = (r  * math.cos(((i+1)/N_z)*math.pi)) * 1
            snap.particles.diameter[N_N] = d + 0.000001

            N_N += 1
    print(f'number of particles are {N_N}')
    output_file_name = f'length{N_N}_{ddd}.gsd'
    snap.bonds.N = N_N -1  
    bond_gorup = []
    bond_type = []
    for i in range(snap.bonds.N):
        bond_gorup.append([i,i+1])
        bond_type.append('G')
        # snap.bonds.types[i] = 'G'
        # print( i , "  ", snap.bonds.types)
        # snap.bonds.group[i] =  [i,i+1]

    snap.bonds.types = bond_type
    snap.bonds.group[:] = bond_gorup
    
    for i in range(snap.bonds.N):
        snap.bonds.typeid[i] = 2

    sim.create_state_from_snapshot(snap)
    
    all = hoomd.filter.All()
    h1 = hoomd.md.bond.Harmonic()
    h1.params['G'] = dict( k = 1000.0, r0 = d1)#{'k':10000000,'r0':d1}#
    # nl = hoomd.md.nlist.Cell()
    # lj = hoomd.md.pair.LJ( nlist = nl , default_r_cut = d1)
    # lj.params[('P', 'P')] = { 'sigma': d1 * 1**(-1/6) , 'epsilon' : 0.1 }
    # lj.r_cut[('P', 'P')] = d1

    gsd_file = hoomd.write.GSD( trigger = hoomd.trigger.Periodic(50) , filename = output_file_name , filter = hoomd.filter.All() , mode = 'wb' )
    sim.operations.writers.append( gsd_file )

    # integrator.methods.append(nve)
    langevin = hoomd.md.methods.Langevin ( filter = hoomd.filter.All(), kT = 1.0, alpha = 0.0)
    integrator = hoomd.md.Integrator(dt=0.005, forces=[h1] , methods=[langevin])
    sim.operations.integrator = integrator
    sim.run(1e5)

    # fire = hoomd.md.minimize.FIRE(dt=0.05,forces=[h1,lj])
    # fire.force_tol = 1e-2
    # fire.energy_tol = 1e-5
    # fire.methods.append(hoomd.md.methods.NVE(hoomd.filter.All()))
    # sim.operations.integrator = fire
    # while not(fire.converged):
    #    sim.run(100)
