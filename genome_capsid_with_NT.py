param = {'r_ver': 0.1,'r_wall':0.2,'r_terminal':0.05,'r_gen':0.05}
from hoomd import md
import pandas as pd
import gsd.hoomd
import hoomd
import itertools
import numpy
import matplotlib.pyplot as plt
import math
import matplotlib.cm
from matplotlib.colors import Normalize
import numpy as np
from random import seed
from random import random
from hoomd import _hoomd
from hoomd.md import _md
from hoomd.md import force
import time
t0= time.time()

input_file_genome = 'length210_0_05.gsd'
v = 1
vertex_file_namge = '2_0vertexs.csv'
e1 = input_file_genome.replace('.','')
e2 = vertex_file_namge.replace('.','')
scale_2 = 10.0
end_terminal_length = 4
output_file_genome =  f"genoem_{e1[0:len(e1)-3]}_inside_r_capsid{e2[0:len(e2)-3]}_NT_Length_{end_terminal_length}_epsilon_is_{str(scale_2*0.1).replace('.','_')}.gsd"


def count_inside(snap):
    size_p = len(index_p)
    size_v = len(index_c)
    max_r = 0
    dis_1 = 0
    centeral_mass = [0,0,0]
    count_p = 0
    
    for i in range(size_v):
        centeral_mass[0] += snap.particles.position[index_c[i]][0]
        centeral_mass[1] += snap.particles.position[index_c[i]][1]
        centeral_mass[2] += snap.particles.position[index_c[i]][2]
        
    centeral_mass[0] = centeral_mass[0] / size_v
    centeral_mass[1] = centeral_mass[1] / size_v
    centeral_mass[2] = centeral_mass[2] / size_v

    for i in range(size_v):
        dis_1 = distance(centeral_mass[0],centeral_mass[1],centeral_mass[2],snap.particles.position[index_c[i]][0]
                 ,snap.particles.position[index_c[i]][1],snap.particles.position[index_c[i]][2])
        if (dis_1 > max_r):
            max_r = dis_1
            
    for j in range(size_p):
        dis_1 = distance(centeral_mass[0],centeral_mass[1],centeral_mass[2],snap.particles.position[index_p[j]][0]
             ,snap.particles.position[index_p[j]][1],snap.particles.position[index_p[j]][2])
        if (dis_1 < max_r ):
            count_p +=1 
    CGer = '\033[32m'
    CEND = '\033[0m'
    print(CGer+'number of polymer inside the capsid',count_p,CEND)

def density_fun(snap):
    size_p = len(index_p)
    size_v = len(index_c)
    delta_r = 0.04
    max_r = 0
    dis_1 = 0
    step = 0
    centeral_mass = [0,0,0]
    
    for i in range(size_v):
        centeral_mass[0] += snap.particles.position[index_c[i]][0]
        centeral_mass[1] += snap.particles.position[index_c[i]][1]
        centeral_mass[2] += snap.particles.position[index_c[i]][2]
        
    centeral_mass[0] = centeral_mass[0] / size_v
    centeral_mass[1] = centeral_mass[1] / size_v
    centeral_mass[2] = centeral_mass[2] / size_v

    for i in range(size_v):
        dis_1 = distance(centeral_mass[0],centeral_mass[1],centeral_mass[2],snap.particles.position[index_c[i]][0]
                 ,snap.particles.position[index_c[i]][1],snap.particles.position[index_c[i]][2])
        if (dis_1 > max_r):
            max_r = dis_1
    step = math.floor((max_r)/delta_r)
    density = np.zeros(step)
    rad = np.zeros(step)
    for i in range(step):
        for j in range(size_p):
            dis_1 = distance(centeral_mass[0],centeral_mass[1],centeral_mass[2],snap.particles.position[index_p[j]][0]
                 ,snap.particles.position[index_p[j]][1],snap.particles.position[index_p[j]][2])
            if (dis_1 > i *delta_r and dis_1<= (i+1)*delta_r):
                density[i] +=1
        rad[i] = (i+1)* delta_r
        density[i] = 3*density[i]/(4*math.pi*(math.pow((i+1)* delta_r,3)-math.pow(i* delta_r,3)))
           
    print('radius array :',len(rad))
    print('density array :',len(density))
    plt.plot(rad,density,'ro')
    plt.xlabel('distance')
    plt.ylabel('density')
    plt.show()

def distance(x1,y1,z1,x2,y2,z2):
    return math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

def update_wall_position(snap):
    triangle = pd.read_csv ('1_triangle.csv')
    vertex = pd.read_csv (vertex_file_namge)
    
    p1x = 0
    p1y = 0
    p1z = 0
    p2x = 0
    p2y = 0
    p2z = 0
    
    vx1 = vy1 = vz1 = vx2 = vy2 = vz2 = vx3 = vy3 = vz3 = 0 
    
    wall_x = 0
    wall_y = 0
    wall_z = 0
    wall_number = 0
    max_wall_number = len(index_w)

    diameter_of_m = 1 * param['r_wall']
    wall_vertex = pd.DataFrame(columns=['x','y','z','t'])
    
    for i in range(len(triangle)-1):
        vx1 = snap.particles.position[index_c[triangle.tv1[i]-1]][0]
        vy1 = snap.particles.position[index_c[triangle.tv1[i]-1]][1]
        vz1 = snap.particles.position[index_c[triangle.tv1[i]-1]][2]
        vx2 = snap.particles.position[index_c[triangle.tv2[i]-1]][0]
        vy2 = snap.particles.position[index_c[triangle.tv2[i]-1]][1]
        vz2 = snap.particles.position[index_c[triangle.tv2[i]-1]][2]

        number1 = distance(vertex.x[triangle.tv1[i]-1],vertex.y[triangle.tv1[i]-1],vertex.z[triangle.tv1[i]-1],
                          vertex.x[triangle.tv3[i]-1],vertex.y[triangle.tv3[i]-1],vertex.z[triangle.tv3[i]-1])/diameter_of_m
        
        for j in range(int(number1)):
            
            p1x = vertex.x[triangle.tv1[i]-1] + ( j / int(number1) ) * ( vertex.x[triangle.tv3[i]-1] - vertex.x[triangle.tv1[i]-1])
            p1y = vertex.y[triangle.tv1[i]-1] + ( j / int(number1) ) * ( vertex.y[triangle.tv3[i]-1] - vertex.y[triangle.tv1[i]-1])
            p1z = vertex.z[triangle.tv1[i]-1] + ( j / int(number1) ) * ( vertex.z[triangle.tv3[i]-1] - vertex.z[triangle.tv1[i]-1])

            p2x = vertex.x[triangle.tv2[i]-1] + ( j / int(number1) ) * ( vertex.x[triangle.tv3[i]-1] - vertex.x[triangle.tv2[i]-1])
            p2y = vertex.y[triangle.tv2[i]-1] + ( j / int(number1) ) * ( vertex.y[triangle.tv3[i]-1] - vertex.y[triangle.tv2[i]-1])
            p2z = vertex.z[triangle.tv2[i]-1] + ( j / int(number1) ) * ( vertex.z[triangle.tv3[i]-1] - vertex.z[triangle.tv2[i]-1])

            number2 = distance(p1x,p1y,p1z,p2x,p2y,p2z) / diameter_of_m
    
            vx1 = snap.particles.position[index_c[triangle.tv1[i]-1]][0]
            vy1 = snap.particles.position[index_c[triangle.tv1[i]-1]][1]
            vz1 = snap.particles.position[index_c[triangle.tv1[i]-1]][2]

            vx2 = snap.particles.position[index_c[triangle.tv2[i]-1]][0]
            vy2 = snap.particles.position[index_c[triangle.tv2[i]-1]][1]
            vz2 = snap.particles.position[index_c[triangle.tv2[i]-1]][2]
            
            vx3 = snap.particles.position[index_c[triangle.tv3[i]-1]][0]
            vy3 = snap.particles.position[index_c[triangle.tv3[i]-1]][1]
            vz3 = snap.particles.position[index_c[triangle.tv3[i]-1]][2]
            
            p1x = vx1 + ( j / int(number1) ) * (vx3 - vx1)
            p1y = vy1 + ( j / int(number1) ) * (vy3 - vy1)
            p1z = vz1 + ( j / int(number1) ) * (vz3 - vz1)
            p2x = vx2 + ( j / int(number1) ) * (vx3 - vx2)
            p2y = vy2 + ( j / int(number1) ) * (vy3 - vy2)
            p2z = vz2 + ( j / int(number1) ) * (vz3 - vz2)
            
            for k in range(int(number2)+1):
                wall_x = p1x + (k/int(number2)) * (p2x - p1x)
                wall_y = p1y + (k/int(number2)) * (p2y - p1y)
                wall_z = p1z + (k/int(number2)) * (p2z - p1z)
                
                snap.particles.position[index_w[wall_number]][0] = wall_x
                snap.particles.position[index_w[wall_number]][1] = wall_y
                snap.particles.position[index_w[wall_number]][2] = wall_z 
                
                wall_number = wall_number + 1
    return snap

def move_to_origin(snap):
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
    
    for i in range(size):
        snap.particles.position[i][0] = snap.particles.position[i][0] - x_cm 
        snap.particles.position[i][1] = snap.particles.position[i][1] - y_cm 
        snap.particles.position[i][2] = snap.particles.position[i][2] - z_cm
        
def create_shell():
   
    vertex = pd.read_csv (vertex_file_namge)
    line = pd.read_csv ('1_lines.csv')
    triangle = pd.read_csv ('1_triangle.csv')
    wall_vertex = pd.DataFrame(columns=['x','y','z','t'])
    
    p1x = 0
    p1y = 0
    p1z = 0
    p2x = 0
    p2y = 0
    p2z = 0
    wall_x = 0
    wall_y = 0
    wall_z = 0
    number1 = 0
    number2 = 0
    centeral_mass = [0,0,0]
    
    for i in range(len(vertex)):
        centeral_mass[0] += vertex.x[i]
        centeral_mass[1] += vertex.y[i]
        centeral_mass[2] += vertex.z[i]
        
    centeral_mass[0] = centeral_mass[0] / len(vertex)
    centeral_mass[1] = centeral_mass[1] / len(vertex)
    centeral_mass[2] = centeral_mass[2] / len(vertex)
    
    diameter_of_m = 1 * param['r_wall']
    wall_vertex = pd.DataFrame(columns=['x','y','z','t'])
    for i in range(len(triangle)-1):
        number1 = distance(vertex.x[triangle.tv1[i]-1],vertex.y[triangle.tv1[i]-1],vertex.z[triangle.tv1[i]-1],
                          vertex.x[triangle.tv3[i]-1],vertex.y[triangle.tv3[i]-1],vertex.z[triangle.tv3[i]-1])/diameter_of_m
        
        for j in range(int(number1)):
            
            p1x = vertex.x[triangle.tv1[i]-1] + ( j / int(number1) ) * ( vertex.x[triangle.tv3[i]-1] - vertex.x[triangle.tv1[i]-1])
            p1y = vertex.y[triangle.tv1[i]-1] + ( j / int(number1) ) * ( vertex.y[triangle.tv3[i]-1] - vertex.y[triangle.tv1[i]-1])
            p1z = vertex.z[triangle.tv1[i]-1] + ( j / int(number1) ) * ( vertex.z[triangle.tv3[i]-1] - vertex.z[triangle.tv1[i]-1])

            p2x = vertex.x[triangle.tv2[i]-1] + ( j / int(number1) ) * ( vertex.x[triangle.tv3[i]-1] - vertex.x[triangle.tv2[i]-1])
            p2y = vertex.y[triangle.tv2[i]-1] + ( j / int(number1) ) * ( vertex.y[triangle.tv3[i]-1] - vertex.y[triangle.tv2[i]-1])
            p2z = vertex.z[triangle.tv2[i]-1] + ( j / int(number1) ) * ( vertex.z[triangle.tv3[i]-1] - vertex.z[triangle.tv2[i]-1])

            number2 = distance(p1x,p1y,p1z,p2x,p2y,p2z) / diameter_of_m

            for k in range(int(number2)+1):
                wall_x = p1x + (k / int(number2)) * (p2x - p1x)
                wall_y = p1y + (k / int(number2)) * (p2y - p1y)
                wall_z = p1z + (k / int(number2)) * (p2z - p1z)
                wall_vertex.loc[len(wall_vertex.index)] = [wall_x,wall_y,wall_z,i]
    return shell(vertex,line,triangle,wall_vertex,centeral_mass)
class shell():
    """docstring for Shell"""
#     def __init__(self,param,vertex,line,triangle):
    def __init__(self,vertex,line,triangle,wall_vertex,centeral_mass):
        # super(polymer, self).__init__()
        
        self.rp = param['r_ver'] #param['Rx']
        self.R = 1.1 #param['R']
        
        self.l0 = 1.0 #param['l0']
        self.theta = 0.644#param['theta0']
        self.theta_upper = 40.0#param['theta_upper']
        self.theta_lower = 15.0#param['theta_lower']
        self.gamma = 120 #param['gamma']
        self.sigma = 0.0001 #param['sigma']
        self.seed = [1] #param['s_seed']
        self.typeid = [1] #param['s_typeid']
        self.bodyid = [-1,-2,-3] #param['s_bodyid']
        self.bondid = [1,2] #param['s_bondid']
        self.dihedralid = [0] #param['s_dihedralid']
        self.vertex = vertex
        self.centeral_mass = centeral_mass
        self.wall_vertex = wall_vertex
        self.line = line
        self.triangle = triangle
        self.n = len(self.vertex)
        self.n_wall = len(self.wall_vertex)
        self.particles = self.vertex[list('xyz')].values
        self.pids = [1]*self.n
        self.particles_wall = self.wall_vertex[list('xyz')].values
        self.pids_wall = [2] * self.n_wall
        self.pids_spik = [3] * (self.n * end_terminal_length)
        self.bonds = np.int_(self.line[['lv1','lv2']])
        self.bonds_spik = self.n * end_terminal_length
        self.dihedrals=[]
        
    def shellprint(self):
        cprint('Shell Structure:\n','yellow')
        
    def update_shell_info(self):
        self.n=len(self.vertex)
        self.update_pids()
        self.particles=np.array(self.vertex[list('xyz')])
        self.bonds=np.int_(self.line[['v0','v1']])
        self.update_dihedrals()


def initial2():
    
    # sim.create_state_from_gsd(filename = input_file_genome,frame = 249)
    f = gsd.hoomd.open(name=input_file_genome,mode='rb')
    snap = hoomd.Snapshot()
    # f = gsd.hoomd.Snapshot(input_file_genome)#, mode='rb'
    snap.particles.N = f[249].particles.N
    for i in range(f[249].particles.N):
        for j in range(3):
            snap.particles.position[i][j] = f[249].particles.position[i][j]
        snap.particles.typeid[i] = f[249].particles.typeid[i]
        snap.particles.body[i] = f[249].particles.body[i]
        snap.particles.diameter[i] = f[249].particles.diameter[i]

    snap.bonds.N = f[249].bonds.N
    
    bond_type = []
    for i in range(f[249].bonds.N):
        bond_type.append("G")
        snap.bonds.group[i] = f[249].bonds.group[i]
        snap.bonds.typeid[i] = f[249].bonds.typeid[i]

    snap.bonds.types = bond_type
    
    move_to_origin(snap)
    snap.dihedrals.types = ['capsomer']
    snap.particles.types = ['P','X','W','S'] #.append('X') ['P','X']
    snap.bonds.types  = ['G','shell','spik']
    snap.angles.types = ['ds']
    return snap

def assign_snap_particles(snap,obj,cover):
    # particles assign in range:[pstart,pend)
    if cover: 
        # pstart=int(np.nonzero(snap.particles.typeid==obj.typeid)[0][0])
        for i,pid in enumerate(snap.particles.typeid):
            if pid in obj.typeid: 
                pstart=i
                break
    else: pstart=snap.particles.N
    pend = pstart + obj.n + obj.n*end_terminal_length
    # snap.particles.N = (pend + obj.n_wall)
    snap.particles.N = (pend + obj.n_wall)
    for i, [cor,pid] in enumerate(zip(obj.particles,obj.pids),start=pstart):
        snap.particles.position[i] = cor
        snap.particles.typeid[i] = pid
        snap.particles.body[i] = obj.bodyid[0]
        snap.particles.diameter[i] = 2 * obj.rp
        
    count = 0
    
    bstart_bond = snap.bonds.N
    # snap.bonds.N = (bstart_bond )
    snap.bonds.N = (bstart_bond + end_terminal_length * obj.n)
    
    for i, [cor,pid] in enumerate(zip(obj.particles,obj.pids),start=pstart+obj.n):
        
        for j in range(end_terminal_length):
            snap.particles.diameter[pstart+obj.n + end_terminal_length*count + j] = 2*param['r_terminal']
        x = -cor[0] + obj.centeral_mass[0]
        y = -cor[1] + obj.centeral_mass[1]
        z = -cor[2] + obj.centeral_mass[2]  
#         print('radius of capsid : ', math.sqrt(x*x+y*y+z*z))
        for j in range(end_terminal_length):
            snap.particles.position[pstart+obj.n + end_terminal_length*count+j] = [ cor[0] + x*2*param['r_terminal']*(j)/math.sqrt(x**2+y**2+z**2)+x*(param['r_ver']+param['r_terminal']) / math.sqrt(x**2+y**2+z**2),
                                                                  cor[1] + y*2*param['r_terminal']*(j)/math.sqrt(x**2+y**2+z**2) + y*(param['r_ver']+param['r_terminal']) / math.sqrt(x**2+y**2+z**2),
                                                                  cor[2] + z*2*param['r_terminal']*(j)/math.sqrt(x**2+y**2+z**2) + z*(param['r_ver']+param['r_terminal']) / math.sqrt(x**2+y**2+z**2)]
        for j in range(end_terminal_length):
            snap.particles.typeid[pstart+obj.n + end_terminal_length*count+j] = obj.pids_spik[end_terminal_length*(i-pstart-obj.n)+j] #+ 3
        for j in range(end_terminal_length):
            snap.particles.body[pstart+obj.n + end_terminal_length*count+j] = obj.bodyid[2] #-3
            
        if (count==0):
            snap.bonds.group[bstart_bond+end_terminal_length*count] = [i-obj.n,i+end_terminal_length*count]
            snap.bonds.typeid[bstart_bond+end_terminal_length*count] = obj.bondid[1]
            for j in range(0,end_terminal_length-1):
                snap.bonds.group[bstart_bond+end_terminal_length*count+j+1] = [i+end_terminal_length*count+j,i+end_terminal_length*count+j+1]
                snap.bonds.typeid[bstart_bond+end_terminal_length*count+j+1] = obj.bondid[1]
        else:
            snap.bonds.group[bstart_bond+end_terminal_length*count] = [i-obj.n,i+(end_terminal_length-1)*count]
            snap.bonds.typeid[bstart_bond+end_terminal_length*count] = obj.bondid[1]
            for j in range(0,end_terminal_length-1):
                snap.bonds.group[bstart_bond+end_terminal_length*count+j+1] = [i+(end_terminal_length-1)*count+j,i+(end_terminal_length-1)*count+j+1]
                snap.bonds.typeid[bstart_bond+end_terminal_length*count+j+1] = obj.bondid[1]
        count += 1
    final_index = i

    for i,[cor_wall,pid_wall] in enumerate(zip(obj.particles_wall,obj.pids_wall),start=pstart+(end_terminal_length+1)*obj.n):
        snap.particles.diameter[i] = 0.4 * param['r_wall']
        snap.particles.position[i] = cor_wall
        snap.particles.typeid[i] = pid_wall
        snap.particles.body[i] = obj.bodyid[2]

    return pstart,pend
       
def assign_snap_bonds(snap,obj,pstart,cover):
    if cover: bstart = int(np.nonzero(snap.bonds.typeid==obj.bondid))#[0][0]
    else: bstart=snap.bonds.N

    # snap.bonds.resize(bstart+len(obj.bonds))
    snap.bonds.N = (bstart+len(obj.bonds))

    for i, bond in enumerate(obj.bonds,start=bstart):
        snap.bonds.group[i] = pstart-1+bond#pstart-1+bond
        snap.bonds.typeid[i] = obj.bondid[0]
        
        
def assign_snap_dihedrals(snap,obj,pstart,cover):
    if cover: dstart = int(np.nonzero(snap.dihedrals.typeid==obj.dihedralid))#[0][0]
    else: dstart = snap.dihedrals.N

    # snap.dihedrals.resize(dstart+len(obj.dihedrals))
    snap.dihedrals.N = (dstart+len(obj.dihedrals))
    for i, dihedral in enumerate(obj.dihedrals,start=dstart):
        snap.dihedrals.group[i] = pstart-1+dihedral
        snap.dihedrals.typeid[i] = obj.dihedralid

def assign_obj_to_snap(snap,obj,cover=True,particle=True,bond=True,dihedral=False):
    #!!assign bond and dihedral before particle
    #particle number needed for re-index bond and dihedral
    if particle: pstart,pend=assign_snap_particles(snap,obj,cover)
    if bond: assign_snap_bonds(snap,obj,pstart,cover)
    if dihedral: assign_snap_dihedrals(snap,obj,pstart,cover)
    return pstart,pend

def assign_shell(snap):
    shell = create_shell()
    shstart,shend = assign_obj_to_snap(snap,shell,cover=False,dihedral=True)
    return shell

def dihedral_harmonic(theta, kappa, theta0):
    V = 0.5 * kappa * (1 - np.cos( theta - theta0 ))
    F = -0.5 * kappa * np.sin ( theta - theta0 )
    return (V, F)

number_of_bid = 0
hoomd.device.GPU()
sim = hoomd.Simulation(device=hoomd.device.GPU(),seed=1)
snap_1 = initial2()

number_of_bid = snap_1.particles.N
object_of_capsid = assign_shell(snap_1)
BOX = hoomd.Box( Lx=50, Ly=50, Lz=50)
snap_1.configuration.box = BOX
print('object_of_capsid.n_wall' , object_of_capsid.n_wall)
print("Total number of particles : ", snap_1.bonds.N)
sim.create_state_from_snapshot(snap_1)
new_snap = sim.state.get_snapshot()
    
index_p = []
index_c = []
index_w = []
index_s = []
for i in range(len(new_snap.particles.position)):
    if(new_snap.particles.typeid[i] == 0 ):
        index_p.append(i)
    if(new_snap.particles.typeid[i] == 1 ):
        index_c.append(i)
    if(new_snap.particles.typeid[i] == 2 ):
        index_w.append(i)
    if(new_snap.particles.typeid[i] == 3 ):
        index_s.append(i)
        
    
# harmonic = hoomd.md.bond.harmonic()
# harmonic.bond_coeff.set('shell', k = 100*np.sqrt(10) , r0 = 2.1 )
# harmonic.bond_coeff.set('G', k = 50, r0 = 2*param['r_gen'])
# harmonic.bond_coeff.set('spik', k = 100.0, r0 = 2*param['r_ver'])
harmonic = hoomd.md.bond.Harmonic()
harmonic.params['shell'] = dict( k = 100*np.sqrt(10) , r0 = 2.1 )
harmonic.params['G'] = dict( k = 10000, r0 = 2 * param['r_gen'])
harmonic.params['spik'] = dict( k = 100.0, r0 = 2 * param['r_ver'])

scale = 1
# scale_2 = ( 380 / object_of_capsid.n_wall ) * 0.2
nl=hoomd.md.nlist.Cell()
lj = hoomd.md.pair.LJ( nlist = nl )
lj.params[('P', 'P')] = { 'sigma': 2*param['r_gen'] * 1**(-1/6) , 'epsilon' : 0.1*0.1 }
lj.r_cut[('P', 'P')] = 2*param['r_gen']
# excluded volume
lj.params[('P', 'W')] = {'sigma': ((param['r_gen']+param['r_wall'])/1)*2**(-1/6) , 'epsilon' : scale * 0.1}
lj.r_cut[('P', 'W')] = (param['r_gen']+param['r_wall'])

lj.params[('W', 'S')] = {'sigma': (param['r_wall']/1+param['r_terminal'])*2**(-1/6) , 'epsilon' : scale*0.1}
lj.r_cut[('W', 'S')] = (param['r_wall']+param['r_terminal'])

lj.params[('X', 'S')] = {'sigma': (param['r_ver']+param['r_terminal'])*2**(-1/6) , 'epsilon' : scale*0.1 }
lj.r_cut[('X', 'S')] = param['r_ver']+param['r_terminal']

lj.params[('P', 'X')] = {'sigma': (param['r_gen']+param['r_ver'])*2**(-1/6) , 'epsilon' : scale*0.1}
lj.r_cut[('P', 'X')] = (param['r_gen']+param['r_ver'])


# lj.pair_coeff.set('P', 'X', r_cut = (param['r_gen']+param['r_ver']) ,epsilon=scale*0.1, sigma=(param['r_gen']+param['r_ver'])*2**(-1/6),alpha =1)
# lj.pair_coeff.set('S', 'S', r_cut = 2*param['r_ver'], epsilon=0.1, sigma=param['r_ver']*2**(-1/6),alpha =0)
# lj.pair_coeff.set('P', 'P', epsilon=0.0, sigma=0.,alpha =0)
# lj.pair_coeff.set('X', 'X', epsilon=0.2, sigma=2**(1/6),alpha =1)
# non_interaction

lj.params[('S', 'S')] = {'sigma': 0  , 'epsilon' : 0}
lj.r_cut[('S', 'S')] = 0.01
lj.params[('X', 'W')] = {'sigma': 0  , 'epsilon' : 0}
lj.r_cut[('X', 'W')] = 0.01
lj.params[('W', 'W')] = {'sigma': 0  , 'epsilon' : 0}
lj.r_cut[('W', 'W')] = 0.01
lj.params[('X', 'X')] = {'sigma': 0  , 'epsilon' : 0}
lj.r_cut[('X', 'X')] = 0.01

# interaction
# lj.pair_coeff.set('P', 'S', epsilon=10, sigma=()*(2**(-1/6)),alpha =0)
lj.params[('P', 'S')] = {'sigma': ((param['r_gen']+param['r_terminal']))*(2**(-1/6))  , 'epsilon' : scale_2*0.1}
lj.r_cut[('P', 'S')] = 5



# dtable = hoomd.md.dihedral.table(width=10)
# dtable.dihedral_coeff.set('capsomer', func=dihedral_harmonic, coeff=dict(kappa=1/np.sqrt(10), theta0=0.644))

all = hoomd.filter.All()
groupA = hoomd.filter.Type('X')
groupB = hoomd.filter.Type('P')
groupC = hoomd.filter.Type('W')
groupD = hoomd.filter.Type('S')
# groupE = hoomd.group.type('X','W','S')

new_snap = sim.state.get_snapshot()
density_fun(new_snap)

gsd_file = hoomd.write.GSD( trigger = hoomd.trigger.Periodic(2000) , filename = output_file_genome , filter = hoomd.filter.All() , mode = 'wb' )
sim.operations.writers.append( gsd_file )

langevin = hoomd.md.methods.Langevin ( filter = hoomd.filter.Type('P'), kT = 0.1, alpha = 0.2)
integrator = hoomd.md.Integrator(dt=0.001, forces=[harmonic,lj] , methods=[langevin])
sim.operations.integrator = integrator

sim.run(5000000)
new_snap = sim.state.get_snapshot()
count_inside(new_snap)

# integrator = hoomd.md.Integrator(dt=0.01, forces=[h1,lj] , methods=[langevin])
# sim.operations.integrator = integrator
# sim.run(1000)
# new_snap = sim.state.get_snapshot()
# count_inside(new_snap)
t1 = time.time() - t0
print("Time elapsed: ", t1) 
print(sim.tps)
