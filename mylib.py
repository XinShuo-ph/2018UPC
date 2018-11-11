import numpy as np

phi0=0; sigp=1.25; sigb=0.0;


#define parameter:
pi = np.pi
sail_mass = 0.01 #kg
c = 3e8 #Light Speed
P = 50000000000 #Laser array power
S = 10 #area of space craft
R = np.sqrt(S/pi)
alpha = np.pi/6 #angle of space craft
h = R*np.tan(alpha) #height of space craft
LaserDist_r = np.linspace(0,1-sigp/R,5)
LaserDist_r = LaserDist_r*R
LaserDis_theta = np.linspace(0,1,6)
LaserDis_theta = LaserDis_theta*2*pi
N = len(LaserDis_theta) * len(LaserDist_r) #Laser beam number

##define caculating function:
def PowerFlux(r,beta,sigma):
    P0 = P/N
    total_flux = 0
    r_m,theta_m = np.meshgrid(LaserDist_r,LaserDis_theta)
    distance = np.sqrt(r**2+r_m**2-2*r*r_m*np.cos(beta-theta_m))
    flux = P0*np.exp(-distance**2/(2*sigma)**2)  /np.sqrt(2*pi*sigma**2)
    total_flux = np.sum(np.sum(flux))
    return total_flux


def n(alpha,phi,gamma):
    nx = np.cos(alpha)*np.sin(phi)+np.sin(alpha)*np.cos(gamma)*np.cos(phi)
    ny = np.sin(alpha)*np.sin(gamma)
    nz = np.cos(alpha)*np.cos(phi)-np.sin(alpha)*np.cos(gamma)*np.sin(phi)
    vect = np.array([nx,ny,nz])
    return vect

def b(sigma):
    bx = 0 + np.random.normal(0,sigma,1)
    bx = bx[0]
    by = 0 + np.random.normal(0,sigma,1)
    by = by[0]
    bz = 1 + np.random.normal(0,sigma,1)
    bz = bz[0]
    len_b = np.sqrt(bx**2 + by**2 + bz**2)
    vec_b = np.array([bx,by,bz])/len_b
    return vec_b

def Force(x,y,z,K,phi,sigma_P,sigma_b):

    int_rho = np.linspace(0,R,100)
    delta_rho = int_rho[1] - int_rho[0]
    int_gamma = np.linspace(0,2*pi,100)
    delta_gamma = int_gamma[1] - int_gamma[0]
    total_force = np.array([0,0,0])
    Gamma,Rho = np.meshgrid(int_gamma,int_rho)
    x_bias = x + Rho*np.sin(Gamma)
    y_bias = y + Rho*np.cos(alpha)*np.cos(phi)- (h/2-Rho*np.tan(alpha))*np.sin(phi)
    r_int = np.sqrt(x_bias**2 + y_bias**2)
    theta_int = pi/2*np.ones((len(x_bias),len(x_bias)))
    P_xyz = np.zeros((len(r_int),len(theta_int)))
    for mm in range(len(r_int)):
        for nn in range(len(theta_int)):
            if x_bias[mm,nn] != 0:
                theta_int[mm,nn] = np.arctan(y_bias[mm,nn]/x_bias[mm,nn])
            P_xyz[mm,nn] = PowerFlux(r_int[mm,nn],theta_int[mm,nn],sigma_P)
    n_xyz = n(alpha,phi,Gamma)
    b_xyz = b(sigma_b)
    len_K = np.sqrt(np.dot(K,K))
    alpha1 = np.tensordot(K,n_xyz,axes=[[0],[0]])
    alpha1 = np.arccos(alpha1/(np.dot(K,K)))
    beta1 = np.tensordot(b_xyz,n_xyz,axes=[[0],[0]])
    beta1 = np.arccos(beta1/(np.dot(b_xyz,b_xyz)))
    Pf = (np.sqrt(sail_mass**2*c**4+len_K**2*c**2)-len_K*np.cos(alpha1-beta1))/((1-np.cos(2*beta1))+np.sqrt(sail_mass**2*c**4+len_K**2*c**2)/P_xyz+len_K*(np.cos(alpha1+beta1)/P_xyz))
    Force_parr = (P_xyz+Pf)/c
    Force_parr = Force_parr*np.cos(2*beta1)*Rho*delta_rho*delta_gamma/np.cos(alpha)
    totalForce_parr = np.sum(np.sum(Force_parr))
    Force_perp = Pf*np.sin(2*beta1)/c
    Force_perp = Force_perp * Rho*delta_rho*delta_gamma/np.cos(alpha)
    totalForce_perp = np.sum(np.sum(Force_perp))
    total_force = np.array([totalForce_parr,totalForce_perp])
    return total_force



def equation(var,t):
    diff=np.array([0.,0.,0.,0.,0.,0.])
    diff[0:3]=var[3:6]
    diff[2]=var[5]/sail_mass/np.sqrt( 1+(var[3]**2+var[4]**2+var[5]**2)/c**2 )
    F=Force(var[0],var[1],var[2],[var[3],var[4],var[5]],phi0,sigp,sigb)
    diff[3]=F[1]
    diff[5]=F[0]
    return diff