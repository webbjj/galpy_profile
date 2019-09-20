import matplotlib
matplotlib.use('Agg')

from amuse.lab import *
from amuse.couple import bridge
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.datamodel import Particles

from galpy import potential
from galpy.util import bovy_conversion
import numpy as np

import matplotlib.pyplot as plt

class galpy_profile(LiteratureReferencesMixIn):
    """
    User-defined potential from GALPY
    
    .. [#] Bovy, J; ApJSS, Volume 216, Issue 2, article id. 29, 27 pp. (2015)
    
    """
    def __init__(self,pot, t = 0. | units.Gyr, tgalpy = 0. | units.Gyr, ro=8, vo=220., reverse=False):
        LiteratureReferencesMixIn.__init__(self)

        #Galpy scaling parameters
        self.ro=ro
        self.vo=vo

        #Galpy potential
        self.pot = pot

        #Initialize model time
        self.model_time=t

        #Initialize galpy time
        #Depending on how the galpy potential is setup, one may want to define a different initial galpy time
        #from which to evolve the potential (e.g. negative times) 
        self.tgalpy=tgalpy.value_in(units.Gyr)/bovy_conversion.time_in_Gyr(ro=self.ro,vo=self.vo)
        #Move galpy potential backwards in time?
        self.reverse=reverse

    def evolve_model(self,t_end):

        print('EVOLVE: ',self.model_time.value_in(units.Gyr),self.tgalpy,t_end.value_in(units.Gyr))
        dt=t_end-self.model_time
        self .model_time=t_end  

        if self.reverse:
            self.tgalpy-=dt.value_in(units.Gyr)/bovy_conversion.time_in_Gyr(ro=self.ro,vo=self.vo)
        else:
            self.tgalpy+=dt.value_in(units.Gyr)/bovy_conversion.time_in_Gyr(ro=self.ro,vo=self.vo)

    def get_potential_at_point(self,eps,x,y,z):

        R=np.sqrt(x.value_in(units.kpc)**2.+y.value_in(units.kpc)**2.)
        zed=z.value_in(units.kpc)
        phi=np.arctan2(y.value_in(units.kpc),x.value_in(units.kpc))

        pot=potential.evaluatePotentials(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy,ro=self.ro,vo=self.vo) | units.km**2*units.s**-2
                    
        return pot

    def get_gravity_at_point(self,eps,x,y,z):
    
        R=np.sqrt(x.value_in(units.kpc)**2.+y.value_in(units.kpc)**2.)
        zed=z.value_in(units.kpc)
        phi=np.arctan2(y.value_in(units.kpc),x.value_in(units.kpc))

        Rforce=potential.evaluateRforces(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy)
        phiforce=potential.evaluatephiforces(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy)/(R/self.ro)
        zforce=potential.evaluatezforces(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy)

        ax=(Rforce*np.cos(phi)-phiforce*np.sin(phi))*bovy_conversion.force_in_kmsMyr(ro=self.ro,vo=self.vo) | units.kms * units.myr**-1
        ay=(Rforce*np.sin(phi)+phiforce*np.cos(phi))*bovy_conversion.force_in_kmsMyr(ro=self.ro,vo=self.vo) | units.kms * units.myr**-1
        az=zforce*bovy_conversion.force_in_kmsMyr(ro=self.ro,vo=self.vo) | units.kms * units.myr**-1
        
        return ax,ay,az

    def mass_density(self,x,y,z):
        R=np.sqrt(x.value_in(units.kpc)**2.+y.value_in(units.kpc)**2.)
        zed=z.value_in(units.kpc)
        phi=np.arctan2(y.value_in(units.kpc),x.value_in(units.kpc))

        dens=potential.evaluateDensities(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy,ro=self.ro,vo=self.vo) | units.MSun/(units.parsec**3.)
                        
        return dens

    def circular_velocity(self,r):
        vcirc=potential.vcirc(self.pot,r.value_in(units.kpc)/self.ro,phi=0,ro=self.ro,vo=self.vo) | units.kms
        return vcirc
        
    def enclosed_mass(self,r):
        grav=4.302e-6 #kpc (km/s)^2/Msun
        vc2=potential.vcirc(self.pot,r.value_in(units.kpc)/self.ro,phi=0,t=self.tgalpy,ro=self.ro,vo=self.vo)**2.
        menc= vc2*r.value_in(units.kpc)/grav | units.MSun
        return menc

    def stop(self):
        pass

def setup_cluster(N,Mcluster,Rcluster,Rinit,Vinit):

    #Create star cluster with origin at 0,0,0 and no net velocity
    converter=nbody_system.nbody_to_si(Mcluster,Rcluster)
    stars=new_plummer_sphere(N,converter)

    #Place star cluster in Galactocentric position
    stars.x += Rinit[0]
    stars.y += Rinit[1]
    stars.z += Rinit[2]

    stars.vx += Vinit[0]
    stars.vy += Vinit[1]
    stars.vz += Vinit[2]

    return stars,converter

def evolve_cluster_in_galaxy(N,Mcluster,Rcluster,Rinit,Vinit, galaxy_code, dt, dtout, tend, epsilon):

    #Setup cluster
    stars,converter=setup_cluster(N,Mcluster,Rcluster,Rinit,Vinit)
    stars.scale_to_standard(convert_nbody=converter, smoothing_length_squared = epsilon**2)

    cluster_code=BHTree(converter,number_of_workers=1)     #Change number of workers depending on CPUS available
    cluster_code.parameters.epsilon_squared = epsilon**2
    cluster_code.parameters.opening_angle=0.6
    cluster_code.parameters.timestep=dt
    cluster_code.particles.add_particles(stars)

    #Setup channels between stars particle dataset and the cluster code
    channel_from_stars_to_cluster_code=stars.new_channel_to(cluster_code.particles, attributes=["mass", "x", "y", "z", "vx", "vy", "vz"])    
    channel_from_cluster_code_to_stars=cluster_code.particles.new_channel_to(stars, attributes=["mass", "x", "y", "z", "vx", "vy", "vz"])

    #Setup gravity bridge
    gravity=bridge.Bridge(use_threading=False)
    #stars in cluster_code depend on gravity from galaxy_code
    gravity.add_system(cluster_code, (galaxy_code,))
    #galaxy_code still needs to be added to system so it evolves with time
    gravity.add_system(galaxy_code,)
    #Set how often to update external potential
    gravity.timestep = cluster_code.parameters.timestep/2.

    time=0.0 | tend.unit
    while time<tend:
        print(time.value_in(units.Myr),tend.value_in(units.Myr))
        gravity.evolve_model(time+dt)

        #You need to copy stars from cluster_code to output or analyse:

        #channel_from_cluster_code_to_stars.copy()
        #Output/Analyse
            
        #If you edited the stars particle set, lets say to remove stars from the array because they have 
        #been kicked far from the cluster, you need to copy the array back to cluster_code:
            
        #channel_from_stars_to_cluster_code.copy()
        
        time = gravity.model_time

    #Copy back to stars for final dataset
    channel_from_cluster_code_to_stars.copy()
    gravity.stop()

    return stars

if __name__ == "__main__":
    #Set initial cluster parameters
    N=1000
    Mcluster=1000. | units.MSun
    Rcluster = 10. | units.parsec
    Rinit=[10.,0.,0.] | units.kpc
    Vinit=[0.,220.,0.] | units.km/units.s
    epsilon=0,.75 | units.parsec

    #Setup star cluster simulation
    #Simulation end time
    tend = 1000.0 | units.Myr
    #Frequency of data output
    dtout=5.0 | units.Myr
    #Frequency of star cluster gravity calculation
    dt = 1.0 | units.Myr

    #Set Galactic Potential 
    #- note that initial galpy time can be set to a different value than model_time.
    #- also can tell galpy potential to evolve backwards with reverse flag
    pot=potential.MWPotential2014
    galaxy_code = galpy_profile(pot, tgalpy = 0. | units.Gyr, reverse=False)

    #Evolve star cluster
    stars = evolve_cluster_in_galaxy(N,Mcluster,Rcluster,Rinit,Vinit,galaxy_code, dt, dtout, tend, epsilon)

    #Plot final snapshot
    plt.plot(stars.x.value_in(units.kpc),stars.y.value_in(units.kpc),'.')
    plt.xlabel('X (kpc)')
    plt.ylabel('Y (kpc)')
    plt.savefig('test_cluster.png')
