import matplotlib
matplotlib.use('Agg')

from amuse.lab import *
from amuse.couple import bridge
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.datamodel import Particles
from amuse.community.sse.interface import SSE
from amuse.ext.galactic_potentials import *

from galpy import potential
import numpy as np

import matplotlib.pyplot as plt

def evolve_cluster_in_galaxy(N,Rcluster,Rinit,Vinit, galaxy_code, dt, dtout, tend, epsilon):

    #Setup cluster
    masses = new_kroupa_mass_distribution(N)
    Mcluster=masses.sum()

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

    stars.mass=masses

    stars.scale_to_standard(convert_nbody=converter, smoothing_length_squared = epsilon**2)

    #Gravity code
    cluster_code=BHTree(converter,number_of_workers=1)     #Change number of workers depending on CPUS available
    cluster_code.parameters.epsilon_squared = epsilon**2
    cluster_code.parameters.opening_angle=0.6
    cluster_code.parameters.timestep=dt
    cluster_code.particles.add_particles(stars)

    #DEBUG: Apparently all is well if the below line is included:
    Etot_init = cluster_code.kinetic_energy + cluster_code.potential_energy

    #Setup channels between stars particle dataset and the cluster code
    channel_from_stars_to_cluster_code=stars.new_channel_to(cluster_code.particles, attributes=["mass", "x", "y", "z", "vx", "vy", "vz"])    
    channel_from_cluster_code_to_stars=cluster_code.particles.new_channel_to(stars, attributes=["mass", "x", "y", "z", "vx", "vy", "vz"])

    #Stellar evolution code
    stellar_evolution=SSE()
    stellar_evolution.particles.add_particle(stars)

    #Setup channels between stars particle dataset and stellar evolution code
    channel_from_stars_to_stellar_evolution=stars.new_channel_to(stellar_evolution.particles, attributes=["mass"])    
    channel_from_stellar_evolution_to_stars=stellar_evolution.particles.new_channel_to(stars, attributes=["mass"])

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
        print(time.value_in(units.Myr),tend.value_in(units.Myr),stars.mass.sum().value_in(units.MSun))

        stellar_evolution.evolve_model(time+dt/2.)
        channel_from_stellar_evolution_to_stars.copy_attributes(["mass"])
        channel_from_stars_to_cluster_code.copy_attributes(["mass"])

        gravity.evolve_model(time+dt)
        channel_from_cluster_code_to_stars.copy()

        stellar_evolution.evolve_model(time+dt)
        channel_from_stellar_evolution_to_stars.copy_attributes(["mass"])
        channel_from_stars_to_cluster_code.copy_attributes(["mass"])
       
        time = gravity.model_time

    #Copy back to stars for final dataset
    channel_from_cluster_code_to_stars.copy()
    gravity.stop()

    return stars

if __name__ == "__main__":
    #Set initial cluster parameters
    N=1000
    Rcluster = 10. | units.parsec
    Rinit=[10.,0.,0.] | units.kpc
    Vinit=[0.,220.,0.] | units.km/units.s
    epsilon=0.75 | units.parsec

    #Setup star cluster simulation
    #Simulation end time
    tend = 1000.0 | units.Myr
    #Frequency of data output
    dtout=5.0 | units.Myr
    #Frequency of star cluster gravity calculation
    dt = 1.0 | units.Myr

    #Set Galactic Potential 
    galaxy_code = MWpotentialBovy2015()

    #Evolve star cluster
    stars = evolve_cluster_in_galaxy(N,Rcluster,Rinit,Vinit,galaxy_code, dt, dtout, tend, epsilon)

    #Plot final snapshot
    plt.plot(stars.x.value_in(units.kpc),stars.y.value_in(units.kpc),'.')
    plt.xlabel('X (kpc)')
    plt.ylabel('Y (kpc)')
    plt.savefig('test_cluster.png')
