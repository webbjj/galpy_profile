# galpy_profile

The purpose of the galpy_profile is to be able to use any potential defined in GALPY as a galaxy code in AMUSE. The current version, while not thoroughly tested, appears to work for static, time dependent, and moving GALPY potentials. An example of how one would initialize the galpy potential in amuse is as follows:

#Define GALPY potential

from galpy.potential import MWPotential2014
pot=MWPotential2014

#Define a galaxy potential in AMUSE
galaxy_code = galpy_profile(pot,t = 0. | units.Gyr, tgalpy = 0. | units.Gyr, ro=8, vo=220.)

Note that galpy_profile has three additional inputs than a standard AMUSE potential:

tgalpy - the time at which your GALPY potential is being defined. This input is important for time dependent potentials that are setup such that they require a negative value of tgalpy, but the AMUSE N-body integrator you are using requires t>0. 

ro,vo - these are the standard GALPY scaling relations for distance and velocity. Default values are 8. and 220.

#For use with an N-body star cluster simulation
#Assuming one has already initialized a code for simulation a star cluster (cluster_code), the star cluster code and galpy potential can be bridged together such that the GALPY potential is taken to be the external potential that cluster stars experience:
from amuse.couple import bridge

gravity=bridge.Bridge()
gravity.add_system(cluster_code, (galaxy_code,))
gravity.add_system(galaxy_code,)

And example can be found in test_cluster.py
