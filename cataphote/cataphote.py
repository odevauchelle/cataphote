
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Olivier Devauchelle 2022

import numpy as np
from matplotlib.pyplot import gca
from sympy import Ray, Segment, intersection

# from shapely.geometry import LineString

class cataphote :

    def __init__( self, reflector ) :
        self.set_reflector( reflector )

    def set_reflector( self, reflector ) :

        self.reflector = np.array( reflector )
        self.segments = [ Segment( self.reflector[i], self.reflector[i+1] ) for i in range( len(self.reflector)- 1 )  ]

        dx = np.diff( self.reflector, axis = 0 )

        self.tangents = np.array( [ dx/np.sqrt( np.sum( dx**2 ) ) for dx in dx ] )
        self.normals = np.array( [ np.array(  [ [ 0, 1], [ -1, 0 ] ] )@t for t in self.tangents ] )

    def plot_reflector( self, ax = None, **kwargs ):

        if ax is None :
            ax = gca()

        return ax.plot( *np.array( self.reflector ).T, **kwargs )

    def get_intersections( self, ray ) :

        intersections = {}

        for i, segment in enumerate( self.segments ) :
            try :
                if not segment.contains( ray.p1 ) :
                    intersections[i] = intersection( ray, segment )[0]
            except :
                pass

        return intersections

    def get_reflection( self, ray ) :

        distances = []
        segment_indices = []
        intersections = self.get_intersections( ray )

        if len(intersections) > 0 :

            for i, point in intersections.items() :
                distances += [ point.distance( ray.p1 ) ]
                segment_indices += [i]

            segment_index = segment_indices[ np.argmin( distances ) ]

            u = ray.p2 - ray.p1
            n = self.normals[ segment_index ]
            t = self.tangents[ segment_index ]
            u = - np.dot( u, n )*n + np.dot( u, t )*t
            x = intersections[ segment_index ]

            return Ray( x, x + u )

        else :
            return None

    def get_trajectory( self, ray, max_iterations = 30 ) :

        traj = []

        while len(traj) <= max_iterations :

            traj += [ ray.p1 ]

            new_ray = self.get_reflection( ray )

            if new_ray is None :
                break
            else :
                ray = new_ray

        return traj, ray

####################
#
# Try it out
#
####################

if __name__ == '__main__' :

    from pylab import *

    c = cataphote( [ (-.5,0), (-.9,-.7), (0,-.4), (1.2,-.8), (.5,0), (-.5,0) ] )


    ax = gca()
    c.plot_reflector( ax = ax, color = 'k', lw = 2 )

    for angle in -0.2*array( [1]  )*np.pi/2 :

        x = np.array([0,0])

        ray = Ray( x, angle = angle )

        traj, ray = c.get_trajectory( ray, max_iterations = 150 )

        # plot( *np.array( traj + [ ray.p2 ] ).T , label = str(len(traj)))
        plot( *np.array( traj ).T, label = str( len(traj) ), alpha = 1, lw = .3 )

    ax.axis('equal')
    ax.axis('off')
    ax.legend()
    ax.set_xticks([])
    ax.set_yticks([])

    savefig('billiard.svg', bbox_inches = 'tight')

    show()
