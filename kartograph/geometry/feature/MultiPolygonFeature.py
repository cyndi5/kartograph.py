from Feature import Feature
from kartograph.errors import KartographError
from kartograph.simplify.unify import unify_rings
from shapely.geometry import LinearRing, Polygon, MultiPolygon


class MultiPolygonFeature(Feature):
    """ wrapper around shapely.geometry.MultiPolygon """

    def __repr__(self):
        return 'MultiPolygonFeature(), geometry={0}'.format(self.geometry)

    def project_geometry(self, proj):
        """ project the geometry """
        #print 'Projecting geometry with {0}'.format(proj)
        self.geometry = proj.plot(self.geometry)

    def compute_topology(self, point_store, precision=None):
        """
        converts the MultiPolygon geometry into a set of linear rings
        and 'unifies' the points in the rings in order to compute topology
        """
        rings = []
        num_holes = []
        for polygon in self._geoms:
            num_holes.append(len(polygon.interiors))  # store number of holes per polygon
            ext = polygon.exterior.coords
            rings.append(ext)
            for hole in polygon.interiors:
                rings.append(hole.coords)
        self._topology_rings = unify_rings(rings, point_store, precision=precision, feature=self)
        self._topology_num_holes = num_holes


    # Returns the offset from the first canonical point' 
    def get_scale_offset(self, scale_factor=1.):
        """
        returns the offset from some point
        """
      
        lines=[self.geometry.exterior]
        for hole in self.geometry.interiors:
            lines.append(hole)
        offset={'x':0,'y':0}
        temp_offset={}
        new_x=0
        new_y=0
        print 'in offset, len(lines)={0}'.format(len(lines))
        if len(lines)>0 and len(lines[0].coords)>0:
            print '@@doing len stuff, len(lines[0].coords)={0}'.format(len(lines[0].coords))
            for coord in lines[0].coords:
                new_x=(coord[0] * 1.*scale_factor)
                new_y=(coord[1] * 1.*scale_factor)
                temp_offset['x']=-new_x+coord[0]
                temp_offset['y']=-new_y+coord[1]
                if temp_offset['x']*temp_offset['x']+temp_offset['y']*temp_offset['y'] >= offset['x']*offset['x']+offset['y']*offset['y']:
                   offset['x']=temp_offset['x']
                   offset['y']=temp_offset['y']
                print 'offset={0}'.format(offset)
            #for line in lines:
               # print 'Entering new line in offset thing'
                #for coord in line.coords:
                 #   temp_offset['x']=coord[0]*(1-scale_factor)
                  #  temp_offset['y']=coord[1]*(1-scale_factor)
                   # print 'doing scale offset comp: resulting offset={0}'.format(temp_offset)
        else:
            print '@@bad len stuff'
        return offset

    def scale_feature(self, scale_factor=1.,offset={'x':0.,'y':0.}):
        """
        scales the feature and adjusts points by offset 
        """

        lines=[self.geometry.exterior]
        for hole in self.geometry.interiors:
            lines.append(hole)

        new_lines=[]
        curr_line_coords=[]
       # offset['x']=offset['y']=0
        print 'scaling at {1}, offset={0}'.format(offset, scale_factor)
       # print 'self.geometry={0}\n\n'.format(self.geometry)
        for line in lines:
            curr_line_coords=[]
            for coord in line.coords:
                # append the scaled and offsetted coordinate
                curr_line_coords.append((coord[0]*scale_factor+offset['x'],coord[1]*scale_factor+offset['y']))
                print 'coord[0]*scale_factor+offset[\'x\']-coord[0]={0}'.format(coord[0]*scale_factor+offset['x']-coord[0])
            #print '\tcurr_line_coords={0}'.format(curr_line_coords)
            #new_lines.append(line)
            new_lines.append(curr_line_coords)

        self.geometry=Polygon(new_lines[0],new_lines[1:])
    #    print 'self.geometry={0}'.format(self.geometry)
            
    def offset_feature(self, offset={'x':0,'y':0}):
        """
        adjusts points by offset 
        """
        lines=[self.geometry.exterior]
        for hole in self.geometry.interiors:
            lines.append(hole)
        new_lines=[]
        curr_line_coords=[]
        for line in lines:
            curr_line_coords=[]
            for coord in line.coords:
                # append the scaled and offsetted coordinate
              #  print 'coord={0}'.format(coord)
                curr_line_coords.append((coord[0]+offset['x'],coord[1]+offset['y']))
            #print '\tcurr_line_coords={0}'.format(curr_line_coords)
            new_lines.append(curr_line_coords)


        self.geometry=Polygon(new_lines[0],new_lines[1:])
     #   print 'self.geometry={0}'.format(self.geometry)
        
    def break_into_lines(self):
        """
        temporarily stores the geometry in a custom representation in order
        to preserve topology during simplification
        """
        # print '\n\n', self.props['NAME_1'],
        lines = []
        lines_per_ring = []
        for ring in self._topology_rings:
            l = 0  # store number of lines per ring
            s = 0
            i = s + 1
            K = len(ring)
            # print '\n\tnew ring (' + str(K) + ')',
            # find first break-point
            while i < K and ring[i].features == ring[i - 1].features:
                i += 1
            if i == len(ring):  # no break-point found at all
                line = ring  # so the entire ring is treated as one line
                # print len(line),
                lines.append(line)  # store it
                lines_per_ring.append(1)  # and remember that this ring has only 1 line
                continue  # proceed to next ring
            s = i  # store index of first break-point
            a = None
            # loop-entry conditions:
            # - 'a' holds the index of the break-point, equals to s
            # - 's' is the index of the first break-point in the entire ring
            while a != s:
                if a is None:
                    a = s  # if no break-point is set, start with the first one
                line = [a]  # add starting brak-point to line
                i = a + 1  # proceed to next point
                if i == K:
                    i = 0  # wrap around to first point if needed
                while i != s and ring[i].features == ring[((i - 1) + len(ring)) % len(ring)].features:  # look for end of this line
                    line.append(i)  # add point to line
                    i += 1  # proceed to next point
                    if i == K:
                        i = 0  # eventually wrap around
                #if i != s:
                #    line.append(i)  # add end point to line
                l += 1  # increase line-per-ring counter
                a = i  # set end point as next starting point
                #if a == s:  # if next starting point is the first break point..
                # line.append(s)  # append
                # print len(line),
                for ll in range(len(line)):
                    line[ll] = ring[line[ll]]  # replace point indices with actual points
                lines.append(line)  # store line
            lines_per_ring.append(l)

        self._topology_lines_per_ring = lines_per_ring
        return lines

    def restore_geometry(self, lines, minArea=0):
        """
        restores geometry from linear rings
        """
        from shapely.geometry import Polygon, MultiPolygon
        # at first we restore the rings
        rings = []
        isIslands = []
        p = 0
        for l in self._topology_lines_per_ring:
            ring = []
            island = True
            for k in range(p, p + l):
                line = []
                for pt in lines[k]:
                    if len(pt.features) > 1:
                        island = False
                    if not pt.deleted:
                        line.append((pt[0], pt[1]))
                ring += line
            p += l
            rings.append(ring)
            isIslands.append(island)
        # then we restore polygons from rings
        ring_iter = iter(rings)
        islands_iter = iter(isIslands)
        polygons = []
        holes_total = 0
        for num_hole in self._topology_num_holes:
            ext = ring_iter.next()
            island = islands_iter.next()
            holes = []
            while num_hole > 0:
                hole = ring_iter.next()
                islands_iter.next()  # skip island flag for holes
                if len(hole) > 3:
                    holes.append(hole)
                holes_total += 1
                num_hole -= 1
            if len(ext) > 3:
                poly = Polygon(ext, holes)
                if minArea == 0 or not island or poly.area > minArea:
                    polygons.append(poly)

        if len(polygons) > 0:
            self.geometry = MultiPolygon(polygons)
        else:
            self.geometry = None

    def is_simplifyable(self):
        return True

    @property
    def _geoms(self):
        """ returns a list of geoms """
        return hasattr(self.geometry, 'geoms') and self.geometry.geoms or [self.geometry]
