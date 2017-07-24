import sys
from shapely.geometry import LinearRing, Polygon, MultiPolygon, Point
from math import sqrt, cos, sin, atan
verbose = False

# hullseg
# convenience class for line segment of the convex hull to use for
#

# For a vertical segment, we define left to be below,

# above is whether or not the SEGMENT is above the hull

class hullseg(object):
    #pointA is the first point, pointB the second point, otherPoint is used to compute which way is out
    # distParam is none for the side hull 
    def __init__(self, point1, point2, otherPoint, distParam=None, point_pos = None):
        self.point_pos = point_pos
        if point2.x<point1.x:
            self.pointA = point2
            self.pointB = point1
        else:
            self.pointA = point1
            self.pointB = point2
        self.otherPoint = otherPoint
        self.distParam = distParam
        self.outPoint = None
        pA = self.pointA
        pB = self.pointB

        self.length = sqrt((pB.x-pA.x)**2+(pB.y-pA.y)**2)
        pO = self.otherPoint
        if self.pointA.x == self.pointB.x:
            self.slope = sgn(pB.y-pA.y)*float('inf')
            self.intercept = self.pointA.x
            self.above = self.pointA.x >= self.otherPoint.x
            self.midPoint = Point((pA.x+pB.x)/2,(pA.y+pB.y)/2)
            if distParam is not None:
                self.inv_slope = 0
                if self.above:
                    self.outPoint = Point(self.midPoint.x+distParam,self.midPoint.y)
                else:
                    self.outPoint = Point(self.midPoint.x-distParam,self.midPoint.y)
                self.inv_intercept = self.outPoint.y
            else:
                print('Error: distparam is None')
        else:
            self.slope = (pB.y-pA.y)/(pB.x-pA.x)   # Not the slope in the way we expect
            self.intercept = pA.y-self.slope*pA.x
            self.above = pO.y > pO.x*self.slope+self.intercept
            self.midPoint = Point((pA.x+pB.x)/2,(pA.y+pB.y)/2)
            if distParam is not None and self.slope == 0:
                self.inv_slope = (1 if self.above else -1)*float('inf')
                # Deal with 0 separately because 1/0 = not a number
                if self.above:
                    self.outPoint = Point(self.midPoint.x,self.midPoint.y-distParam)
                else:
                    self.outPoint = Point(self.midPoint.x,self.midPoint.y+distParam)
                self.inv_intercept = self.outPoint.x
            elif distParam is not None:
                self.inv_slope = inv_slope = -1./self.slope
                x_factor = 1./abs(sqrt(inv_slope**2+1))
                y_factor = abs(inv_slope)/abs(sqrt(inv_slope**2+1))
                if self.above and inv_slope>=0:
                    # up and to the left based on how slope is
                    self.outPoint = Point(self.midPoint.x-distParam*x_factor,
                                              self.midPoint.y-distParam*y_factor)
                elif self.above and inv_slope < 0:
                    self.outPoint = Point(self.midPoint.x+distParam*x_factor,
                                              self.midPoint.y-distParam*y_factor)
                elif not self.above and inv_slope >= 0:
                    self.outPoint = Point(self.midPoint.x+distParam*x_factor,
                                              self.midPoint.y+distParam*y_factor)
                else:
                    self.outPoint = Point(self.midPoint.x-distParam*x_factor,
                                              self.midPoint.y+distParam*y_factor)
                self.inv_intercept = self.outPoint.y - self.inv_slope * self.outPoint.x


                                              
        
        
    
    def __repr__(self):
        return 'hullseg(' + self.geometry.__class__.__name__ + ')'

    def __str__(self):
        return 'hullseg( points=({0},{1}), above={2}, slope={3}, midPoint={4}, outPoint={5},   length={6})'.format(self.pointA, self.pointB, self.above, self.slope, self.midPoint, self.outPoint,self.length)

    # check if point is above (or on) line (left is below) 
    def is_above(self, pt, slope, intercept):
        if abs(slope) == float('inf'):
            return pt.x >= intercept
        else:
            return pt.x*slope+intercept >= pt.y
        
    # check if point is below (or on) line (left is below) 
    def is_below(self, pt, slope, intercept):
        if abs(slope) == float('inf'):
            return pt.x <= intercept
        else:
            return pt.x*slope+intercept <= pt.y
    
            
        
    # Check to see if the point is ``good" relative to this hull segment
    def check_point(self, s_point_pos, m_coords, s_coords):
        s_point = Point(s_coords[s_point_pos][0],s_coords[s_point_pos][1])
        curr_slope=self.slope
        if abs(curr_slope) == float('inf'):
            curr_intercept = s_point.x
        else:
            curr_intercept = s_point.y-curr_slope*s_point.x

        # Check where the points next to it are
        if s_point_pos == 0:
            s_before = len(s_coords)-2
        else:
            s_before=s_point_pos-1
        s_after=(s_point_pos+1+len(s_coords)) % len(s_coords)
        s_point_before = Point(s_coords[s_before][0],s_coords[s_before][1])
        s_point_after = Point(s_coords[s_after][0],s_coords[s_after][1])
        if self.above:
            # Want s_point_before and s_point_after to be above line too
            return self.is_above(s_point_before, curr_slope, curr_intercept) and self.is_above(s_point_after, curr_slope, curr_intercept)
        else:
            return self.is_below(s_point_before, curr_slope, curr_intercept) and self.is_below(s_point_after, curr_slope, curr_intercept)

    def sgn(x):
        if x ==0:
            return 0
        elif x>0:
            return 1
        else:
            return -1
    
            
