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
    def __init__(self, point1, point2, otherPoint, distParam=None):
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
        pO = self.otherPoint
        if self.pointA.y == self.pointB.y:
            self.slope = NaN
            self.intercept = self.pointA.y
            self.above = self.pointA.y > self.otherPoint.y
            self.midPoint = Point((pA.x+pB.x)/2,(pA.y+pB.y)/2)
            if distParam is not None:
                if self.above:
                    self.outPoint = Point(self.midPoint.x+distParam,self.midPoint.y)
                else:
                    self.outPoint = Point(self.midPoint.x-distParam,self.midPoint.y)
        else:
            self.slope = (pB.y-pA.y)/(pB.x-pA.x)   # Not the slope in the way we expect
            self.intercept = pA.y-self.slope*pA.x
            self.above = pO.y > pO.x*self.slope+self.intercept
            self.midPoint = Point((pA.x+pB.x)/2,(pA.y+pB.y)/2)
            if distParam is not None and self.slope == 0:
            # Deal with 0 separately because 1/0 = not a number 
                if self.above:
                    self.outPoint = Point(self.midPoint.x,self.midPoint.y-distParam)
                else:
                    self.outPoint = Point(self.midPoint.x,self.midPoint.y+distParam)
            elif distParam is not None:
                inv_slope = -1./self.slope
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
                else
                    self.outPoint = Point(self.midPoint.x-distParam*x_factor,
                                              self.midPoint.y+distParam*y_factor)


                                              
        
        

    def __repr__(self):
        return 'hullseg(' + self.geometry.__class__.__name__ + ')'

    def __str__(self):
        return 'hullseg('+self.pointA+','+self.pointB+')'
