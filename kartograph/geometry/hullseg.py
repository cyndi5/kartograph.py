import sys
from shapely.geometry import LinearRing, Polygon, MultiPolygon, Point
verbose = False

# hullseg
# convenience class for line segment of the convex hull to use for
#

# For a vertical segment, we define left to be below,

# above is whether or not the SEGMENT is above the hull

class hullseg(object):

    #pointA is the first point, pointB the second point, otherPoint is used to compute which way is out
    def __init__(self, point1, point2, otherPoint):
        if point2.x<point1.x:
            self.pointA = point2
            self.pointB = point1
        else:
            self.pointA = point1
            self.pointB = point2
        self.otherPoint = otherPoint

        pA = self.pointA
        pB = self.pointB
        pO = self.otherPoint
        if self.pointA.y = self.pointB.y:
            self.slope = NaN
            self.intercept = self.pointA.y
            self.above = self.pointA.y > self.otherPoint.y
        else:
            self.slope = (pB.y-pA.y)/(pB.x-pA.x)   # Not the slope in the way we expect
            self.intercept = pA.y-self.slope*pA.x
            
            #TODO: self.above = 
        
        
        

    def __repr__(self):
        return 'hullseg(' + self.geometry.__class__.__name__ + ')'

    def __str__(self):
        return 'hullseg('+self.pointA+','+self.pointB+')'
