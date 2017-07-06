
from point import Point


class BBox(object):
    """ 2D bounding box """
    def __init__(self, width=None, height=None, left=0, top=0):
        import sys
        if width == None:
            self.xmin = sys.maxint
            self.xmax = sys.maxint * -1
           # self.width = self.xmax - self.xmin
        else:
            self.xmin = self.left = left
            self.xmax = self.right = left + width
            self.width = width
        if height == None:
            self.ymin = sys.maxint
            self.ymax = sys.maxint * -1
            #self.height = self.ymax - self.ymin
        else:
            self.ymin = self.top = top
            self.ymax = self.bottom = height + top
            self.height = height
    def __repr__(self):
        return '(xmin, xmax, ymin, ymax)=({0},{1},{2},{3})'.format(self.xmin, self.xmax, self.ymin,self.ymax)
    def update(self, pt):
        if not isinstance(pt, Point):
            pt = Point(pt[0], pt[1])
        self.xmin = min(self.xmin, pt.x)
        self.ymin = min(self.ymin, pt.y)
        self.xmax = max(self.xmax, pt.x)
        self.ymax = max(self.ymax, pt.y)


        self.left = self.xmin
        self.top = self.ymin
        self.right = self.xmax
        self.bottom = self.ymax
        self.width = self.xmax - self.xmin
        self.height = self.ymax - self.ymin

    def intersects(self, bbox):
        """ returns true if two bounding boxes overlap """
        return bbox.left < self.right and bbox.right > self.left and bbox.top < self.bottom and bbox.bottom > self.top

    def check_point(self, pt):
        """ check if a point is inside the bbox """
        return pt[0] > self.xmin and pt[0] < self.xmax and pt[1] > self.ymin and pt[1] < self.ymax

    def get_center(self):
        return ((self.xmax+self.xmin)/2, (self.ymax+self.ymin)/2)

    def __str__(self):
        return 'BBox(x=%.6f, y=%.6f, w=%.6f, h=%.6f)' % (self.left, self.top, self.width, self.height)

    def join(self, bbox):
        self.update(Point(bbox.left, bbox.top))
        self.update(Point(bbox.right, bbox.bottom))

    def inflate(self, inflate=1, pad_dict={}):
        d = min(self.width, self.height)
        left_amount = d * inflate * pad_dict["left"]
        top_amount = d * inflate * pad_dict["top"]
        right_amount = d * inflate * pad_dict["right"]
        bottom_amount = d * inflate * pad_dict["bottom"]
 #       print 'amount={0}'.format(amount)
 #       print self.xmin,self.ymin, self.xmax, self.ymax
        self.xmin -= left_amount
        self.ymin -= top_amount
        self.xmax += right_amount
        self.ymax += bottom_amount

#        print self.xmin,self.ymin, self.xmax, self.ymax


        self.left = self.xmin
        self.top = self.ymin
        self.right = self.xmax
        self.bottom = self.ymax
        self.width = self.xmax - self.xmin
        self.height = self.ymax - self.ymin

    # New more finely-grained inflate, allows inflation of each side

    def scale(self, scale_factor=1, offset={'x':0,'y':0}):
        d = min(self.width, self.height)
 #       print 'amount={0}'.format(amount)
 #       print self.xmin,self.ymin, self.xmax, self.ymax
        self.xmin = self.xmin*scale_factor+offset['x']
        self.ymin = self.ymin*scale_factor+offset['y']
        self.xmax = self.xmax*scale_factor+offset['x']
        self.ymax = self.ymax*scale_factor+offset['y']

#        print self.xmin,self.ymin, self.xmax, self.ymax


        self.left = self.xmin
        self.top = self.ymin
        self.right = self.xmax
        self.bottom = self.ymax
        self.width = self.xmax - self.xmin
        self.height = self.ymax - self.ymin       

    def __getitem__(self, k):
        if k == 0:
            return self.xmin
        if k == 1:
            return self.ymin
        if k == 2:
            return self.xmax
        if k == 3:
            return self.ymax
        return None
