from point import Point


# PolyPoint
# ---
#
# polypo

class PolyPoint(Point):
    def __init__(self,x,y,prev_x,prev_y,next_x,next_y):
        super(PolyPoint,self).__init__(x,y)
        self.prev_x = prev_x
        self.prev_y = prev_y
        self.next_x = next_x
        self.next_y = next_y

    def __getitem__(self, k):
        if k<2:
            return super(PolyPoint,self).__getitem__(k)
        else:
            k_vals = (self.x,self.y,self.prev_x,self.prev_y,self.next_x,self.next_y)
            return k_vals[k]

    def __setitem__(self, k, value):
        if k < 2:
            super(PolyPoint,self).__setitem__(k, value)
        elif k==2:
            self.prev_x = value
        elif k==3:
            self.prev_y = value
        elif k==4:
            self.next_x = value
        elif k==5:
            self.next_y = value
    
