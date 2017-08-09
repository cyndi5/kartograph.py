"""
geometry utils
"""

from hullseg import hullseg
from feature import create_feature
from shapely.geometry import Point, Polygon, LineString, MultiPolygon
from math import atan, cos
from copy import deepcopy
def is_clockwise(pts):
    """ returns true if a given linear ring is in clockwise order """
    s = 0
    for i in range(len(pts) - 1):
        if 'x' in pts[i]:
            x1 = pts[i].x
            y1 = pts[i].y
            x2 = pts[i + 1].x
            y2 = pts[i + 1].y
        else:
            x1, y1 = pts[i]
            x2, y2 = pts[i + 1]
        s += (x2 - x1) * (y2 + y1)
    return s >= 0

# returns 1 if positive, 0 if 0, -1 otherwise
def sgn(x):
    if x ==0:
        return 0
    elif x>0:
        return 1
    else:
        return -1

def bbox_to_polygon(bbox):
    from shapely.geometry import Polygon
    s = bbox
    poly = Polygon([(s.left, s.bottom), (s.left, s.top), (s.right, s.top), (s.right, s.bottom)])
    return poly


def geom_to_bbox(geom, min_area=0):
    from kartograph.geometry import BBox
    from shapely.geometry import MultiPolygon
    if min_area == 0 or not isinstance(geom, MultiPolygon):
        # if no minimum area ratio is set or the geometry
        # is not a multipart geometry, we simply use the
        # full bbox
        minx, miny, maxx, maxy = geom.bounds
        return BBox(width=maxx - minx, height=maxy - miny, left=minx, top=miny)
    else:
        # for multipart geometry we use only the bbox of
        # the 'biggest' sub-geometries, depending on min_area
        bbox = BBox()
        areas = []
        bb = []
        for polygon in geom.geoms:
            areas.append(polygon.area)
        max_a = max(areas)
        for i in range(len(geom.geoms)):
            a = areas[i]
            if a < max_a * min_area:
                # ignore this sub polygon since it is too small
                continue
            bb.append(geom.geoms[i].bounds)
    for b in bb:
        bbox.update((b[0], b[2]))
        bbox.update((b[1], b[2]))
        bbox.update((b[0], b[3]))
        bbox.update((b[1], b[3]))
    return bbox


def join_features(features, props, buf=False):
    """ joins polygonal features
    """
    from feature import MultiPolygonFeature, MultiLineFeature
    from shapely.ops import linemerge

    if len(features) == 0:
        return features

    joined = []
    polygons = []
    lines = []

    for feat in features:
        if isinstance(feat, MultiPolygonFeature):
            polygons.append(feat.geom)
        elif isinstance(feat, MultiLineFeature):
            lines.append(feat.geom)
        else:
            joined.append(feat)  # cannot join this

    polygons = filter(lambda x: x is not None, polygons)
    if len(polygons) > 0:
        poly = polygons[0]
        if buf is not False:
            poly = poly.buffer(buf, 4)
        for poly2 in polygons[1:]:
            if buf is not False:
                poly2 = poly2.buffer(buf, 4)
            poly = poly.union(poly2)
        joined.append(MultiPolygonFeature(poly, props))

    if len(lines) > 0:
        rings = []
        for line in lines:
            geoms = hasattr(line, 'geoms') and line.geoms or [line]
            rings += geoms
        joined.append(MultiLineFeature(linemerge(rings), props))
    return joined

# Deal with the offset coords in old way that doesn't look good for some places
def get_offset_coords(mainbbox, sidebbox, position_factor):
        m_bounds = mainbbox
        s_bounds = sidebbox
        left_good = s_bounds.width + m_bounds.width  <= max(m_bounds.height,s_bounds.height)*position_factor
        next_side_offset={'x': 0, 'y': 0}
        if left_good:
            # Add a little breathing room on the left
            print 'Adding on left'
            next_side_offset['x'] = -s_bounds.left+m_bounds.left-m_bounds.width/10.-s_bounds.width
            next_side_offset['y'] = -s_bounds.top+m_bounds.top+m_bounds.height/2-s_bounds.height/2.
        else:
            print 'Adding on top'
            next_side_offset['x'] = -s_bounds.left+m_bounds.left+m_bounds.width/2.-s_bounds.width/2.
            next_side_offset['y'] = -s_bounds.top+m_bounds.top-m_bounds.height/10.-s_bounds.height

        return next_side_offset['x'], next_side_offset['y']

# Deal with offset coords in a more sophisticated way
def get_offset_coords_complex(mainbbox, sidebbox, main_geom, side_geom, position_factor, the_map):
    opts = the_map.options
    data = opts['bounds']['data']
    sidelayer = the_map.layersById[data['sidelayer']]
    m_hull = main_geom#geom_to_bbox(main_geom,min_area=0)
    s_hull = side_geom#geom_to_bbox(side_geom,min_area=0)
    dist_param = min(mainbbox.width, mainbbox.height)/25.

    

    best_x=0
    best_y=0
    min_area=float('inf')
    max_diff=0
    min_max_side=float('inf')
    
    m_coords = m_hull.exterior.coords
    s_coords = s_hull.exterior.coords

    hull_list=[]
    side_hull_list=[]
    main_perim=0
    for i in range(len(m_coords)-1):
        pA=Point(m_coords[i % len(m_coords)][0], m_coords[i % len(m_coords)][1])
        b_pt_val=i+1
        o_pt_val = i+2 if i<len(m_coords)-2 else i+3
        pB=Point(m_coords[(b_pt_val) % len(m_coords)][0], m_coords[(b_pt_val) % len(m_coords)][1])
        pO=Point(m_coords[(o_pt_val) % len(m_coords)][0], m_coords[(o_pt_val) % len(m_coords)][1])
        #print('{0}, {1}, {2}'.format(pA,pB,pO))
        temp_hull_seg=hullseg(pA,pB,pO,m_hull,distParam=dist_param, min_area=data['min-area'])
        # print('current segment slope={0}'.format(temp_hull_seg.slope))
        main_perim+=temp_hull_seg.length
        hull_list.append(temp_hull_seg)

 
    best_my_box = None
    best_seg = None
    for seg in hull_list:
        # Find the rotated square bounding box with slope of some side according to seg
        if seg.length < 0 * main_perim/len(hull_list):
            continue
        my_box = seg.get_box_pts(s_hull)
        sidelayer = the_map.layersById[data['sidelayer']]
        box_coords = my_box.exterior.coords
        # coordinate order is left top right bottom
        if seg.above and seg.slope >= 0:
            temp_x_offset = seg.outPoint.x - (box_coords[3][0]+box_coords[0][0])/2.
            temp_y_offset = seg.outPoint.y - (box_coords[3][1]+box_coords[0][1])/2.
        elif seg.above and seg.slope < 0:
            temp_x_offset = seg.outPoint.x - (box_coords[3][0]+box_coords[0][0])/2.
            temp_y_offset = seg.outPoint.y - (box_coords[3][1]+box_coords[0][1])/2.
        elif not seg.above and seg.slope >= 0:
            temp_x_offset = seg.outPoint.x - (box_coords[2][0]+box_coords[1][0])/2.
            temp_y_offset = seg.outPoint.y - (box_coords[2][1]+box_coords[1][1])/2.
        elif not seg.above and seg.slope < 0:
            temp_x_offset = seg.outPoint.x - (box_coords[1][0]+box_coords[2][0])/2.
            temp_y_offset = seg.outPoint.y - (box_coords[1][1]+box_coords[2][1])/2.
        else:
            continue

        my_box_2=Polygon([(box_coords[0][0]+temp_x_offset,box_coords[0][1]+temp_y_offset),(box_coords[1][0]+temp_x_offset,box_coords[1][1]+temp_y_offset),(box_coords[2][0]+temp_x_offset,box_coords[2][1]+temp_y_offset),(box_coords[3][0]+temp_x_offset,box_coords[3][1]+temp_y_offset),(box_coords[0][0]+temp_x_offset,box_coords[0][1]+temp_y_offset)])
        temp_bbox=geom_to_bbox(my_box_2,data['min-area'])
     #   print('temp_bbox={0}'.format(temp_bbox))
        temp_bbox.join(mainbbox)
      #  print('now temp_bbox={0}'.format(temp_bbox))

        #print('hullseg={0}\n\tarea={1}\n\n'.format(seg,temp_bbox.area()))
        temp_max_side = max(temp_bbox.height,temp_bbox.width)
        if temp_max_side < min_max_side: #temp_bbox.area() < min_area:
            x_offset = temp_x_offset
            y_offset = temp_y_offset
            min_area = temp_bbox.area()
            min_max_side = temp_max_side
            #print('New temp_bbox = {0}, seg.above={1}, seg.slope={2}, min_area={3}, min_max_side={4}'.format(temp_bbox,seg.above,seg.slope, min_area, min_max_side))
            best_my_box = deepcopy(my_box)
            best_pointA = deepcopy(seg.pointA)#.buffer(dist_param/4.)
            best_pointB = deepcopy(seg.pointB)
            best_out_point = deepcopy(seg.outPoint)
            
    # Deal with checking how good it is offset 
    
    return x_offset, y_offset

def get_complex_hull(mainbbox, sidebbox, geom, position_factor, the_map):
    opts = the_map.options
    data = opts['bounds']['data']
    sidelayer = the_map.layersById[data['sidelayer']]
    m_hull = geom.convex_hull#geom_to_bbox(geom,min_area=0)
 #   s_hulls = side_geom#geom_to_bbox(side_geom,min_area=0)
    max_area=0
    best_poly = None
    coord_list=[]
    poly_list=[]

    ret_poly_list=[]
    
    if not isinstance(geom, MultiPolygon):
        poly_list = [geom]
    else:
        poly_list = geom
    max_area = max(poly.area for poly in poly_list)
    big_poly_list = [poly for poly in poly_list if poly.area*50000 >= max_area]
    next_poly_list=[]
    
    for poly in big_poly_list:
        pt_0=poly.exterior.coords[0]
        #coord_list=[]
        for coord in poly.exterior.coords[:-1]:
            coord_list.append(coord)
        
    #curr_poly = Polygon(coord_list).convex_hull
    #ret_poly_list.append(curr_poly)

    curr_poly = Polygon(convex_hull_jacob(coord_list, big_poly_list))

    print 'Bounds(curr_poly)={0}'.format(curr_poly.bounds)
    ret_poly_list.append(curr_poly)
        

    if len(ret_poly_list)>1:
        return MultiPolygon(ret_poly_list)
    else:
        return ret_poly_list[0]


def convex_hull_graham(points):
    '''
    Returns points on convex hull in CCW order according to Graham's scan algorithm. 
    By Tom Switzer <thomas.switzer@gmail.com>.
    '''
    TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

    MIN_DIFF=1e-7
    def near_eq(a,b):
        return abs(a-b) <= MIN_DIFF
    
    def cmp(a, b):
        return (a > b) - (a < b)

    def turn(p, q, r):
        return cmp((q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]), 0)

    def _keep_left(hull, r):
        while len(hull) > 1 and turn(hull[-2], hull[-1], r) != TURN_LEFT:
            hull.pop()
        if not len(hull) or hull[-1] != r:
            hull.append(r)
        return hull

    points = sorted(points)
    l = reduce(_keep_left, points, [])
    u = reduce(_keep_left, reversed(points), [])
    return l.extend(u[i] for i in range(1, len(u) - 1)) or l
    
def convex_hull_jacob(points, big_poly_list):
    '''
    Returns points on convex hull in CCW order according to Graham's scan algorithm. 
    By Jacob Alperin-Sheriff
    '''
    from math import atan2, sqrt, pi
    from bisect import bisect_left, bisect, bisect_right
    TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

    MIN_DIFF=1e-7
    def near_eq(a,b):
        return abs(a-b) <= MIN_DIFF

    def near_lt(a,b):
        return a + MIN_DIFF <= b

    def near_gt(a,b):
        return a - MIN_DIFF >= b

    def near_le(a,b):
        return not near_gt(a,b)

    def near_ge(a,b):
        return not near_lt(a,b)
    
    def index(a, x):
        'Locate the leftmost value exactly equal to x'
        i = bisect_left(a, x)
        if i != len(a) and a[i] == x:
            return i
        raise ValueError

    def find_lt(a, x):
        'Find index of rightmost value less than x'
        i = bisect_left(a, x)
        if i:
            return i-1
        raise ValueError

    def find_le(a, x):
        'Find index of rightmost value less than or equal to x'
        i = bisect_right(a, x)
        if i:
            return i-1
        raise ValueError

    def find_gt(a, x):
        'Find index of leftmost value greater than x'
        i = bisect_right(a, x)
        if i != len(a):
            return i
        raise ValueError

    def find_ge(a, x):
        'Find index of leftmost item greater than or equal to x'
        i = bisect_left(a, x)
        if i != len(a):
            return i
        raise ValueError
        
    def length(a,b):
        return sqrt((a[1]-b[1])**2+(a[0]-b[0])**2)


    def cmp(a, b):
        return (a > b) - (a < b)

    def turn(p, q, r):
        return cmp((q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]), 0)
    
    def _keep_left(hull, r):
        while len(hull) > 1 and turn(hull[-2], hull[-1], r) != TURN_LEFT:
            hull.pop()
        if not len(hull) or hull[-1] != r:
            hull.append(r)
        return hull

    def turn_real(p, q, r):
        return ((q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]))/(length(p,q)*length(q,r))

    def turn_angle(p,q,r):
        return atan2(r[1]-q[1],r[0]-q[0])-atan2(q[1]-p[1],q[0]-p[0])

    def pt_eval(ptA, ptB, to_eval):
        y_int = get_y_intersection(to_eval[0], ptA,ptB)
        x_int = get_x_intersection(to_eval[1],ptA,ptB)
        return abs(y_int - to_eval[1])*abs(x_int-to_eval[0])
        
        
    
    def print_angles(points):
        pt0 = points[0]
        for i in range(0,len(points)-1):
            print('({0: 8.2f}, {1: 8.2f}),({2: 8.2f},{3: 8.2f}) = {4: 8.2f}'.format(points[i][0],points[i][1],points[i+1][0],points[i+1][1], atan2(points[i+1][1]-points[i][1],points[i+1][0]-points[i][0]) )  )
    
    def get_best_pt(keys, to_search, ptA, ptB, is_above, is_left, ang):
        this_coord=1
        other_coord=0
        if is_above:
            this_func = lambda a,b: min(a,b)
            other_func = lambda a,b: max(a,b)
        else:
            this_func = lambda a,b: max(a,b)
            other_func = lambda a,b: min(a,b)
        lo = find_gt(keys,min(ptA[other_coord],ptB[other_coord]))
        hi= find_lt(keys,max(ptA[other_coord],ptB[other_coord]))
        extremum = other_func(ptA[this_coord],ptB[this_coord])
        best_pt = (ptB[0],ptB[1])
        best_eval = pt_eval(ptA, ptB, best_pt)
        # print('ptA={0}, ptB={1}'.format(ptA,ptB))
        # print('extremum={0}'.format(extremum))
        the_iter=range(lo,hi+1,1) if bool(is_left) else range(hi,lo-1,-1)
        good_pts=[]
        cross_bool = bool(not is_left) ^ bool(is_above)
        for i in the_iter:
            if (is_above and to_search[i][this_coord] < extremum) or (not is_above and to_search[i][this_coord] > extremum):
                good_pts.append((to_search[i][0], to_search[i][1]))
                extremum = to_search[i][this_coord]
            elif (not is_above and to_search[i][this_coord] < extremum):
               good_pts=[pt for pt in good_pts if pt[1] > to_search[i][this_coord]]
                    
        # the_iter=range(lo,hi+1,1)#if bool(not is_left) ^ bool(is_above) else range(hi,lo-1,-1)
        extremum = other_func(ptA[this_coord],ptB[this_coord])           
        for curr_pt in good_pts:

            temp_eval = pt_eval(ptA,ptB,curr_pt)

            if temp_eval > best_eval:
                best_pt = (curr_pt[0], curr_pt[1])
                best_eval = temp_eval
        # print('(is_above={0},is_left={1}, ang={4}), best_pt={2}, temp_eval={3}'.format(is_above,is_left,best_pt,temp_eval, ang/pi))

            extremum = curr_pt[this_coord]
        # print('')
        return best_pt

    
    def _line_angle(a,b):
        return atan2(b[1]-a[1],b[0]-a[0])

    def get_x_intersection(y,ptA,ptB):
        if ptA[0]==ptB[0]:
            raise ValueError
        slope = (ptB[1]-ptA[1])/(ptB[0]-ptA[0])
        intercept = ptB[1]-slope*ptB[0]
        return y/slope - intercept/slope

    def get_y_intersection(x,ptA,ptB):
        if ptA[0]==ptB[0]:
            raise ValueError
        slope = (ptB[1]-ptA[1])/(ptB[0]-ptA[0])
        intercept = ptB[1]-slope*ptB[0]
        return x*slope + intercept

    
    def _add_hull(points, hull, indiv_hulls):
        # Need to check if two points on convex hull are on same simple polygon or not??
        ret_hull=[]
        x_points = points
        y_points = sorted(deepcopy(points), key = lambda pt: pt[1])
        
        x_keys = [pt[0] for pt in points]
        y_keys = [pt[1] for pt in y_points]
        for i in range(0,len(hull)):
            j=(i+1) % len(hull)
            ang = _line_angle(hull[i],hull[j])
            ret_hull.append(hull[i])
            temp_pt1 = temp_pt2 = temp_pt3 = None
            print('hull[{0}]={1}, hull[{2}]={3}'.format(i,hull[i],j,hull[j]))
            if near_lt(0,ang) and near_lt(ang,pi/2):
                print('0 < ang < pi/2')
                min_pt=get_best_pt(x_keys, x_points,hull[i],hull[j],True,False,ang)
                print '\tmin_pt={0}'.format(min_pt)

                if near_ge(min_pt[0],max(hull[i][0],hull[j][0])):
                    print '\t near_ge'
                    temp_pt1=(hull[i][0],min_pt[1])
                    temp_pt2 = (get_x_intersection(min_pt[1],hull[i],hull[j]),min_pt[1])
                else:
                    print '\t not near_ge'
                    temp_pt1 = (min_pt[0], get_y_intersection(min_pt[0], hull[i],hull[j]))
                    temp_pt2 = (min_pt[0], min_pt[1])
                    temp_pt3 = (get_x_intersection(min_pt[1],hull[i],hull[j]),min_pt[1])

            elif near_lt(-pi/2,ang) and near_lt(ang,0):
                print('-pi/2 < ang < 0')

                min_pt=get_best_pt(x_keys, x_points,hull[i],hull[j],True,True,ang)
                print '\tmin_pt={0}'.format(min_pt)

                if near_ge(min_pt[0],hull[j][0]):
                    temp_pt1 = (get_x_intersection(min_pt[1],hull[i],hull[j]),min_pt[1])
                    temp_pt2=(hull[j][0],min_pt[1])
                else:
                    temp_pt1 = (get_x_intersection(min_pt[1],hull[i],hull[j]),min_pt[1])
                    temp_pt2 = (min_pt[0], min_pt[1])
                    temp_pt3=(min_pt[0], get_y_intersection(min_pt[0], hull[i], hull[j]))


            elif near_lt(-pi,ang) and near_lt(ang,-pi/2):
                print('-pi < ang < -pi/2')
                
                min_pt=get_best_pt(x_keys, x_points,hull[i],hull[j],False,True,ang)
                print '\tmin_pt={0}'.format(min_pt)
                if near_ge(min_pt[0],hull[i][0]):
                    temp_pt1=(hull[i][0],min_pt[1])
                    temp_pt2 = (get_x_intersection(min_pt[1],hull[i],hull[j]),min_pt[1])
                else:
                    temp_pt1 = (min_pt[0], get_y_intersection(min_pt[0], hull[i],hull[j]))
                    temp_pt2 = (min_pt[0], min_pt[1])
                    temp_pt3 = (get_x_intersection(min_pt[1],hull[i],hull[j]),min_pt[1])

                
            elif near_lt(pi/2,ang) and near_lt(ang,pi):
                continue
                temp_y=get_max_y(x_keys, x_points,hull[i],hull[j])
                temp_pt2=(hull[j][0],temp_y)
                temp_pt1 = (get_x_intersection(temp_y,hull[i],hull[j]),temp_y)

            else:
                continue

            temp_pt1 = temp_pt1# if is_a else temp_ptb1
            temp_pt2 = temp_pt2# if is_a else temp_ptb2
            # temp_pt3 = None if
            # if temp_pta3 is not None:# and temp_ptb3 is not None:
            #     temp_pt3 = temp_pta3 #if is_a else temp_ptb3
            if length(temp_pt1,hull[i])>0 and length(temp_pt2, hull[j])>0:
                ret_hull.append(temp_pt1)
                ret_hull.append(temp_pt2)
                if temp_pt3 is not None:
                    ret_hull.append(temp_pt3)

        return ret_hull
    
    pt0 = min(points, key = lambda pt: pt[0])
    points.remove(pt0)
#    pt_array = [pt0]
    points = sorted(points, key = lambda pt: atan2(pt[1] - pt0[1],pt[0] - pt0[0] ))
    points.insert(0,pt0)

    # = reduce(_keep_left, points, [])
    
    l = reduce(_keep_left, points, [])
#    print('len(l)={0}'.format(len(l)))

#    print_angles(l)
   # pt0 = min(points, key = lambda pt: pt[0])
   # points.remove(pt0)
#    pt_array = [pt0]
    points = sorted(points, key = lambda pt: pt[0])
    big_hull_list = [p.exterior.coords[:-1] for p in big_poly_list]
    ret_hull = _add_hull(points,l, big_hull_list)
    print('len(points)={0}'.format(len(points)))
    return ret_hull
    #return l.extend(u[i] for i in range(1, len(u) - 1)) or l


def ccw(p1, p2, p3):
    return (p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1] - p1[1])*(p3[0] - p1[0])




def pt_on_line(pt, slope, intercept):
    if abs(slope)==float('inf'):
        # Check if x coord of pt is good
        return (pt.x == intercept)
    else:
        return (pt.y == pt.x*slope+intercept)


def pt_colinear(target_pt, pt1, pt2):
    if pt1.x == pt2.x:
        return pt_on_line(target_pt,float('inf'),pt1.x)
    else:
        slope = (pt2.y-pt1.y)/(pt2.x-pt1.x)
        intercept = pt2.y - slope*pt2.x
        return pt_on_line(target_pt,slope,intercept)

