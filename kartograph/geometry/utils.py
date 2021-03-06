"""
geometry utils
"""

from hullseg import hullseg
from feature import create_feature
from shapely.geometry import Point, Polygon, LineString, MultiPolygon
from math import atan, cos, sin, atan2, pi,sqrt, hypot, floor, ceil
from copy import deepcopy
from polypoint import PolyPoint

TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

MIN_DIFF=1e-4


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

# Get optimal side geometry offsets 

def get_offset_coords_super_complex(mainbbox, sidebbox, geom, side_geom, position_factor, the_map):
    opts = the_map.options
    data = opts['bounds']['data']
    sidelayer = the_map.layersById[data['sidelayer']]
    #m_hull = geom.convex_hull#geom_to_bbox(geom,min_area=0)

    the_dist = mainbbox.width * mainbbox.height / 500.

    print '\nMain'
    m_hull = get_complex_hull(mainbbox, sidebbox, geom, position_factor, the_map).exterior.coords[:-1]
    print '\nSide'
    s_hull = get_complex_hull(mainbbox, sidebbox, side_geom, position_factor, the_map).exterior.coords[:-1]
    m_out=[]
    s_out=[]
    m_ang=[]
    s_ang=[]
    for l in range(len(m_hull)):
        i=(l-1) % len(m_hull)
        j=l
        k=(l+1) % len(m_hull)
        ang_1 = atan2(round(m_hull[j][1]-m_hull[i][1],7),round(m_hull[j][0]-m_hull[i][0],7))
        ang_2 = atan2(round(m_hull[k][1]-m_hull[j][1],7),round(m_hull[k][0]-m_hull[j][0],7))
        turn_ang = turn_angle(m_hull[i],m_hull[j],m_hull[k])/pi

        sign_mult = 1 if turn_ang > 0 and turn_ang < 1 else -1
        # print 'm[{2}]: ang_1={0:.5f}, ang_2={1:.5f}, ang_3={3:.5f}'.format(ang_1/pi, ang_2/pi,l,turn_ang)
        next_pt = (m_hull[j][0] + sign_mult*the_dist*(cos(ang_1)-cos(ang_2)), m_hull[j][1] + sign_mult*the_dist*(sin(ang_1)-sin(ang_2)))
        
        m_ang.append(atan2(next_pt[1]-m_hull[j][1],next_pt[0]-m_hull[j][0]))
        m_out.append(next_pt)
    for l in range(len(s_hull)):
        i=(l-1) % len(s_hull)
        j=l
        k=(l+1) % len(s_hull)
        ang_1 = atan2(s_hull[j][1]-s_hull[i][1],s_hull[j][0]-s_hull[i][0])
        ang_2 = atan2(s_hull[k][1]-s_hull[j][1],s_hull[k][0]-s_hull[j][0])
        turn_ang = turn_angle(s_hull[i],s_hull[j],s_hull[k])/pi

        sign_mult = 1 if turn_ang > 0 and turn_ang < 1 else -1
        # print 's[{2}]: ang_1={0:.5f}, ang_2={1:.5f}'.format(ang_1/pi, ang_2/pi,l) 

        next_pt = (s_hull[j][0] + sign_mult*the_dist*(cos(ang_1)-cos(ang_2)), s_hull[j][1] + sign_mult*the_dist*(sin(ang_1)-sin(ang_2)))
        s_ang.append(atan2(next_pt[1]-s_hull[j][1],next_pt[0]-s_hull[j][0]))
        
        s_out.append(next_pt)

    # print_hull(m_hull)
    # print '\n'
    # print_hull(m_out)
    # print 'm_out={0}, s_out={1}'.format(m_out, s_out)
    new_x = m_out[1][0]
    new_y = m_out[1][1]
    min_length=float('inf')
    if near_eq(ang_diff(m_ang[1],s_ang[0]),pi):
        offset_x = new_x - s_hull[0][0]
        offset_y = new_y - s_hull[0][1]
    else:
        offset_x=0
        offset_y=0
    best_i=0
    best_j=0
    min_length=float('inf')
    for i,j in [(i,j) for i in range(len(s_hull)) for j in range(len(m_hull))]:
        temp_new_x = m_out[j][0]
        temp_new_y = m_out[j][1]
        tempoffset_x = temp_new_x - s_hull[i][0]
        tempoffset_y = temp_new_y - s_hull[i][1]
        
        the_diff=ang_diff(m_ang[j],s_ang[i])
        temp_min_max_length = min_max_length(m_hull, s_hull, tempoffset_x, tempoffset_y)
        if temp_min_max_length < min_length and adequate_spacing(m_hull, s_hull, tempoffset_x, tempoffset_y, the_dist/4.):
            
            # print('m_hull[{5}]={0}, s_hull[{1}]={2}, m_ang[{5}]={3:.4f}, s_ang[{1}]={4:.4f}, len={6}'.format(m_hull[j],i,s_hull[i],m_ang[j]/pi,s_ang[i]/pi,j, temp_min_max_length))
            offset_x = tempoffset_x
            offset_y = tempoffset_y
            min_length = temp_min_max_length

            best_i=i
            best_j=j

    s_hull=[(s_hull[i][0],s_hull[i][1]) for i in range(len(s_hull))]

    
    return Polygon(m_hull), Polygon(s_hull), offset_x, offset_y,LineString([m_hull[best_j],(s_hull[best_i][0]+offset_x, s_hull[best_i][1]+offset_y)])

''' Return the length of the larger side of the rectangular region containing m_hull and the
    offset by tempoffset_x, tempoffset_y of s_hull'''
def min_max_length(m_hull, s_hull, tempoffset_x, tempoffset_y):
    maxx=max([pt[0] for pt in m_hull]+[pt[0]+tempoffset_x for pt in s_hull])
    minx=min([pt[0] for pt in m_hull]+[pt[0]+tempoffset_x for pt in s_hull])
    maxy=max([pt[1] for pt in m_hull]+[pt[1]+tempoffset_y for pt in s_hull])
    miny=min([pt[1] for pt in m_hull]+[pt[1]+tempoffset_y for pt in s_hull])

    return max((maxx-minx),(maxy-miny))

def adequate_spacing(m_hull, s_hull, tempoffset_x, tempoffset_y, the_dist):
    s_temp = list(map(lambda x: (x[0]+tempoffset_x, x[1]+tempoffset_y), s_hull))
    poly1 = Polygon(m_hull)
    poly2 = Polygon(s_temp)
    return poly1.distance(poly2) >= the_dist
    

    
    
    
def print_hull(the_hull):
    to_print='({0:.5f}, {1:.5f})'.format(the_hull[0][0], the_hull[0][1])

    for i in range(1,len(the_hull)):
        to_print+=',({0:.5f},{1:.5f})'.format(the_hull[i][0], the_hull[i][1])

    print '{0}'.format(to_print)

def ang_diff(ang1, ang2):
    my_diff = ang1 - ang2
    if my_diff < 0:
        my_diff = my_diff + 2*pi
    return my_diff


def get_complex_hull(mainbbox, sidebbox, geom, position_factor, the_map):
    opts = the_map.options
    data = opts['bounds']['data']
    sidelayer = the_map.layersById[data['sidelayer']]
    m_hull = geom.convex_hull#geom_to_bbox(geom,min_area=0)
 #   s_hulls = side_geom#geom_to_bbox(side_geom,min_area=0)
    max_area=0
    best_poly = None
    coord_list=[]
    seg_list=[]
    polypoint_list=[]
    # print('Getting complex hull')

    ret_poly_list=[]

    if not isinstance(geom, MultiPolygon):
        poly_list = [geom]
    else:
        poly_list = geom
    max_area = max(poly.area for poly in poly_list)
    big_poly_list = [poly for poly in poly_list]# if poly.area*50000 >= max_area]
    next_poly_list=[]
    
    for poly in big_poly_list:
        temp_poly=poly.exterior.coords[:-1]
        #coord_list=[]
        for i in range(len(temp_poly)):
            j=(i+1) % len(temp_poly)
            k=(i-1) % len(temp_poly)
            temp_pt=PolyPoint(temp_poly[i][0],temp_poly[i][1],temp_poly[j][0],temp_poly[j][1],temp_poly[k][0],temp_poly[k][1])
            coord_list.append(temp_pt)
        
    #curr_poly = Polygon(coord_list).convex_hull
    #ret_poly_list.append(curr_poly)

    curr_poly = Polygon(convex_hull_jacob(coord_list, big_poly_list))

#    return geom.convex_hull
    return curr_poly
    # print 'Bounds(curr_poly)={0}'.format(curr_poly.bounds)
    # ret_poly_list.append(curr_poly)
        

    # if len(ret_poly_list)>1:
    #     return MultiPolygon(ret_poly_list)
    # else:
    #     return ret_poly_list[0]

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

    
def convex_hull_graham(points):
    '''
    Returns points on convex hull in CCW order according to Graham's scan algorithm. 
    By Tom Switzer <thomas.switzer@gmail.com>.
    '''
    TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

    MIN_DIFF=1e-7

    points = sorted(points)
    l = reduce(_keep_left, points, [])
    u = reduce(_keep_left, reversed(points), [])
    return l.extend(u[i] for i in range(1, len(u) - 1)) or l

def turn_angle(p,q,r):
    ang = atan2(r[1]-q[1],r[0]-q[0])-atan2(q[1]-p[1],q[0]-p[0])
    ang = ang + 2*pi if ang <= -pi else ang
    return ang

def length(a,b):
    return sqrt((a[1]-b[1])**2+(a[0]-b[0])**2)


def convex_hull_jacob(points, big_poly_list):
    '''
    Returns points on convex hull in CCW order according to Graham's scan algorithm. 
    By Jacob Alperin-Sheriff
    '''
    from math import atan2, sqrt, pi
    from bisect import bisect_left, bisect, bisect_right
    
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
        # print('extremum={0}'.format(extremum))
        the_iter=range(lo,hi+1,1) if bool(is_left) and not bool(is_above) else range(hi,lo-1,-1)
        good_pts=[]
        cross_bool = bool(not is_left) ^ bool(is_above)
        for i in the_iter:
            if (is_above and to_search[i][this_coord] < extremum) or (not is_above and to_search[i][this_coord] > extremum):
                good_pts.append((to_search[i][0], to_search[i][1]))
                extremum = to_search[i][this_coord]
            # elif (not is_above and to_search[i][this_coord] < extremum):
            #    good_pts=[pt for pt in good_pts if pt[1] > to_search[i][this_coord]]
                    
        # the_iter=range(lo,hi+1,1)#if bool(not is_left) ^ bool(is_above) else range(hi,lo-1,-1)
        if not is_above and not is_left:
            print('ptA={0}, ptB={1}, good_pts={2}'.format(ptA,ptB,good_pts))
        extremum = other_func(ptA[this_coord],ptB[this_coord])           
        for curr_pt in good_pts:

            temp_eval = pt_eval(ptA,ptB,curr_pt)

            if temp_eval > best_eval:
                best_pt = (curr_pt[0], curr_pt[1])
                best_eval = temp_eval
        # print('(is_above={0},is_left={1}, ang={4}), best_pt={2}, temp_eval={3}'.format(is_above,is_left,best_pt,temp_eval, ang/pi))

            extremum = curr_pt[this_coord]
        # print('')
        if best_pt[0] == ptB[0] and best_pt[1] == ptB[1]:
            best_pt = (ptA[0],ptA[1])
        return best_pt

    ''' Note: keys will be the x points sorted by x low to high, to_search will be sorted by x 
as well '''
    
    def get_best_pt_crosssec(keys, to_search, ptA, ptB, is_above, is_left, ang):
        num_divs=8
        other_coord=0
        cross_bool = bool(not is_left) ^ bool(is_above)
        x_sign = 1 if is_left else -1
        y_sign = 1 if is_above else -1
        z_sign = x_sign
        if not cross_bool:
            inner_pt = (ptA[0], ptB[1])
            outer_pt = (ptB[0], ptA[1])
        else:
            inner_pt = (ptB[0], ptA[1])
            outer_pt = (ptA[0], ptB[1])
        width = abs(inner_pt[0] - outer_pt[0])
        height = abs(inner_pt[1] - outer_pt[1])
        is_full = {(i,j) : False for i in range(num_divs) for j in range(num_divs)}
        lo = find_gt(keys,min(ptA[other_coord],ptB[other_coord]))
        hi= find_lt(keys,max(ptA[other_coord],ptB[other_coord]))
        the_iter=range(lo,hi+1,1)
        map_x = lambda x: int(floor(num_divs*abs(x-outer_pt[0])/width))
        map_y = lambda y: int(floor(num_divs*abs(y-outer_pt[1])/height))
        
        map_pos = lambda x, y: (map_x(x),map_y(y))
        pt_ct=0
        set_lines=[]
        set_linesB=[]
        for i in the_iter:
            pt_ct=pt_ct+1
            set_lines=[]
            set_linesB=[]
            # print 'to_search[i]={0}'.format(to_search[i])
            if not (to_search[i][1] >= min(ptA[1], ptB[1]) and to_search[i][1] <= max(ptA[1], ptB[1])):
                continue
            (x,y)=map_pos(to_search[i][0],to_search[i][1])
            if x < num_divs and y < num_divs and 0 <= x and 0<=y:
                is_full[(x,y)]=True

                (prev_x,prev_y)=map_pos(to_search[i][2],to_search[i][3])
                (next_x,next_y)=map_pos(to_search[i][4],to_search[i][5])
                # if is_left and is_above:
                #     print '(prev_x,prev_y)={0}, (next_x, next_y)={1}'.format((prev_x,prev_y),(next_x,next_y))
                if to_search[i][0]==to_search[i][2]:
                    for temp_y in range(min(y,prev_y),max(y,prev_y)+1):
                        is_full[(x,temp_y)]=True
                else:
                    slope1,int1 = slope_intercept((to_search[i][0],to_search[i][1]),(to_search[i][2],to_search[i][3]))
                    t_func = lambda x: slope1*x+int1
                    for temp_x in range(max(0,min(x,prev_x)),min(num_divs,max(x,prev_x)+1)):
                        ty1 = map_y(t_func(outer_pt[0]+x_sign*temp_x*(width/num_divs)))
                        ty2 = map_y(t_func(outer_pt[0]+x_sign*(temp_x+1.0)*(width/num_divs)))
                        t_miny=min(ty1,ty2)
                        t_maxy=max(ty1,ty2)
                        for temp_y in range(max(t_miny,0), min(t_maxy+1,num_divs)):
                            # print '({0},{1}) set to True'.format(temp_x,temp_y)
                            is_full[(temp_x,temp_y)] = True
                            set_lines.append((temp_x,temp_y))
           

                if to_search[i][0]==to_search[i][4]:
                    for temp_y in range(min(y,next_y),max(y,next_y)+1):
                        is_full[(x,temp_y)]=True
                else:
                    slope1,int1 = slope_intercept((to_search[i][0],to_search[i][1]),(to_search[i][4],to_search[i][5]))

                    t_func = lambda x: slope1*x+int1
                    for temp_x in range(max(0,min(x,next_x)),min(num_divs,max(x,next_x)+1)):
                        ty1 = map_y(t_func(outer_pt[0]+x_sign*temp_x*(width/num_divs)))
                        ty2 = map_y(t_func(outer_pt[0]+x_sign*(temp_x+1.0)*(width/num_divs)))
                        t_miny=min(ty1,ty2)
                        t_maxy=max(ty1,ty2)
                        for temp_y in range(max(t_miny,0), min(t_maxy+1,num_divs)):
                            # print '({0},{1}) set to True'.format(temp_x,temp_y)
                            is_full[(temp_x,temp_y)] = True
                            set_linesB.append((temp_x,temp_y))
                # if is_left and is_above:
                #     print('set_lines={0}\tset_linesB={1}'.format(set_lines, set_linesB))         
                    
            # Need to fill the gap 
        for i in range(num_divs):
            for j in range(num_divs):
                if is_full[(i,j)]:
                    for k in range(i+1,num_divs):
                        is_full[(k,j)]=True
                    for k in range(j+1,num_divs):
                        is_full[(i,k)]=True
        if pt_ct==0:
            print('no pts')
            best_pt=(outer_pt[0],outer_pt[1])
            print 'ptA={0}, ptB={1}, best_pt={2}, outer_pt={3}\n'.format(ptA,ptB,best_pt,outer_pt)
            return best_pt
        try:
            unfull=[(x,y) for (x,y) in is_full if not is_full[(x,y)]]
            full=[(x,y) for (x,y) in is_full if is_full[(x,y)]]
            (x,y) = max(unfull, key = lambda (x,y): x+y)
            best_pt = (outer_pt[0] + x_sign * width * (x+1) / num_divs, outer_pt[1] + y_sign * height * (y+1) / num_divs)
            if is_left and is_above:
                print 'full={0}'.format(full)

                print 'ptA={0}, ptB={1}, (x,y)={2}, best_pt={3}, outer_pt={4}\n'.format(ptA,ptB,(x,y),best_pt,outer_pt)

        except ValueError:
            best_pt=(outer_pt[0], outer_pt[1])
        # best_pt = (outer_pt[0], outer_pt[1])   
        return best_pt
    

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

    def slope_intercept(ptA,ptB):
        if ptA[0]==ptB[0]:
            raise ValueError
        slope = (ptB[1]-ptA[1])/(ptB[0]-ptA[0])
        intercept = ptB[1]-slope*ptB[0]
        return slope, intercept

    
    def _add_hull(points, hull):
        # Need to check if two points on convex hull are on same simple polygon or not??
        ret_hull=[]
        x_points = points
        # y_points = sorted(deepcopy(points), key = lambda pt: pt[1])
        
        x_keys = [pt[0] for pt in points]
        # y_keys = [pt[1] for pt in y_points]
        for i in range(0,len(hull)):
            j=(i+1) % len(hull)
            ang = _line_angle(hull[i],hull[j])
            ret_hull.append(hull[i])
            temp_pt1 = temp_pt2 = temp_pt3 = None
            print('hull[{0}]={1}, hull[{2}]={3}'.format(i,hull[i],j,hull[j]))
            if near_lt(0,ang) and near_lt(ang,pi/2):
                # print('0 < ang < pi/2')
                #continue

                min_pt=get_best_pt_crosssec(x_keys, x_points,hull[i],hull[j],True,False,ang)
                # print '\tmin_pt={0}'.format(min_pt)
                #continue
                if near_ge(min_pt[0],max(hull[i][0],hull[j][0])):
                    # temp_pt1=(hull[i][0],min_pt[1])
                    # temp_pt2 = (get_x_intersection(min_pt[1],hull[i],hull[j]),min_pt[1])
                    temp_pt1=(hull[i][0],min_pt[1])
                    temp_pt2 = (hull[j][0],min_pt[1])

                else:
                    # temp_pt1 = (min_pt[0], get_y_intersection(min_pt[0], hull[i],hull[j]))
                    # temp_pt2 = (min_pt[0], min_pt[1])
                    # temp_pt3 = (get_x_intersection(min_pt[1],hull[i],hull[j]),min_pt[1])
                    temp_pt1 = (min_pt[0], hull[i][1])
                    temp_pt2 = (min_pt[0], min_pt[1])
                    temp_pt3 = (hull[j][0],min_pt[1])

            elif near_lt(-pi/2,ang) and near_lt(ang,0):
                # print('-pi/2 < ang < 0')
                #continue
                min_pt=get_best_pt_crosssec(x_keys, x_points,hull[i],hull[j],True,True,ang)
                print '\tmin_pt={0}'.format(min_pt)

                temp_pt1 = (hull[i][0],min_pt[1])
                temp_pt2 = (min_pt[0], min_pt[1])
                temp_pt3=(min_pt[0], hull[j][1])

                # if near_gt(min_pt[0],max(hull[i][0],hull[j][0])):
                #     # temp_pt1 = (get_x_intersection(min_pt[1],hull[i],hull[j]),min_pt[1])
                #     # temp_pt2=(hull[j][0],min_pt[1])
                #     temp_pt1 = (hull[i][0],min_pt[1])
                #     temp_pt2=(hull[j][0],min_pt[1])

                # else:
                #     temp_pt1 = (hull[i][0],min_pt[1])
                #     temp_pt2 = (min_pt[0], min_pt[1])
                #     temp_pt3=(min_pt[0], hull[j][1])


            elif near_lt(-pi,ang) and near_lt(ang,-pi/2):
                # print('-pi < ang < -pi/2')
                #continue
                min_pt=get_best_pt_crosssec(x_keys, x_points,hull[i],hull[j],False,True,ang)
                # print '\tmin_pt={0}'.format(min_pt)
                if near_ge(min_pt[0],max(hull[i][0],hull[j][0])):
                    temp_pt1=(hull[i][0],min_pt[1])
                    temp_pt2 = (hull[j][0],min_pt[1])
                else:
                    temp_pt1 = (min_pt[0], hull[i][1])
                    temp_pt2 = (min_pt[0], min_pt[1])
                    temp_pt3 = (hull[j][0],min_pt[1])

                
            elif near_lt(pi/2,ang) and near_lt(ang,pi):
                # print('pi/2 < ang < pi')
                min_pt=get_best_pt_crosssec(x_keys, x_points,hull[i],hull[j],False,False,ang)
                #continue

                # print '\tmin_pt={0}'.format(min_pt)
                if near_ge(min_pt[0],max(hull[i][0],hull[j][0])):
                    temp_pt1 = (hull[i][0],min_pt[1])
                    temp_pt2=(hull[j][0],min_pt[1])
                else:
                    temp_pt1 = (hull[i][0],min_pt[1])
                    temp_pt2 = (min_pt[0], min_pt[1])
                    temp_pt3=(min_pt[0], hull[j][1])

            else:
                continue
            # ret_hull.append(min_pt)
            temp_pt1 = temp_pt1# if is_a else temp_ptb1
            temp_pt2 = temp_pt2# if is_a else temp_ptb2
            # temp_pt3 = None if
            # if temp_pta3 is not None:# and temp_ptb3 is not None:
            #     temp_pt3 = temp_pta3 #if is_a else temp_ptb3
            if length(temp_pt1,hull[i])>0:
                ret_hull.append(temp_pt1)
            if length(temp_pt2,hull[j])>0:
                ret_hull.append(temp_pt2)
            if temp_pt3 is not None and length(temp_pt3,hull[j])>0:
                ret_hull.append(temp_pt3)
        return ret_hull

    
    pt0 = min(points, key = lambda pt: pt.x)
    points.remove(pt0)
#    pt_array = [pt0]
    points = sorted(points, key = lambda pt: atan2(pt.y - pt0.y,pt.x - pt0.x ))
    points.insert(0,pt0)

    # = reduce(_keep_left, points, [])
    
    l = reduce(_keep_left, points, [])
#    print('len(l)={0}'.format(len(l)))

#    print_angles(l)
   # pt0 = min(points, key = lambda pt: pt[0])
   # points.remove(pt0)
#    pt_array = [pt0]
    points = sorted(points, key = lambda pt: pt.x)
    print('len(l) before add_hull={0}'.format(len(l)))
    ret_hull = _add_hull(points,l)
    print('len(hull) after add_hull={0}'.format(len(ret_hull)))
    return ret_hull
#return l.extend(u[i] for i in range(1, len(u) - 1)) or l

def _line_angle(a,b):
    return atan2(b[1]-a[1],b[0]-a[0])

    


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

