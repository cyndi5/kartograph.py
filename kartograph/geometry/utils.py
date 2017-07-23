"""
geometry utils
"""

from hullseg import hullseg
from shapely.geometry import Point
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
def get_offset_coords_complex(mainbbox, sidebbox, main_geom, side_geom, position_factor):
    m_hull = main_geom#geom_to_bbox(main_geom,min_area=0)
    s_hull = side_geom#geom_to_bbox(side_geom,min_area=0)
    dist_param = min(mainbbox.width, mainbbox.height)/10.

    m_coords = m_hull.exterior.coords
    s_coords = s_hull.exterior.coords

    best_x=0
    best_y=0
    min_area=float('inf')
    max_diff=0
    hull_list=[]
    side_hull_list=[]
    for i in range(len(m_coords)-1):
        pA=Point(m_coords[i % len(m_coords)][0], m_coords[i % len(m_coords)][1])
        if i<len(m_coords)-2:
            b_pt_val=i+1
            o_pt_val=i+2
        else:
            b_pt_val=i+1
            o_pt_val=i+3
        pB=Point(m_coords[(b_pt_val) % len(m_coords)][0], m_coords[(b_pt_val) % len(m_coords)][1])
        pO=Point(m_coords[(o_pt_val) % len(m_coords)][0], m_coords[(o_pt_val) % len(m_coords)][1])
        #print('{0}, {1}, {2}'.format(pA,pB,pO))
        temp_hull_seg=hullseg(pA,pB,pO,distParam=dist_param)
        print('current segment-{0}'.format(temp_hull_seg))
        hull_list.append(temp_hull_seg)

    for i in range(len(s_coords)-1):
        pA=Point(s_coords[i % len(s_coords)][0], s_coords[i % len(m_coords)][1])
        if i<len(s_coords)-2:
            b_pt_val=i+1
            o_pt_val=i+2
        else:
            b_pt_val=i+1
            o_pt_val=i+3
        pB=Point(s_coords[(b_pt_val) % len(s_coords)][0], s_coords[(b_pt_val) % len(s_coords)][1])
        pO=Point(s_coords[(o_pt_val) % len(s_coords)][0], s_coords[(o_pt_val) % len(s_coords)][1])
        #print('{0}, {1}, {2}'.format(pA,pB,pO))
        temp_hull_seg=hullseg(pA,pB,pO,distParam=dist_param, point_pos=i)
        print('  current segment-{0}'.format(temp_hull_seg))
        side_hull_list.append(temp_hull_seg)

    for seg in hull_list:
        # Find the rotated square box with slope of some side according to seg
         
        
        # for side_seg in side_hull_list:
        #     # check_point should tell if line through s_point
        #     # with slope of seg is internal to s_hull or not
        #     is_a = seg.check_point(side_seg.point_pos, m_coords, s_coords)
        #     is_b = seg.check_point(side_seg.point_pos+1, m_coords, s_coords)

        #     #is_a = True    
        #     if is_a:# or is_b:
        #         # May went to do something with segments on side too
        #         # Fix later to do the annoying checking to see where it should lie
        #         x_offset, y_offset = get_curr_offsets(seg, side_seg, is_a)
        #         temp_sbox = sidebbox.get_offset_box(x_offset, y_offset)
        #         temp_join_bbox = deepcopy(mainbbox)
        #         #print('temp_join_bbox={0}').format(temp_join_bbox)
        #         temp_join_bbox.join(temp_sbox)
        #         temp_diff = min(temp_join_bbox.width/temp_join_bbox.height,
        #                         temp_join_bbox.height/temp_join_bbox.width)
        #         if seg.slope == side_seg.slope:
        #             temp_diff=0
        #         elif seg.slope == float('nan') or side_seg.slope==float('nan'):
        #             temp_diff=-1*float('inf')
        #         else:
        #             temp_diff = abs(atan(seg.slope)-atan(side_seg.slope))
        #         if temp_diff < max_diff:
        #             max_diff = temp_diff
        #             #min_area=temp_join_bbox.area()
        #             best_x = x_offset
        #             best_y = y_offset
        #             print('{0}, max_diff={1}'.format(seg, max_diff))

    return x_offset, y_offset

# get the desired offsets         
def get_curr_offsets(seg, side_seg, is_a):
    sideOutPoint = Point(0,0)
    if side_seg.inv_slope == float('nan'):
        if seg.slope == float('nan'):
            sideOutPoint = deepcopy(side_seg.midPoint)
        else:
            if is_a:
                next_intercept = side_seg.pointA.y - seg.slope*side_seg.pointA.x
            else:
                next_intercept = side_seg.pointB.y - seg.slope*side_seg.pointB.y            
            sideOutPoint = Point(side_seg.inv_intercept,seg.slope*side_seg.inv_intercept + next_intercept)
            return sideOutPoint.x, sideOutPoint.y
    else:
        if seg.slope == float('nan'):
            if is_a:
                next_intercept = side_seg.pointA.x
            else:
                next_intercept = side_seg.pointB.x
            sideOutPoint = Point(next_intercept, side_seg.slope*next_intercept + side_seg.inv_intercept)
            # leave for now
        else:
            if is_a:
                next_intercept = side_seg.pointA.y - seg.slope*side_seg.pointA.x
            else:
                next_intercept = side_seg.pointB.y - seg.slope*side_seg.pointB.x
            sideOutPoint = deepcopy(side_seg.midPoint)
            #sideOutPoint = Point((-next_intercept - side_seg.inv_intercept*1.)/(seg.slope + side_seg.inv_slope), sideOutPoint.x*seg.slope+next_intercept)
    #return -side_seg.pointA.x+seg.outPoint.x, -side_seg.pointA.y+seg.outPoint.y
    return -sideOutPoint.x+seg.outPoint.x, -sideOutPoint.y+seg.outPoint.y
        


    
        

        
            
