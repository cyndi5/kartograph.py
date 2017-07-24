"""
geometry utils
"""

from hullseg import hullseg
from feature import create_feature
from shapely.geometry import Point, Polygon
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
    best_my_box = None
    for seg in hull_list:
        # Find the rotated square bounding box with slope of some side according to seg
        my_box = get_box_pts(seg.slope, s_hull)
        sidelayer = the_map.layersById[data['sidelayer']]
        box_coords = my_box.exterior.coords
        # coordinate order is left top right bottom
        if seg.above and seg.slope >= 0:
            x_offset = seg.outPoint.x - (box_coords[0][0]+box_coords[3][0])/2.
            y_offset = seg.outPoint.y - (box_coords[0][1]+box_coords[3][1])/2.
        elif seg.above and seg.slope < 0:
            x_offset = seg.outPoint.x - (box_coords[2][0]+box_coords[3][0])/2.
            y_offset = seg.outPoint.y - (box_coords[2][1]+box_coords[3][1])/2.
        elif not seg.above and seg.slope >= 0:
            x_offset = seg.outPoint.x - (box_coords[1][0]+box_coords[2][0])/2.
            y_offset = seg.outPoint.y - (box_coords[1][1]+box_coords[2][1])/2.
        else:
            x_offset = seg.outPoint.x - (box_coords[0][0]+box_coords[1][0])/2.
            y_offset = seg.outPoint.y - (box_coords[0][1]+box_coords[1][1])/2.

    # Deal with checking how good it is offset 
  

    return x_offset, y_offset

# get the possibly rotated box holding the hull with a side parallel 
def get_box_pts(slope, hull):
    hull_coords = hull.exterior.coords
    if slope == 0 or abs(slope) == float('inf'):
        temp_bbox=geom_to_bbox(hull.geometry)
        return Polygon([(temp_bbox.left,temp_bbox.top),(temp_bbox.right,temp_bbox.top),(temp_bbox.right,temp_bbox.bottom),(temp_bbox.left,temp_bbox.bottom),(temp_bbox.left,temp_bbox.top)])
    else:
        inv_slope = -1./slope
        min_intcpt = float('inf')
        min_inv_intcpt = float('inf')
        max_intcpt = -1*float('inf')
        max_inv_intcpt = -1*float('inf')
        left_pt=Point(hull_coords[0][0],hull_coords[0][1])
        right_pt=Point(hull_coords[0][0],hull_coords[0][1])
        top_pt=Point(hull_coords[0][0],hull_coords[0][1])
        bottom_pt=Point(hull_coords[0][0],hull_coords[0][1])
        for coord in hull_coords:
            curr_intcpt = coord[1] - slope*coord[0]
            curr_inv_intcpt = coord[1] - inv_slope * coord[0]
            if curr_intcpt < min_intcpt:
                right_pt = Point(coord[0],coord[1])
                min_intcpt = curr_intcpt
            if curr_intcpt > max_intcpt:
                left_pt = Point(coord[0],coord[1])
                max_intcpt = curr_intcpt
            if curr_inv_intcpt < min_inv_intcpt:
                bottom_pt = Point(coord[0],coord[1])
                min_inv_intcpt = curr_inv_intcpt
            if curr_inv_intcpt > max_inv_intcpt:
                top_pt = Point(coord[0],coord[1])
                max_intcpt = curr_inv_intcpt
        return Polygon([(left_pt.x,left_pt.y),(top_pt.x,top_pt.y),(right_pt.x,right_pt.y),(bottom_pt.x,bottom_pt.y),(left_pt.x,left_pt.y)])

        
    
        

        
            
