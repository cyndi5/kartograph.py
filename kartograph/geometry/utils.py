"""
geometry utils
"""

from hullseg import hullseg
from shapely.geometry import Point
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
        temp_hull_seg=hullseg(pA,pB,pO,dist_param)
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
        temp_hull_seg=hullseg(pA,pB,pO,dist_param)
        print('current segment-{0}'.format(temp_hull_seg))
        side_hull_list.append(temp_hull_seg)

    for seg in hull_list:
        for side_seg in side_hull_list:
            # check_point should tell if line through s_point
            # with slope of seg is internal to s_hull or not 
            if seg.check_point(side_seg.pointA, m_coords, s_coords):
                # May went to do something with segments on side too
                # Fix later to do the annoying checking to see where it should lie
                x_offset = (-s_coords[s_point_pos][0]+seg.outPoint.x)
                y_offset = (-s_coords[s_point_pos][1]+seg.outPoint.y)
                temp_sbox = sidebbox.get_offset_box(x_offset, y_offset)
                temp_join_bbox = deepcopy(mainbbox)
                #print('temp_join_bbox={0}').format(temp_join_bbox)
                temp_join_bbox.join(temp_sbox)
                #print('Post join, temp_join_bbox={0}').format(temp_join_bbox)
            #    print('temp_join_bbox.area={0},min_area={1}'.format(temp_join_bbox.area(),min_area))
                if temp_join_bbox.area() < min_area:
                    min_area=temp_join_bbox.area()
                    best_x = x_offset
                    best_y = y_offset
                    print('{0}, {1},min_area={2}'.format(seg, s_coords[s_point_pos],min_area))
            elif seg.check_point(side_seg.pointB, m_coords, s_coords):
                # May went to do something with segments on side too
                x_offset = (-s_coords[s_point_pos][0]+seg.outPoint.x)
                y_offset = (-s_coords[s_point_pos][1]+seg.outPoint.y)
                temp_sbox = sidebbox.get_offset_box(x_offset, y_offset)
                temp_join_bbox = deepcopy(mainbbox)
                #print('temp_join_bbox={0}').format(temp_join_bbox)
                temp_join_bbox.join(temp_sbox)
                #print('Post join, temp_join_bbox={0}').format(temp_join_bbox)
            #    print('temp_join_bbox.area={0},min_area={1}'.format(temp_join_bbox.area(),min_area))
                if temp_join_bbox.area() < min_area:
                    min_area=temp_join_bbox.area()
                    best_x = x_offset
                    best_y = y_offset
                    print('{0}, {1},min_area={2}'.format(seg, s_coords[s_point_pos],min_area))
    return x_offset, y_offset
            
        
    

    
        

        
            
