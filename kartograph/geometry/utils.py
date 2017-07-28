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

    # temp_STATEFP=sidelayer.features[0].props['STATEFP']
    # temp_feat=create_feature(best_my_box,{'NAME': 'Side Box', 'LSAD': '01', 'STATEFP': temp_STATEFP, 'PLACEFP': '00000'})
    # sidelayer.options['specialstyle']+='#countylayer[PLACEFP=00000] { fill: none; }\n'
    # sidelayer.features.append(temp_feat)

    # statelayer = the_map.layersById['statelayer']
    # temp_feat = create_feature(best_out_point.buffer(dist_param/4.), {'NAME': 'Place Point', 'LSAD': '01', 'STATEFP': temp_STATEFP, 'PLACEFP': '99999'})
    # temp_feat2 = create_feature(LineString([(best_pointA.x,best_pointA.y), (best_pointB.x,best_pointB.y)] ), {'NAME': 'Best Seg', 'LSAD': '01', 'STATEFP': temp_STATEFP, 'PLACEFP': '99998'})
    # statelayer.options['specialstyle']=''
    # statelayer.options['specialstyle']+='#statelayer[PLACEFP=99999]\n{\n\tfill: black;\n}\n'
    # statelayer.options['specialstyle']+='#statelayer[PLACEFP=99998]\n{\n\tstroke: red;\nstroke-width: 5px;\n}\n'
    # statelayer.features.append(temp_feat)
    # statelayer.features.append(temp_feat2)
    
    return x_offset, y_offset

def get_offset_coords_super_complex(mainbbox, sidebbox, main_geom, side_geom, position_factor, the_map):
    opts = the_map.options
    data = opts['bounds']['data']
    sidelayer = the_map.layersById[data['sidelayer']]
    m_hull = main_geom#geom_to_bbox(main_geom,min_area=0)
    s_hull = side_geom#geom_to_bbox(side_geom,min_area=0)
    dist_param = min(mainbbox.width, mainbbox.height)/25.

    mainbbox = geom_to_bbox(main_geom,data['min-area'])

    best_x=0
    best_y=0
    min_area=float('inf')
    max_diff=0
    min_max_side=float('inf')
    
    #m_coords = m_hull.exterior.coords
    #s_coords = s_hull.exterior.coords

    hull_list=[]
    side_hull_list=[]
    main_perim=0
    polygons = isinstance(main_geom, MultiPolygon) and main_geom or [main_geom]
    for poly in polygons:
        m_coords = poly.exterior.coords
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
        #if seg.length < 0 * main_perim/len(hull_list):
        #    continue
        my_box = seg.get_box_pts(side_geom.convex_hull)
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
        my_feat = create_feature(deepcopy(side_geom.convex_hull),{})
        #print('my_feat.geometry={0}'.format(my_feat.geometry))
        my_feat.offset_feature(offset={'x': temp_x_offset, 'y': temp_y_offset})
        #print('my_feat.geometry={0}'.format(my_feat.geometry))
        print 'mainbbox={0}'.format(mainbbox)
        temp_bbox=deepcopy(mainbbox)
        temp_bbox.join(geom_to_bbox(my_feat.geometry,data['min-area']))
     
        print('now temp_bbox={0}'.format(temp_bbox))

        #print('hullseg={0}\n\tarea={1}\n\n'.format(seg,temp_bbox.area()))
        temp_max_side = max(temp_bbox.height,temp_bbox.width)
        if temp_max_side < min_max_side and my_feat.geometry.distance(main_geom) >= dist_param: #temp_bbox.area() < min_area:
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

    # temp_STATEFP=sidelayer.features[0].props['STATEFP']
    # temp_feat=create_feature(best_my_box,{'NAME': 'Side Box', 'LSAD': '01', 'STATEFP': temp_STATEFP, 'PLACEFP': '00000'})
    # sidelayer.options['specialstyle']+='#countylayer[PLACEFP=00000] { fill: none; }\n'
    # sidelayer.features.append(temp_feat)

    # statelayer = the_map.layersById['statelayer']
    # temp_feat = create_feature(best_out_point.buffer(dist_param/4.), {'NAME': 'Place Point', 'LSAD': '01', 'STATEFP': temp_STATEFP, 'PLACEFP': '99999'})
    # temp_feat2 = create_feature(LineString([(best_pointA.x,best_pointA.y), (best_pointB.x,best_pointB.y)] ), {'NAME': 'Best Seg', 'LSAD': '01', 'STATEFP': temp_STATEFP, 'PLACEFP': '99998'})
    # statelayer.options['specialstyle']=''
    # statelayer.options['specialstyle']+='#statelayer[PLACEFP=99999]\n{\n\tfill: black;\n}\n'
    # statelayer.options['specialstyle']+='#statelayer[PLACEFP=99998]\n{\n\tstroke: red;\nstroke-width: 5px;\n}\n'
    # statelayer.features.append(temp_feat)
    # statelayer.features.append(temp_feat2)
    
    return x_offset, y_offset



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

