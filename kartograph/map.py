from shapely.geometry import Polygon, MultiPolygon
from shapely.geometry.base import BaseGeometry
from maplayer import MapLayer
from geometry.utils import join_features
from geometry import create_feature
from geometry.utils import geom_to_bbox
from geometry.utils import bbox_to_polygon
from geometry.feature import MultiPolygonFeature
from math import sqrt


from copy import deepcopy
from geometry import BBox, View
from proj import projections
from filter import filter_record
from errors import KartographError
import sys
import re
# Map
# ---
#
# This class performs like 80% of the functionality of Kartograph. It
# loads the features for each layer, processes them and passes them
# to a renderer at the end.

verbose = False


class Map(object):

    def __init__(me, options, layerCache, format='svg', src_encoding=None):
        me.options = options
        me.format = format
#        print 'map.init : me.options={0}'.format(me.options)
        # List and dictionary references to the map layers.
        me.layers = []
        me.layersById = {}
        # We will cache the bounding geometry since we need it twice, eventually.
        me._bounding_geometry_cache = False
        me._unprojected_bounds = None
        me._side_bounding_geometry = None
        me._projected_bounds = None
        me._side_offset = {'x':0,'y':0}
        me._next_side_offset={'x':0,'y':0}
        me._init_offset = None # the initial offset as a result of scaling
        # The **source encoding** will be used as first guess when Kartograph tries to decode
        # the meta data of shapefiles etc. We use Unicode as default source encoding.
        if not src_encoding:
            src_encoding = 'utf-8'
        me._source_encoding = src_encoding

        # Construct [MapLayer](maplayer.py) instances  and store references
        # to the layers in a list and a dictionary. 
        # I wish I remember why I split it, should've commented it ...
        for layer_cfg in options['layers']:
            #print 'layer_cfg={0}\n'.format(layer_cfg)
            layer_id = layer_cfg['id']
            #if layer_id=='statelayer':
            layer = MapLayer(layer_id, layer_cfg, me, layerCache)
            me.layers.append(layer)
            me.layersById[layer_id] = layer

        layer = me.layersById['placelayer']
        layerMainFilter = lambda rec: filter_record(layer.options['main-filter'], rec)
        me._main_geom = layer.source.get_main_geom(main_filter=layerMainFilter)
        #print 'main_geom={0}'.format(me._main_geom)
        # Initialize the projection that will be used in this map. This sounds easier than
        
        # it is since we need to compute lot's of stuff here.
        #print 'init projection'
        me.proj = me._init_projection()
        me.side_proj = me._init_side_projection()
        #print '**me.side_proj={0}'.format(me.side_proj)
        # Compute the bounding geometry for the map.
        me.bounds_poly = me._init_bounds()
        #print '**init bounds, me._projected_bounds={0}, me.bounds_poly={1}'.format(me._projected_bounds, me.bounds_poly)
    
    #       First, load the main_geometry of the place we want
      

        # Load all features that could be visible in each layer. The feature geometries will
        # be projected 

        print 'me._side_bounding_geometry.hash={0}'.format(hash(str(me._side_bounding_geometry)))
        for layer in me.layers:
            if "sidelayer" in layer.options:
                layer.options["init_offset"]=(0,0) #me._init_offset
               
            else:
                layer.options["init_offset"]=(0, 0)
#            print "getting features for layer={0}".format(layer)
            if layer.id==me.options['bounds']['data']['sidelayer']:
                layer.get_features(contained_geom=me._main_geom)
                #me._side_bounding_geometry=layer.features[0].geometry
               # print 'side bounding geometry bbox={0}'.format(geom_to_bbox(me._side_bounding_geometry))
                print "done getting features for layer={0}".format(layer.id)

       
        print 'Next me._side_bounding_geometry.hash={0}'.format(hash(str(me._side_bounding_geometry)))
        for layer in me.layers:
            if "sidelayer" in layer.options:
                layer.options["init_offset"]=(0,0) #me._init_offset
                #print 'layer.id={0}, init_offset={1}'.format(layer.id,me._init_offset)
            else:
                #print "No sidelayer"
                layer.options["init_offset"]=(0, 0)
#            print "getting features for layer={0}".format(layer)
            if layer.id!=me.options['bounds']['data']['sidelayer']:
                layer.get_features()
                #print "done getting features for layer={0}".format(layer.id)

       
        # initialize the projected bounds of the main layer and the sidelayer
        me._auto_scale_factor=me._init_projected_bounds()
        
        # scale and offset the side features
        me._scale_offset_side_features()
        # Do the view AFTER projecting 
        me.view = me._get_view()
        # Get the polygon (in fact it's a rectangle in most cases) that will be used
        # to clip away unneeded geometry unless *cfg['export']['crop-to-view']* is set to false.
        me.view_poly = me._init_view_poly()

        me._project_to_view()

       
        # In each layer we will join polygons.
        me._join_features()
        # Eventually we crop geometries to the map bounding rectangle.
        if options['export']['crop-to-view']:
            me._crop_layers_to_view()
        # Here's where we apply the simplification to geometries.
        me._simplify_layers()
        # Also we can crop layers to another layer, useful if we need to limit geological
        # geometries such as tree coverage to a political boundary of a country.
        me._crop_layers()
        # Or subtract one layer from another (or more), for instance to cut out lakes
        # from political boundaries.
        me._subtract_layers()

    # Initialize the projected bounds after getting layers and projecting 
    def _init_projected_bounds(self):
        opts=self.options
        sidelayer_bbox=BBox()
        layer_bbox=BBox()
        data = opts['bounds']['data']
        auto_scale_factor=1
        # Set up initial projected bounding box
        if data['sidelayer'] in self.layersById:
            layer=self.layersById[data['sidelayer']]
            for feature in layer.features:
                sidelayer_bbox.join(geom_to_bbox(feature.geometry, data["min-area"]))
        self._side_projected_bounds=geom_to_bbox(self._side_bounding_geometry)
        #sidelayer_bbox
        if data['layer'] in self.layersById:
            layer=self.layersById[data['layer']]
            for feature in layer.features:
                layer_bbox.join(geom_to_bbox(feature.geometry, data["min-area"]))
            self._projected_bounds=layer_bbox
        if opts['bounds']['scale-sidelayer']=='auto':
            auto_scale_factor=opts['bounds']['scale-sidelayer-factor']*sqrt(layer_bbox.width/sidelayer_bbox.width*layer_bbox.height/sidelayer_bbox.height)
        return auto_scale_factor
            
    # Project features to view coordinates
    def _project_to_view(self):
        for layer in self.layers:
            for feature in layer.features:
                feature.project_view(self.view)
    
    # Scale and offset the side features
    def _scale_offset_side_features(self):
        opts=self.options
        self._side_offset=self._get_side_offset()
        mybbox=BBox()
        statebbox=BBox()
        data = opts['bounds']['data']
        if data['layer'] in self.layersById:
            layer=self.layersById[data['layer']]
            for feature in layer.features:
                statebbox.join(geom_to_bbox(feature.geometry, data["min-area"]))
        # Scale the features in the sidelayer if it's countylayer
        for layer in self.layers:
            if "sidelayer" in layer.options and layer.options['sidelayer']=='countylayer':
                for feature in layer.features:
                    if isinstance(feature,MultiPolygonFeature):
                        #print 'scaling feature {0}'.format(feature.props['NAME'])
                        if opts['bounds']['scale-sidelayer']=='auto':
                            feature.scale_feature(scale_factor=self._auto_scale_factor,offset=self._side_offset)
                        else:
                            feature.scale_feature(scale_factor=layer.options['scale'],offset=self._side_offset)
                    if layer.id==layer.options['sidelayer']:
                            mybbox.join(geom_to_bbox(feature.geometry, self.options['bounds']['data']["min-area"]))


        # We need to set the view up and project AFTER all this

    
        id = data['sidelayer']
        #TODO: ensure this is countylayer?
        # Check that the layer exists.
        if id not in self.layersById:
            #print 'sidelayer not found'
            raise KartographError('layer not found "%s"' % id)
            #return (0,0)
        temp_layer = self.layersById[id]
        #print 'statebbox={0}'.format(statebbox)
        self._side_projected_bounds=mybbox
        self._projected_bounds=statebbox
        #print 'Pre-first offsetting: self._side_projected_bounds={0}'.format(self._side_projected_bounds)
      

        #print 'Pre-second offsetting: self._side_projected_bounds={0}'.format(self._side_projected_bounds)

        #Choose whether to put sidelayer on side or below, depending
        
        if self._projected_bounds.width <= self._projected_bounds.height:
            # Add a little breathing room on the left
            #print 'Adding on left'
            self._next_side_offset['x'] = -self._side_projected_bounds.left+self._projected_bounds.left-self._projected_bounds.width/6.-self._side_projected_bounds.width
            self._next_side_offset['y'] = -self._side_projected_bounds.top+self._projected_bounds.top+self._projected_bounds.height/2-self._side_projected_bounds.height/2.
        else:
            # Add some breathing room on the bottom
            #print 'Adding on bottom'
            self._next_side_offset['y'] = -self._side_projected_bounds.top+self._projected_bounds.bottom+self._projected_bounds.height/6.
            self._next_side_offset['x'] = -self._side_projected_bounds.left+self._projected_bounds.left+self._projected_bounds.width/2.-self._side_projected_bounds['width']/2.
        #print 'self._next_side_offset={0}'.format(self._next_side_offset)
         # transform to offset the sidelayers
        new_proj_bbox=BBox()
        for layer in self.layers:
            if "sidelayer" in layer.options and layer.options['sidelayer']=='countylayer':
                for feature in layer.features:
                    if isinstance(feature,MultiPolygonFeature):
                        #print 'offsetting feature {0}'.format(feature.props['NAME'])
                        feature.offset_feature(self._next_side_offset)
                    if layer.id==layer.options['sidelayer']:
                        new_proj_bbox=geom_to_bbox(feature.geometry, self.options['bounds']['data']["min-area"])
        self._side_projected_bounds=new_proj_bbox
  
    def _init_projection(self):
        """
        ### Initializing the map projection
        """
        opts = self.options
        # If either *lat0* or *lon0* were set to "auto", we need to
        # compute a nice center of the projection and update the
        # projection configuration.
        autoLon = 'lon0' in opts['proj'] and opts['proj']['lon0'] == 'auto'
        autoLat = 'lat0' in opts['proj'] and opts['proj']['lat0'] == 'auto'
        if autoLon or autoLat:
            map_center = self.__get_map_center()
            #print('main map_center={0}'.format(map_center))
           
            if autoLon:
                opts['proj']['lon0'] = map_center[0]
            if autoLat:
                opts['proj']['lat0'] = map_center[1]

        # Load the projection class, if the id is known.
        if opts['proj']['id'] in projections:
            projC = projections[opts['proj']['id']]
            #print '*****projC={0}'.format(projC)
        else:
            raise KartographError('projection unknown %s' % opts['proj']['id'])
        # Populate a dictionary of projection properties that
        # will be passed to the projection constructor as keyword
        # arguments.
        p_opts = {}
        for prop in opts['proj']:
            if prop != "id":
                p_opts[prop] = opts['proj'][prop]
        #print "p_opts={0}".format(p_opts)
        return projC(**p_opts)

    def __get_map_center(self):
        """
        ### Determining the projection center
        """
        #print 'map.__get_map_center'
        # To find out where the map will be centered to we need to
        # know the geographical boundaries.
        opts = self.options
        mode = opts['bounds']['mode']
        data = opts['bounds']['data']

        #print('bound mode={0}'.format(mode))
        # If the bound mode is set to *bbox* we simply
        # take the mean latitude and longitude as center.
        if mode == 'bbox':
            lon0 = data[0] + 0.5 * (data[2] - data[0])
            lat0 = data[1] + 0.5 * (data[3] - data[1])

        # If the bound mode is set to *point* we average
        # over all latitude and longitude coordinates.
        elif mode[:5] == 'point':
            lon0 = 0
            lat0 = 0
            m = 1 / len(data)
            for (lon, lat) in data:
                lon0 += m * lon
                lat0 += m * lat

        # The computationally worst case is the bound mode
        # *polygon* since we need to load the shapefile geometry
        # to compute its center of mass. However, we need
        # to load it anyway and cache the bounding geometry,
        # so this comes at low extra cost.
        elif mode[:4] == 'poly':
            features = self._get_bounding_geometry()

            print 'len(features in bounding)={0}'.format(len(features))
            if len(features) > 0:
                if isinstance(features[0].geom, BaseGeometry):
                    #print 'MOOO'
                    (lon0, lat0) = features[0].geom.representative_point().coords[0]
                    #print 'lon0, lat0={0}'.format((lon0, lat0))
            else:
                
                lon0 = 0
                lat0 = 0
        else:
            if verbose:
                sys.stderr.write("unrecognized bound mode", mode)
        return (lon0, lat0)

    def _init_side_projection(self):
        """
        ### Initializing the map projection
        """
        opts = self.options
        # If either *lat0* or *lon0* were set to "auto", we need to
        # compute a nice center of the projection and update the
        # projection configuration.
        autoLon = 'lon0' in opts['sideproj'] and opts['sideproj']['lon0'] == 'auto'
        autoLat = 'lat0' in opts['sideproj'] and opts['sideproj']['lat0'] == 'auto'
        if autoLon or autoLat:
            map_center = self.__get_side_map_center()
            #print('side map_center={0}'.format(map_center))

            if autoLon:
                opts['sideproj']['lon0'] = map_center[0]
            if autoLat:
                opts['sideproj']['lat0'] = map_center[1]
        else:
            print 'no autolon, autolat in side_map_center'
        # Load the projection class, if the id is known.
        if opts['sideproj']['id'] in projections:
            projC = projections[opts['proj']['id']]
            #print '*****projC={0}'.format(projC)
        else:
            raise KartographError('projection unknown %s' % opts['proj']['id'])
        # Populate a dictionary of projection properties that
        # will be passed to the projection constructor as keyword
        # arguments.
        p_opts = {}
        for prop in opts['sideproj']:
            if prop != "id":
                p_opts[prop] = opts['sideproj'][prop]
        #print "p_opts={0}".format(p_opts)
        return projC(**p_opts)

    def __get_side_map_center(self):
        """
        ### Determining the projection center
        """
        print 'map.__get_map_center'
        # To find out where the map will be centered to we need to
        # know the geographical boundaries.
        opts = self.options
        mode = opts['bounds']['mode']
        data = opts['bounds']['data']

        print('bound mode={0}'.format(mode))
        # If the bound mode is set to *bbox* we simply
        # take the mean latitude and longitude as center.
        if mode == 'bbox':
            lon0 = data[0] + 0.5 * (data[2] - data[0])
            lat0 = data[1] + 0.5 * (data[3] - data[1])

        # If the bound mode is set to *point* we average
        # over all latitude and longitude coordinates.
        elif mode[:5] == 'point':
            lon0 = 0
            lat0 = 0
            m = 1 / len(data)
            for (lon, lat) in data:
                lon0 += m * lon
                lat0 += m * lat

        # The computationally worst case is the bound mode
        # *polygon* since we need to load the shapefile geometry
        # to compute its center of mass. However, we need
        # to load it anyway and cache the bounding geometry,
        # so this comes at low extra cost.
        elif mode[:4] == 'poly':
            features = self._get_side_bounding_geometry()
            print 'len(features in side bounding)={0}'.format(len(features))
            if len(features) > 0:
                if isinstance(features[0].geom, BaseGeometry):
                    print 'len(features)={0}'.format(len(features))
                    (lon0, lat0) = features[0].geom.representative_point().coords[0]
                    print 'lon0, lat0={0}'.format((lon0, lat0))
            else:
                
                lon0 = 0
                lat0 = 0
        else:
            if verbose:
                sys.stderr.write("unrecognized bound mode", mode)
        return (lon0, lat0)

    def _init_bounds(self):
        """
        ### Initialize bounding polygons and bounding box
        ### Compute the projected bounding box
        """

        opts = self.options
        proj = self.proj
        side_proj = self.side_proj
        print '  proj={0}, side_proj={1}'.format(proj, side_proj)
        mode = opts['bounds']['mode'][:]
        data = opts['bounds']['data']
        padding_dict = opts['bounds']['padding-dict']
        print 'init_bounds: mode={0}, data={1}'.format(mode, data)
        if 'padding' not in opts['bounds']:
            padding = 0. # CHANGED from 0
        else:
            padding = opts['bounds']['padding']

        # If the bound mode is set to *bbox* we simply project
        # a rectangle in lat/lon coordinates.
        if mode == "bbox":  # catch special case bbox
            sea = proj.bounding_geometry(data, projected=True)
            sbbox = geom_to_bbox(sea)
            sbbox.inflate(pad_dict=padding_dict)
            return bbox_to_polygon(sbbox)

        bbox = BBox()
        side_bbox=BBox()

        # If the bound mode is set to *points* we project all
        # points and compute the bounding box.
        if mode[:5] == "point":
            ubbox = BBox()
            for lon, lat in data:
                pt = proj.project(lon, lat)
                bbox.update(pt)
                ubbox.update((lon, lat))
            self._unprojected_bounds = ubbox

        # In bound mode *polygons*, which should correctly be
        # named gemetry, we compute the bounding boxes of every
        # geometry.
        if mode[:4] == "poly":
            features = self._get_bounding_geometry()

            ubbox = BBox()
            side_ubbox=BBox()
            if len(features) > 0:
                for feature in features:
                    ubbox.join(geom_to_bbox(feature.geometry))
                    feature.project(self.side_proj)
                    fbbox = geom_to_bbox(feature.geometry, data["min-area"])
                    bbox.join(fbbox)
  #              # Save the unprojected bounding box for later to
  #              # determine what features can be skipped.
                print 'ubbox={0}'.format(ubbox)
            else:
                raise KartographError('no features found for calculating the map bounds')


           
            side_features = self._get_side_bounding_geometry()
            if len(side_features) > 0:
                print "Found_side features, len={0}".format(len(side_features))
                self._side_bounding_geometry = deepcopy(side_features[0].geom)
                side_ubbox.join(geom_to_bbox(side_features[0].geometry))
                side2=deepcopy(side_features[0])
                side2.project(self.side_proj)
                side_fbbox=geom_to_bbox(side2.geometry, data["min-area"])
                side_bbox.join(side_fbbox)
                self._side_projected_bounds=side_bbox
                print 'self._side_projected_bounds={0}'.format(self._side_projected_bounds)
            else:
                print "Cannot find side features\n"
        # If we need some extra geometry around the map bounds, we inflate
        # the bbox according to the set *padding*.
                           
            ubbox.inflate(pad_dict=padding_dict)
            self._unprojected_bounds = ubbox


        # Instead of this, expand it by the size of the second geometry thing?
        bbox.inflate(pad_dict=padding_dict)
        # At the end we convert the bounding box to a Polygon because
        # we need it for clipping tasks.
        self._projected_bounds=bbox
        return bbox_to_polygon(bbox)

    
    def _get_bounding_geometry(self):
        """
        ### Get bounding geometry
        For bounds mode "*polygons*" this helper function
        returns a list of all geometry that the map should
        be cropped to.
        """
        print ''
    #    proj = self.proj

        # Use the cached geometry, if available.
        if self._bounding_geometry_cache:
            return self._bounding_geometry_cache

        opts = self.options
        features = []
        data = opts['bounds']['data']
        id = data['layer']
        #TODO: ensure this is statelayer?
        # Check that the layer exists.
        if id not in self.layersById:
            raise KartographError('layer not found "%s"' % id)
        layer = self.layersById[id]
        
        print 'map._get_bounding_geometry,\tid={0}'.format(id)
        # Construct the filter function of the layer, which specifies
        # what features should be excluded from the map completely.
        if layer.options['filter'] is False:
            layerFilter = lambda a: True
        else:
            layerFilter = lambda rec: filter_record(layer.options['filter'], rec)

        # Construct the filter function of the boundary, which specifies
        # what features should be excluded from the boundary calculation.
        # For instance, you often want to exclude Alaska and Hawaii from
        # the boundary computation of the map, although a part of Alaska
        # might be visible in the resulting map.
        if data['filter']:
            boundsFilter = lambda rec: filter_record(data['filter'], rec)
        else:
            boundsFilter = lambda a: True

        # Combine both filters to a single function.
        filter = lambda rec: layerFilter(rec) and boundsFilter(rec)
        # Load the features from the layer source (e.g. a shapefile).

        print 'Getting features for id={0}'.format(id)
        features.extend(layer.source.get_features(
                filter=filter,
                min_area=data["min-area"],
                charset=layer.options['charset'],
                bounding=True)
                        )

    # Omit tiny islands, if needed.
        if layer.options['filter-islands']:
            features = [f for f in features
                        if f.geometry.area > layer.options['filter-islands']]
        # Store computed boundary in cache.
#        self._bounding_geometry_cache = features
        bound_feat=join_features(features,[])
        return bound_feat#features

    def _get_side_bounding_geometry(self):
        """
        ### Get side bounding geometry 

        For bounds mode "*polygons*" this helper function
        returns the offset required as a result of scaling the side layer
        """
        print 'map._get_side_bounding_geometry'
#        proj = self.proj


        opts = self.options
        features = []
        data = opts['bounds']['data']
        id = data['sidelayer']
        #TODO: ensure this is countylayer?
        # Check that the layer exists.
        if id not in self.layersById:
            print 'sidelayer not found'
            return (0,0)
#            raise KartographError('layer not found "%s"' % id)
        layer = self.layersById[id]
        
        #print 'layer={0},\tid={1}'.format(layer,id)
        # Construct the filter function of the layer, which specifies
        # what features should be excluded from the map completely.
        if layer.options['filter'] is False:
            layerFilter = lambda a: True
        else:
            layerFilter = lambda rec: filter_record(layer.options['filter'], rec)

        # Construct the filter function of the boundary, which specifies
        # what features should be excluded from the boundary calculation.
        # For instance, you often want to exclude Alaska and Hawaii from
        # the boundary computation of the map, although a part of Alaska
        # might be visible in the resulting map.
        if data['filter']:
            boundsFilter = lambda rec: filter_record(data['filter'], rec)
        else:
            boundsFilter = lambda a: True

        # Combine both filters to a single function.
        filter = lambda rec: layerFilter(rec) and boundsFilter(rec)
        # Load the features from the layer source (e.g. a shapefile).

        features.extend(layer.source.get_features(
                filter=filter,
                min_area=data["min-area"],
                charset=layer.options['charset'],
                offset=layer.options['offset'],
                contained_geom=self._main_geom,
                init_offset=self._init_offset,
                scale=layer.options['scale'],
                bounding=True)
                        )
        print 'Done getting side bounding, features.len={0}'.format(len(features))
        temp_feat=features[0].geom
        for curr_feat in features:
            temp_feat=temp_feat.union(curr_feat.geom)
        
        bound_feat=create_feature(temp_feat,{})
        return [bound_feat]#features

    # Hopefully get a good offset to move the county table to
    

    
    def _get_side_offset(self):
  
        print 'map._get_side_offset'
#        proj = self.proj


        opts = self.options
        feature = None
        data = opts['bounds']['data']
        id = data['sidelayer']
        #TODO: ensure this is countylayer?
        # Check that the layer exists.
        if id not in self.layersById:
            print 'sidelayer not found'
            return (0,0)
#            raise KartographError('layer not found "%s"' % id)
        layer = self.layersById[id]

        ret_offset={'x':0,'y':0}

        print '&&&&layer.id={0}'.format(layer.id)
        # Get the first feature and determine the offset as a result of scaling it
        if len(layer.features)>0:
            feature=layer.features[0]
            print 'len(layer.features)={0}'.format(len(layer.features))
        else:
            print 'len(layer.features)={0}'.format(len(layer.features))
        if isinstance(feature,MultiPolygonFeature):
            print 'Computing scale factor with layer.id={0},layer.options[\'scale\']={1}'.format(layer.id,layer.options['scale'])
            if opts['bounds']['scale-sidelayer']=='auto':
                ret_offset=feature.get_scale_offset(scale_factor=self._auto_scale_factor)
            else:
                ret_offset=feature.get_scale_offset(scale_factor=layer.options['scale'])
        else:
            print 'type of feature: {0}'.format(type(feature))
                        

        return ret_offset

    def _get_view(self):
        """
        ### Initialize the view
        """
        # Compute the bounding box of the bounding polygons.
        bbox = deepcopy(self._projected_bounds)
        bbox.join(self._side_projected_bounds)      #geom_to_bbox(self.projec)
        self.src_bbox = bbox
        print 'self._projected_bounds={0}, self._side_projected_bounds={1}, bbox={2}'.format(self._projected_bounds,self._side_projected_bounds, bbox)

        exp = self.options["export"]
        w = exp["width"]
        h = exp["height"]
        padding = exp["padding"]
        ratio = exp["ratio"]

        # Compute ratio from width and height.
        if ratio == "auto":
            ratio = bbox.width / float(bbox.height)

        # Compute width or heights from ratio.
        if h == "auto":
            h = w / ratio
        elif w == "auto":
            w = h * ratio
        return View(bbox, w, h, padding)

    def _init_view_poly(self):
        """
        ### Initialize the output view polygon

        Creates a polygon that represents the rectangular view bounds
        used for cropping the geometries to not overlap the view
        """
        w = self.view.width
        h = self.view.height
        return Polygon([(0, 0), (0, h), (w, h), (w, 0)])

    def _simplify_layers(self):
        """
        ### Simplify geometries
        """
        from simplify import create_point_store, simplify_lines

        # We will use a glocal point cache for all layers. If the
        # same point appears in more than one layer, it will be
        # simplified only once.
        point_store = create_point_store()

        # Compute topology for all layers. That means that every point
        # is checked for duplicates, and eventually replaced with
        # an existing instance.
        for layer in self.layers:
            if layer.options['simplify'] is not False:
                for feature in layer.features:
                    if feature.is_simplifyable():
                        feature.compute_topology(point_store, layer.options['unify-precision'])

        # Now we break features into line segments, which makes them
        # easier to simplify.
        for layer in self.layers:
            if layer.options['simplify'] is not False:
                for feature in layer.features:
                    if feature.is_simplifyable():
                        feature.break_into_lines()

        # Finally, apply the chosen line simplification algorithm.
        total = 0
        kept = 0
        for layer in self.layers:
            if layer.options['simplify'] is not False:
                for feature in layer.features:
                    if feature.is_simplifyable():
                        lines = feature.break_into_lines()
                        lines = simplify_lines(lines, layer.options['simplify']['method'], layer.options['simplify']['tolerance'])
                        for line in lines:
                            total += len(line)
                            for pt in line:
                                if not pt.deleted:
                                    kept += 1
                        # ..and restore the geometries from the simplified line segments.
                        feature.restore_geometry(lines, layer.options['filter-islands'])
        return (total, kept)

    def _crop_layers_to_view(self):
        """
        cuts the layer features to the map view
        """
        for layer in self.layers:
            #out = []
            for feat in layer.features:
                if not feat.geometry.is_valid:
                    pass
                    #print feat.geometry
                    #_plot_geometry(feat.geometry)
                feat.crop_to(self.view_poly)
                #if not feat.is_empty():
                #    out.append(feat)
            #layer.features = out

    def _crop_layers(self):
        """
        handles crop-to
        """
        for layer in self.layers:
            if layer.options['crop-to'] is not False:
                cropped_features = []
                for tocrop in layer.features:
                    cbbox = geom_to_bbox(tocrop.geom)
                    crop_at_layer = layer.options['crop-to']
                    if crop_at_layer not in self.layersById:
                        raise KartographError('you want to substract '
                            + 'from layer "%s" which cannot be found'
                            % crop_at_layer)
                    for crop_at in self.layersById[crop_at_layer].features:
                        # Sometimes a bounding box may not exist, so get it
                        if not hasattr(crop_at.geom,'bbox'):
                            crop_at.geom.bbox = geom_to_bbox(crop_at.geom)
                        if crop_at.geom.bbox.intersects(cbbox):
                            tocrop.crop_to(crop_at.geom)
                            cropped_features.append(tocrop)
                layer.features = cropped_features

    def _subtract_layers(self):
        """
        ### Subtract geometry
        """
        # Substract geometry of a layer from the geometry
        # of one or more different layers. Added mainly
        # for excluding great lakes from country polygons.
        for layer in self.layers:
            if layer.options['subtract-from']:
                for feat in layer.features:
                    if feat.geom is None:
                        continue
                    cbbox = geom_to_bbox(feat.geom)
                    # We remove it from multiple layers, if wanted.
                    for subid in layer.options['subtract-from']:
                        if subid not in self.layersById:
                            raise KartographError('you want to subtract'
                                + ' from layer "%s" which cannot be found'
                                % subid)
                        for s in self.layersById[subid].features:
                            if s.geom and geom_to_bbox(s.geom).intersects(cbbox):
                                s.subtract_geom(feat.geom)
                # Finally, we don't want the subtracted features
                # to be included in our map.
                layer.features = []

    def _join_features(self):
        """
        ### Joins features within a layer.

        Sometimes you want to merge or join multiple features (say polygons) into
        a single feature. Kartograph uses the geometry.union() method of shapely
        to do that.
        """
        from geometry.utils import join_features

        for layer in self.layers:
            if layer.options['join'] is not False:
                unjoined = 0
                join = layer.options['join']
                # The property we want to group the features by.
                groupBy = join['group-by']
                groups = join['groups']
                if groupBy is not False and not groups:
                    # If no groups are defined, we'll create a group for each
                    # unique value of the ``group-by` property.
                    groups = {}
                    for feat in layer.features:
                        fid = feat.props[groupBy]
                        groups[fid] = [fid]

                groupFeatures = {}

                # Group all features into one group if no group-by is set
                if groupBy is False:
                    groupFeatures[layer.id] = []
                    groups = [layer.id]

                res = []
                # Find all features for each group.
                for feat in layer.features:
                    if groupBy is False:
                        groupFeatures[layer.id].append(feat)
                    else:
                        found_in_group = False
                        for g_id in groups:
                            if g_id not in groupFeatures:
                                groupFeatures[g_id] = []
                            if feat.props[groupBy] in groups[g_id] or str(feat.props[groupBy]) in groups[g_id]:
                                groupFeatures[g_id].append(feat)
                                found_in_group = True
                                break
                        if not found_in_group:
                            unjoined += 1
                            res.append(feat)

                for g_id in groups:
                    # Make a copy of the input features properties.
                    props = {}
                    for feat in groupFeatures[g_id]:
                        fprops = feat.props
                        for key in fprops:
                            if key not in props:
                                props[key] = fprops[key]
                            else:
                                if props[key] != fprops[key]:
                                    props[key] = "---"
                    # If ``group-as``was set, we store the group id as
                    # new property.
                    groupAs = join['group-as']
                    if groupAs is not False:
                        props[groupAs] = g_id

                    # group.attributes allows you to keep or define
                    # certain attributes for the joined features
                    #
                    # attributes:
                    #    FIPS_1: code    # use the value of 'code' stored in one of the grouped features
                    #    NAME:           # define values for each group-id
                    #       GO: Gorenjska
                    #       KO: Koroka
                    if 'attributes' in join:
                        attrs = join['attributes']
                        for key in attrs:
                            if key not in layer.options['attributes']:
                                # add key to layer attributes to ensure
                                # that it's being included in SVG
                                layer.options['attributes'].append({'src': key, 'tgt': key})
                            if isinstance(attrs[key], dict):
                                if g_id in attrs[key]:
                                    props[key] = attrs[key][g_id]
                            else:
                                props[key] = groupFeatures[g_id][0].props[attrs[key]]  # use first value

                    # Finally join (union) the feature geometries.
                    if g_id in groupFeatures:
                        if 'buffer' in join:
                            buffer_polygons = join['buffer']
                        else:
                            buffer_polygons = 0
                        res += join_features(groupFeatures[g_id], props, buf=buffer_polygons)

                # Export ids as JSON dict, if requested
                if join['export-ids']:
                    exp = {}
                    for g_id in groups:
                        exp[g_id] = []
                        for feat in groupFeatures[g_id]:
                            exp[g_id].append(feat.props[join['export-ids']])
                    import json
                    print json.dumps(exp)

                layer.features = res

    def compute_map_scale(me):
        """
        computes the width of the map (at the lower boundary) in projection units (typically meters)
        """
        print 'Computing map scale'
        p0 = (0, me.view.height)
        p1 = (me.view.width, p0[1])
        p0 = me.view.project_inverse(p0)
        p1 = me.view.project_inverse(p1)
        from math import sqrt
        dist = sqrt((p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2)
        return dist / me.view.width

    def scale_bar_width(me):
        print 'in scale_bar_width'
        from math import log
        scale = me.compute_map_scale()
        w = (me.view.width * 0.2) * scale
        exp = int(log(w, 10))
        nice_w = round(w, -exp)
        bar_w = nice_w / scale
        return (nice_w, bar_w)

    # Add extra CSS from special elements in layers
    def add_styling(self):
        ret_style=''
        temp_style=''
        for layer in self.layers:
            opts=layer.options
            if opts['specialstyle'] is not None and layer.special_fips is not None and len(layer.special_fips)>0:
                for curr_special in layer.special_fips:
                    temp_style=re.sub('<SPECIAL_FIPS>',curr_special,opts['specialstyle'])
                    ret_style = ret_style+ temp_style
            elif opts['specialstyle'] is not None:
                ret_style+=opts['specialstyle']
        return ret_style
        
