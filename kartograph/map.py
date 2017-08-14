from shapely.geometry import Polygon, MultiPolygon
from shapely.geometry.base import BaseGeometry
from maplayer import MapLayer
from geometry.utils import join_features
from geometry import create_feature
from geometry.utils import geom_to_bbox
from geometry.utils import bbox_to_polygon
from geometry.utils import get_offset_coords, get_offset_coords_complex, get_complex_hull, get_offset_coords_super_complex
from geometry.feature import Feature, MultiPolygonFeature
from math import sqrt, atan2
from options import parse_curr_layer


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

    def __init__(me, options, layerCache, format='svg', src_encoding=None,boundCache = None, cache_bounds=False, viewCache = {}, cache_view = False, cache_union = False, unionCache={}, verbose = False):
        me.verbose = verbose
        me.options = options
        me.format = format
  
        # We will cache the projections to the view as they will be the same whenever state and county are the same 
        me.viewCache = viewCache
        me.cache_view=cache_view
        me.boundCache = boundCache
        me.cache_bounds = cache_bounds
        me.unionCache = unionCache
        me.cache_union = cache_union
        print('len of caches: viewCache={0}, boundCache={1}, unionCache={2}'.format(len(me.viewCache), len(me.boundCache), len(me.unionCache)))
#        print 'me.cache_union={0}'.format(me.cache_union)
#        print 'map.init : me.options={0}'.format(me.options)
        # List and dictionary references to the map layers.
        me.layers = []
        me.layersById = {}
        me._position_factor = 1.5
        
        # We will cache the bounding geometry since we need it twice, eventually.
        if 'bounding_geometry' in boundCache and cache_bounds and boundCache['bounding_geometry'] is not None:
            me._bounding_geometry_cache = boundCache['bounding_geometry']
           # print 'me._bounding_geometry_cache={0}'.format(me._bounding_geometry_cache)
        else:
           me._bounding_geometry_cache = False
        me._side_bounding_geometry_cache = False
        me._unprojected_bounds = None
        me._side_bounding_geometry = None
        me._projected_bounds = None
        me._side_offset = {'x':0,'y':0}
        me._n_side_off={'x':0,'y':0}
        me._init_offset = None # the initial offset as a result of scaling
        # The **source encoding** will be used as first guess when Kartograph tries to decode
        # the meta data of shapefiles etc. We use Unicode as default source encoding.
        if not src_encoding:
            src_encoding = 'utf-8'
        me._source_encoding = src_encoding

        # Construct [MapLayer](maplayer.py) instances  and store references
        # to the layers in a list and a dictionary (I guess both for many layers)
        for layer_cfg in options['layers']:
            layer_id = layer_cfg['id']
            me.print_debug('layer_id={0}'.format(layer_id))
            layer = MapLayer(layer_id, layer_cfg, me, layerCache)
            me.layers.append(layer)
            me.layersById[layer_id] = layer

        data = me.options['bounds']['data']
        cache_str = None
        # Check if everything is cached 
        if 'sidelayer' in data:
            cache_str = me._init_with_sidelayer()
        else:
            # only handles with sidelayers for now
            raise KartographError('sidelayer not found')

        if not me.cache_view or cache_str is None or cache_str not in me.viewCache:
            # Do the view AFTER projecting 
            me.view = me._get_view()
            # Get the polygon (in fact it's a rectangle in most cases) that will be used
            # to clip away unneeded geometry unless *cfg['export']['crop-to-view']* is set to false.
        else:
            me.view = me.viewCache[cache_str]['{{VIEW}}']
            me.src_bbox = me.view.bbox
            
        me.view_poly = me._init_view_poly()

        if me.cache_view:
            me._project_to_view_cached()
        else:
            me._project_to_view()

        # In each layer we will join polygons.
 #       me._join_features()
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

    # Initialize things assuming 
    def _init_with_sidelayer(self):
        opts=self.options
        data = opts['bounds']['data']               
        # Add red circle around place if it's too small? Wait, how do we know?
        self._main_feat = None
        mainFilterLayer = None
        
        for layer in self.layers:
            #print 'layer.options={0}'.format(str(layer.options))
            if "main-filter" in layer.options and layer.options["main-filter"] is not False:
                #print 'layer.options={0}'.format(str(layer.options))
                mainFilterLayer = layer
                layerMainFilter = lambda rec: filter_record(layer.options['main-filter'], rec)
        #       First, load the geometry of the place we want to highlight 
                self._main_feat = layer.source.get_main_feat(main_filter=layerMainFilter)
                self._main_place_geom = self._main_feat.geometry
    
        # Initialize the projections that will be used in this map.
        # There will be one projection focusing on the main layer 
        # 
        # This sounds easier than
        # it is since we need to compute lots of stuff here.
        self.proj = self._init_projection()
        self.side_proj = self._init_side_projection()
        #print '**self.side_proj={0}'.format(self.side_proj)

        # cache the bounding geometry for state to avoid recomputing
        if 'bounding_geometry' not in self.boundCache and self.cache_bounds:
            self.boundCache['bounding_geometry'] =  deepcopy(self._bounding_geometry_cache)

        # Compute the bounding geometry for the map.
        self.bounds_poly = self._init_bounds()
   

        # Load all features that could be visible in each layer. The feature geometries will
        # be projected 

        # Specifically load the regular sidelayer based on whether it intersects with the main geometry
        sidelayer = self.layersById[self.options['bounds']['data']['sidelayer']]
        sidelayer.get_features(contained_geom=self._main_place_geom)
        self.print_debug("done getting {0} features for layer={1}".format(len(sidelayer.features),sidelayer.id))

        cache_str = self._get_cache_str()
        if cache_str not in self.viewCache:
            self._finish_init_sidelayer(mainFilterLayer, sidelayer)
            return None
        else:
            return cache_str

        
    
    def _finish_init_sidelayer(self, mainFilterLayer, sidelayer):
        opts=self.options
        data = opts['bounds']['data']               

        # Now load the rest of the layers 
        for layer in self.layers:
            if layer.id!=data['sidelayer']:
                layer.get_features()
                self.print_debug("done getting {0} features for layer={1}".format(len(layer.features),layer.id))
        # initialize the projected bounds of the main layer and the sidelayer
        self._auto_scale_factor=self._init_projected_bounds()
        
        # scale and offset the side features
        self._scale_offset_side_features()

        # add a new highlight layer
        main_feat = mainFilterLayer.find_feature(self._main_feat.props)
        highlighter_opts = {"id": "highlightlayer","sidelayer": "countylayer",
                                "special": None, "precedence": 0}
        parse_curr_layer(highlighter_opts)
        #parse_layer_offset(highlighter_opts)
        #parse_layer_scale(highlighter_opts)
        # Create the new layer
        highlight_layer = MapLayer(highlighter_opts['id'], highlighter_opts, self, None)
        highlight_layer.add_highlight(sidelayer, main_feat)
        self.layers.append(highlight_layer)
        
    # Initialize the projected bounds after getting layers and projecting 
    def _init_projected_bounds(self):
        opts=self.options
        proj = self.proj
        side_proj = self.side_proj
        data = opts['bounds']['data']
        sidelayer_bbox=BBox()
        layer_bbox=BBox()
        auto_scale_factor=1
        # Set up initial projected bounding box
        for layer in self.layers:
            if 'sidelayer' in layer.options: #data['sidelayer'] in self.layersById:
                #layer=self.layersById[data['sidelayer']]
                for feature in layer.features:
                    #feature.project(side_proj)
                    sidelayer_bbox.join(geom_to_bbox(feature.geometry, data["min-area"]))
        #self._side_projected_bounds=geom_to_bbox(self._side_bounding_geometry)
        #sidelayer_bbox
        for layer in self.layers:
            if layer.id == data['layer']:
                for feature in layer.features:
                    #print 'feature.props={0}'.format(feature.props)
                    #feature.project(proj)
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

    def _get_cache_str(self):
        # Note if we do one state at a time, we only need the sidelayer loaded to check for caching
        opts=self.options
        data = opts['bounds']['data']
        #layer_str='{{Main}}:'
        #layer = self.layersById[data['layer']]
        #for feature in layer.features:
        #    layer_str+=feature.props['NAME'] #will be rather long, hope its not too slow ...
        sidelayer_str='{{Side}}:'
        sidelayer = self.layersById[data['sidelayer']]
        for feature in sidelayer.features:
            sidelayer_str+=feature.props['NAME'] #will be rather long, hope its not too slow ...
        cache_str=sidelayer_str # layer_str + ";" + sidelayer_str
        return cache_str

    # Project features to view coordinates using cache
    def _project_to_view_cached(self):
        cache_str = self._get_cache_str()
        if cache_str not in self.viewCache:
            # add a viewCache for this bbox
            self.viewCache[cache_str]={}
            self.viewCache[cache_str]['{{VIEW}}'] = self.view
            for layer in self.layers:
                if layer.id not in self.viewCache[cache_str]:
                    self.viewCache[cache_str][layer.id]={}
                for feature in layer.features:
                    if feature.props['NAME'] not in self.viewCache[cache_str][layer.id]:
                        feature.project_view(self.view)
                        self.viewCache[cache_str][layer.id][feature.props['NAME']]=deepcopy(feature)

                    else:
                        #print('Cached {0}'.format(feature.props['NAME']))
                        feature.geometry=deepcopy(self.viewCache[cache_str][layer.id][feature.props['NAME']].geometry)
                        #print 'cached feature {0}'.format(self.viewCache[cache_str][layer.id][feature.props['NAME']])
        else:
            # It's already been cached, load it all from cache
            for layer in self.layers:
                if len(layer.features)==0:
                    # nothing there just append everything
                    for feat_name in self.viewCache[cache_str][layer.id]:
                        layer.features.append(deepcopy(self.viewCache[cache_str][layer.id][feat_name]))
                else:
                    for feature in layer.features:
                        feature.geometry=deepcopy(self.viewCache[cache_str][layer.id][feature.props['NAME']].geometry)
    
                        
                       
            
    
    # Scale and offset the side features
    def _scale_offset_side_features(self):
        opts=self.options
        self._side_offset=self._get_side_offset() # get the amount scaling requires an offset to fix
        sidebbox=BBox()
        mainbbox=BBox()
        main_geom = None
        side_geom = None
        data = opts['bounds']['data']
        if data['layer'] in self.layersById:
            layer=self.layersById[data['layer']]
            layer_str='Main:'
            is_cached = False
            for feature in layer.features:
                layer_str+=feature.props['NAME'] #will be rather long, hope its not too slow ...
            if self.cache_union and layer_str in self.unionCache:
                is_cached = True
                main_geom = deepcopy(self.unionCache[layer_str])
            
            for feature in layer.features:
                #print 'feature.geometry={0}\nself.proj={1}'.format(feature.geometry,self.proj)
                #feature.project(self.proj)
                #print 'feature.geometry={0}'.format(feature.geometry)
                mainbbox.join(geom_to_bbox(feature.geometry,data['min-area']))
                if not is_cached and main_geom is not None:
                    main_geom=main_geom.union(feature.geometry)
                elif not is_cached:
                    main_geom = deepcopy(feature.geometry)
            if not is_cached:
                # Cache it
                self.unionCache[layer_str] = deepcopy(main_geom)
        # Scale the features in the layers associated with the sidelayer

        side_layer_str='Side:'
        for layer in [layer for layer in self.layers if "sidelayer" in layer.options]:
            for feat in [f for f in layer.features if isinstance(f,MultiPolygonFeature)]:
                side_layer_str+=feat.props['NAME']
        is_cached = False
        if self.cache_union and side_layer_str in self.unionCache:
            is_cached = True
            side_geom = deepcopy(self.unionCache[side_layer_str])

        for layer in [layer for layer in self.layers if "sidelayer" in layer.options]:
            #self.print_debug('scaling {0}'.format(layer.id))
           # for feat in layer.features:
            #    self.print_debug('(A) scaling feature {0}, {1}'.format(feat, feat.props['NAME']))
            for feat in [f for f in layer.features if isinstance(f,MultiPolygonFeature)]:
                #feat.project(self.side_proj)
                #self.print_debug('scaling feature {0}'.format(feat.props['NAME'].encode('utf-8','replace')))
                if opts['bounds']['scale-sidelayer']=='auto':
                    feat.scale_feature(scale_factor=self._auto_scale_factor,offset=self._side_offset)
                else:
                    feat.scale_feature(scale_factor=layer.options['scale'],offset=self._side_offset)

                # Everything in the sidelayer should be used to create the
                # bounding bbox 
                temp_bbox=geom_to_bbox(feat.geometry, data["min-area"])
                sidebbox.join(temp_bbox)
                if not is_cached and side_geom is not None:
                    side_geom=side_geom.union(feat.geometry)
                elif not is_cached:
                    side_geom = deepcopy(feat.geometry)
                #self.print_debug('sidebbox={0}'.format(sidebbox))
        if not is_cached:
            self.unionCache[side_layer_str] = deepcopy(side_geom)
        # Create a dummy feature to scale the side bounding geometry

        temp_feat=create_feature(self._side_bounding_geometry,{})
        temp_feat.scale_feature(scale_factor=layer.options['scale'],offset=self._side_offset)
        self._side_bounding_geometry=temp_feat.geometry
        # We need to set the view up and project AFTER all this

    
        id = data['sidelayer']
        #TODO: ensure this is countylayer?
        # Check that the layer exists.
        if id not in self.layersById:
            #print 'sidelayer not found'
            raise KartographError('layer not found "%s"' % id)
            #return (0,0)
        temp_layer = self.layersById[id]
        #print 'mainbbox={0}'.format(mainbbox)
        self._side_projected_bounds=sidebbox
        self._projected_bounds=mainbbox
#        self.print_debug('self._side_projected_bounds={0}'.format(self._side_projected_bounds))
        #print 'Pre-first offsetting: self._side_projected_bounds={0}'.format(self._side_projected_bounds)
      
        #print 'Pre-second offsetting: self._side_projected_bounds={0}'.format(self._side_projected_bounds)
        layer=self.layersById[data['sidelayer']]
        #Choose where to position and main side relatively
        self._n_side_off['x'], self._n_side_off['y'] = get_offset_coords_complex(self._projected_bounds, self._side_projected_bounds, main_geom.convex_hull, side_geom.convex_hull, self._position_factor, self)

        # temp_geom, temp_geom2 = get_offset_coords_super_complex(self._projected_bounds, self._side_projected_bounds, main_geom, side_geom, self._position_factor, self)
        # temp_geom = get_complex_hull(self._projected_bounds, self._side_projected_bounds, main_geom,  self._position_factor, self)

        temp_geom, temp_geom2, self._n_side_off['x'], self._n_side_off['y'], temp_geom3 = get_offset_coords_super_complex(self._projected_bounds, self._side_projected_bounds, main_geom, side_geom, self._position_factor,self)
        
#         temp_STATEFP='26'
#         layer=self.layersById[data['layer']]
#         if len(layer.features)>0:
#             temp_STATEFP=layer.features[0].props['STATEFP']
#             temp_feat=create_feature(temp_geom,{'NAME': 'Hull', 'LSAD': '01', 'STATEFP': temp_STATEFP, 'PLACEFP': '00000'})
#     #         #print 'Adding temp_feat={0}'.format(temp_feat)
#     # #        layer.features = [temp_feat]
#             layer.features.append(temp_feat)

#         # temp_geom2 = get_complex_hull(self._projected_bounds, self._side_projected_bounds, side_geom, self._position_factor, self)

#         layer=self.layersById[data['sidelayer']]
#         if len(layer.features)>0:
#             temp_STATEFP=layer.features[0].props['STATEFP']
#             temp_feat=create_feature(temp_geom2,{'NAME': 'Side Hull', 'LSAD': '01', 'STATEFP': temp_STATEFP, 'PLACEFP': '00000'})
#             temp_feat2=create_feature(temp_geom3,{'NAME': 'Side HullLine', 'LSAD': '01', 'STATEFP': temp_STATEFP, 'PLACEFP': '99999'})
            
#             layer.features.append(temp_feat)
#             layer.features.append(temp_feat2)
#         print 'self._n_side_off={0}'.format(self._n_side_off)
# #         transform to offset the sidelayers
        new_proj_bbox=BBox()
        for layer in self.layers:
            if "sidelayer" in layer.options and layer.options['sidelayer']==data['sidelayer']:
                for feature in layer.features:
                    if isinstance(feature,MultiPolygonFeature):
                        #print 'offsetting feature {0}'.format(feature.props['NAME'])
                        feature.offset_feature(self._n_side_off)
                    if layer.id==layer.options['sidelayer']:
                        new_proj_bbox.join(geom_to_bbox(feature.geometry,data['min-area']))  
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
            print('main map_center={0}'.format(map_center))
           
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
            temp_bbox=BBox()
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
            features = self._get_side_bounding_geometry()
            #print 'len(features in side bounding)={0}'.format(len(features))
            if len(features) > 0:
                if isinstance(features[0].geom, BaseGeometry):
                    #print 'len(features)={0}'.format(len(features))
                    (lon0, lat0) = features[0].geom.representative_point().coords[0]
                    #print 'lon0, lat0={0}'.format((lon0, lat0))
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
        #print '  proj={0}, side_proj={1}'.format(proj, side_proj)
        mode = opts['bounds']['mode'][:]
        data = opts['bounds']['data']
        padding_dict = opts['bounds']['padding-dict']
        #print 'init_bounds: mode={0}, data={1}'.format(mode, data)
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
            features = deepcopy(self._get_bounding_geometry())

            ubbox = BBox()
            side_ubbox=BBox()
            if len(features) > 0:
                for feature in features:
                    #feature.project(self.proj)
                    ubbox.join(geom_to_bbox(feature.geometry))
                    fbbox = geom_to_bbox(feature.geometry, data["min-area"])
                    bbox.join(fbbox)
  #              # Save the unprojected bounding box for later to
  #              # determine what features can be skipped.
                #print 'ubbox={0}'.format(ubbox)
            else:
                raise KartographError('no features found for calculating the map bounds')


           
            side_features = self._get_side_bounding_geometry()
            if len(side_features) > 0:
                #print "Found_side features, len={0}".format(len(side_features))
                self._side_bounding_geometry = deepcopy(side_features[0].geometry)
                side2=deepcopy(side_features[0])
                
                for feat in side_features:
                    y=0
                    #feat.project(self.side_proj)
                    self._side_bounding_geometry=self._side_bounding_geometry.union(feat.geometry)

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
        #print ''
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
        
        #print 'map._get_bounding_geometry,\tid={0}'.format(id)
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

        #print 'Getting features for id={0}'.format(id)
        features.extend(layer.source.get_features(
                filter=filter,
                min_area=data["min-area"],
                charset=layer.options['charset'],
                bounding=True)
                        )

        # Store computed boundary in cache.
#        self._bounding_geometry_cache = features
        bound_feat=join_features(features,[])
        self._bounding_geometry_cache=bound_feat
        return bound_feat#features

    def _get_side_bounding_geometry(self):
        """
        ### Get side bounding geometry 

        For bounds mode "*polygons*" this helper function
        returns the offset required as a result of scaling the side layer
        """
        if self._side_bounding_geometry_cache:
            return self._side_bounding_geometry_cache

        opts = self.options
        features = []
        data = opts['bounds']['data']
        id = data['sidelayer']
        #TODO: ensure this is countylayer?
        # Check that the layer exists.
        if id not in self.layersById:
            #print 'sidelayer not found'
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
                contained_geom=self._main_place_geom,
                bounding=True)
                        )
        #print 'Done getting side bounding, features.len={0}'.format(len(features))
        temp_feat=features[0].geom
        for curr_feat in features:
            temp_feat=temp_feat.union(curr_feat.geom)
        
        bound_feat=create_feature(temp_feat,{})
        self._side_bounding_geometry_cache=[bound_feat]
        return self._side_bounding_geometry_cache#features

    # Hopefully get a good offset to move the county table to
    

    
    def _get_side_offset(self):
  
        #print 'map._get_side_offset'
#        proj = self.proj


        opts = self.options
        feature = None
        data = opts['bounds']['data']
        id = data['sidelayer']
        #TODO: ensure this is countylayer?
        # Check that the layer exists.
        if id not in self.layersById:
            #print 'sidelayer not found'
            return (0,0)
#            raise KartographError('layer not found "%s"' % id)
        layer = self.layersById[id]

        ret_offset={'x':0,'y':0}

        #print '&&&&layer.id={0}'.format(layer.id)
        # Get the first feature and determine the offset as a result of scaling it
        if len(layer.features)>0:
            feature=layer.features[0]
            #print 'len(layer.features)={0}'.format(len(layer.features))
        #else:
        #    print 'len(layer.features)={0}'.format(len(layer.features))
        if isinstance(feature,MultiPolygonFeature):
            #print 'Computing scale factor with layer.id={0},layer.options[\'scale\']={1}'.format(layer.id,layer.options['scale'])
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
        opts=self.options
        # Compute the bounding box of the bounding polygons.
        bbox = deepcopy(self._projected_bounds)
        bbox.join(deepcopy(self._side_projected_bounds))
        self.src_bbox = bbox#geom_to_bbox(temp_bounds,opts["bounds"]["data"]["min-area"])
        #print 'bbox.width={0}, bbox.height={1}, prod={2}'.format(bbox.width,bbox.height,bbox.width*bbox.height)
        #print 'self._projected_bounds={0}, self._side_projected_bounds={1}, bbox={2}'.format(self._projected_bounds,self._side_projected_bounds, bbox)

        exp = opts["export"]
        w = exp["width"]
        h = exp["height"]
        padding = exp["padding"]
        ratio = exp["ratio"]

        # Compute ratio from width and height.
        if ratio == "auto":
            ratio = self.src_bbox.width / float(self.src_bbox.height)

        # Compute width or heights from ratio.
        if h == "auto":
            h = w / ratio
        elif w == "auto":
            w = h * ratio
        return View(self.src_bbox, w, h, padding)

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
        #print 'Computing map scale'
        p0 = (0, me.view.height)
        p1 = (me.view.width, p0[1])
        p0 = me.view.project_inverse(p0)
        p1 = me.view.project_inverse(p1)
        from math import sqrt
        dist = sqrt((p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2)
        return dist / me.view.width

    def scale_bar_width(me):
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

    def print_debug(self,msg):
        if self.verbose:
            print(msg)
        
