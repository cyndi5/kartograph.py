import re
from layersource import handle_layer_source
from filter import filter_record
from geometry import BBox, create_feature
from copy import deepcopy
from osgeo.osr import SpatialReference
import pyproj
from shapely.geometry import MultiPoint, Point


_verbose = False


class MapLayer(object):

    """
    MapLayer
    --------

    Represents a layer in the map which contains a list of map features
    """

    def __init__(self, id, options, _map, cache, special_fips = None):
        # Store layer properties as instance properties
        self.id = id
        self.options = options
        self.precedence=options['precedence']
#        print 'self.options={0}'.format(self.options)
        self.map = _map
        self.cache = cache
        self.features=[]
        self.special_fips = special_fips
        self.max_area_for_circle=.001
        self.high_exp_factor=1.75
        if cache is not None and 'features' in cache:
            self.proj_feat_cache=cache['features']
        elif cache is not None:
            cache['features'] = self.proj_feat_cache={}
        if 'class' not in options:
            self.classes = []
        elif isinstance(options['class'], basestring):
            self.classes = options['class'].split(' ')
        elif isinstance(options['class'], list):
            self.classes = options['class']
        # Make sure that the layer id is unique within the map.
        while self.id in self.map.layersById:
            self.id += "_"
        # Instantiate the layer source which will generate features from the source
        # geo data such as shapefiles or virtual sources such as graticule lines.
        self.source = handle_layer_source(self.options, self.cache)
    def get_features(layer, filter=False, min_area=0, contained_geom=None):
        """
        ### get_features()
        Returns a list of projected and filtered features of a layer.
        """
        opts = layer.map.options
     #   print 'Getting features for layer.id={0}'.format(layer.id)
#        print 'layer.map.options={0}'.format(layer.map.options)
        is_projected = False # should this be left?
      #  print 'First Hash of layer.map._side_bounding_geometry={0}'.format(hash(str(layer.map._side_bounding_geometry)))
        bounding_geom=None
        # Let's see if theres a better bounding box than this..
        bbox = [-180, -90, 180, 90]

#        if layer.map.proj:
#            layer.source.proj = layer.map.proj # TODO: Remove?


        # Use the clipping mode defined in the map configuration
        if opts['bounds']['mode'] == "bbox":
            bbox = opts['bounds']['data']
        # The 'crop' property overrides the clipping settings
        if 'crop' in opts['bounds'] and opts['bounds']['crop']:
            # If crop is set to "auto", which is the default behaviour, Kartograph
            # will use the actual bounding geometry to compute the bounding box
            if opts['bounds']['crop'] == "auto":
                if layer.map._unprojected_bounds:
                    bbox = layer.map._unprojected_bounds
                    bbox.inflate(inflate=1,pad_dict = opts['bounds']['padding-dict'])
                elif _verbose:
                    pass
                    #print 'could not compute bounding box for auto-cropping'
            else:
                # otherwise it will use the user defined bbox in the format
                # [minLon, minLat, maxLon, maxLat]
                bbox = opts['bounds']['crop']
            if "sidelayer" in layer.options and opts['bounds']['data']['sidelayer']!=layer.id and layer.map._side_bounding_geometry is not None:
                # We are cropping/removing stuff based on whether it intersects the sidelayer
                bounding_geom=layer.map._side_bounding_geometry
                #print('\tSetting bounding_geom to side_bounding_geometry, hash of which is {0}'.format(hash(str(bounding_geom))))
        # If the layer has the "src" property, it is a **regular map layer** source, which
        # means that there's an exernal file that we load the geometry and meta data from.
        if 'src' in layer.options:
            if layer.options['filter'] is False:
                filter = None
            else:
                filter = lambda rec: filter_record(layer.options['filter'], rec)

            # Now we ask the layer source to generate the features that will be displayed
            # in the map.
#            print 'layer.options["init_offset"]={0}'.format(layer.options['init_offset'])
            features = layer.source.get_features(
                filter=filter,
                bbox=bbox,
                ignore_holes='ignore-holes' in layer.options and layer.options['ignore-holes'],
                charset=layer.options['charset'], 
                bounding_geom=bounding_geom,
                contained_geom=contained_geom
            )
            if _verbose:
                #print 'loaded %d features from shapefile %s' % (len(features), layer.options['src'])
                pass

        # In contrast to regular layers, the geometry for **special (or virtual) layers** is generated
        # by Kartograph itself, based on some properties defined in the layer config.
        elif 'special' in layer.options:
            # The graticule layer generates line features for longitudes and latitudes
            if layer.options['special'] == "graticule":
                lats = layer.options['latitudes']
                lons = layer.options['longitudes']
                features = layer.source.get_features(lats, lons, layer.map.proj, bbox=bbox)

            # The "sea" layer generates a MultiPolygon that represents the entire boundary
            # of the map. Especially useful for non-cylindrical map projections.
            elif layer.options['special'] == "sea":
                features = layer.source.get_features(layer.map.proj)
                is_projected = True
        
        # If we're in the sidelayer main (e.g. countylayer for our application),
        # note that EVERY county should be specially styled

        layer.special_fips=[]
        for feature in features:
            #print 'feature={0}'.format(feature)
            if 'sidelayer' in layer.options:
                #print 'projecting with side_proj'
                feature.project(layer.map.side_proj)
            elif 'COUNTYFP' in feature.props and feature.props['COUNTYFP'] in layer.proj_feat_cache:
                feature.geometry=deepcopy(layer.proj_feat_cache[feature.props['COUNTYFP']])
                #print 'Found cached feature {0}'.format(feature.props['NAME'])
            else:
                feature.project(layer.map.proj)
                layer.proj_feat_cache[feature.props['COUNTYFP']]=deepcopy(feature.geometry)
                #print 'Caching feature {0}'.format(feature.props['NAME'])
            curr_fp=''
            if 'sidelayer' in layer.options and layer.options['sidelayer'] == layer.id:

                for curr_prop in feature.props:
                    temp_fp=(re.search('FP',curr_prop))
                    if temp_fp is not None and len(temp_fp.group(0))>len(curr_fp):
                        layer.special_fips.append(feature.props[curr_prop])
                        curr_fp=temp_fp.group(0)
            
  

            #It's after this point that we want to adjust the features with scaling and such

        # Remove features that don't intersect our clipping polygon
#        if layer.map.view_poly:
#            features = [feature for feature in features
#            if feature.geometry and feature.geometry.intersects(layer.map.view_poly)]
        layer.features = features
        #print 'len of features={0}'.format(len(layer.features))

    # Finds the first feature in the layer whose props match, return None if none found 
    def find_feature(self, props):
        for feat in self.features:
            if feat.props==props:
                return feat
        return None

    # Add a highlighter to this layer, containing areas is a list of
    # features, feat is a single feature that we want to highlight
    def add_highlight(self, containing_areas, feat):
        new_feat_geom = self.get_high_circle(containing_areas.features, feat.geometry)
        if new_feat_geom is not None:
            curr_props={'STATEFP':'00', 'COUNTYFP': '000', 'NAME': 'Highlight of '+feat.props['NAME']}
            high_feat=create_feature(new_feat_geom, curr_props)
            self.features=[high_feat]
        else:
            self.features=[]

        
    def get_high_circle(self, res, contained_geom):
        containing_areas=0
        max_dist=0
        for geom in [feat.geom for feat in res]:
            containing_areas+=geom.area
        if contained_geom.area/containing_areas < self.max_area_for_circle:
            # Find a good bound for the circle
            my_hull=contained_geom.convex_hull
            centroid=contained_geom.centroid
            for curr_pt in [Point(x) for x in my_hull.exterior.coords]:
                temp_dist=curr_pt.distance(centroid)
                if temp_dist>max_dist:
                    max_dist=temp_dist
            return centroid.buffer(max_dist*self.high_exp_factor)
        else:
            return None
        return None
