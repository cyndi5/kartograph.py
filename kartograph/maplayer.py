import re
from layersource import handle_layer_source
from filter import filter_record
from geometry import BBox, create_feature


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
#        print 'self.options={0}'.format(self.options)
        self.map = _map
        self.cache = cache
        self.special_fips = special_fips
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
        print 'Getting features for layer.id={0}'.format(layer.id)
#        print 'layer.map.options={0}'.format(layer.map.options)
        is_projected = False # should this be left?
        print 'First Hash of layer.map._side_bounding_geometry={0}'.format(hash(str(layer.map._side_bounding_geometry)))
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
#                    print 'computing inflate stuff'
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
                # We are cropping based on whether it's in the sidelayer
                bounding_geom=layer.map._side_bounding_geometry
                print('\tSetting bounding_geom to side_bounding_geometry, hash of which is {0}'.format(hash(str(bounding_geom))))
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
                charset=layer.options['charset'], offset=layer.options['offset'], scale=layer.options['scale'], 
                init_offset=layer.options['init_offset'],
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
        
        print 'Finished with shapelayer call'#.format(len(features))

        # Add bounding_geom to features, see what happens */
       # if bounding_geom is not None:
        #    feature = create_feature(bounding_geom, {})
         #   features.append(feature)

        # If we're in the sidelayer main (e.g. countylayer for our application),
        # note that EVERY county should be specially styled

        layer.special_fips=[]
        for feature in features:
            #print 'feature={0}'.format(feature)

            if 'sidelayer' in layer.options:
                #print 'projecting with side_proj'
                feature.project(layer.map.side_proj)
            else:
                feature.project(layer.map.proj)
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
        print 'len of features={0}'.format(len(layer.features))
