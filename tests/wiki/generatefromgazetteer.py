import sys
from kartograph import Kartograph
import argparse
import re
import cProfile, pstats, StringIO
import io

def make_this_place(K, css, curr_year, curr_state, curr_place, cache_view, show_sub):
    cb_var='cb_'+curr_year+'_us_county_500k'
    
    cb_var_place='cb_'+curr_year+'_'+curr_state+'_place_500k'
    cb_var_cousub='cb_'+curr_year+'_'+curr_state+'_cousub_500k'
    cfg={
    "layers": 
    {

  
     "statelayer": {
         "src":  cb_var+"/"+cb_var+".shp",
         "precedence": 4,
         "filter": lambda record: record['STATEFP'] == curr_state,
         "scale": 1,
         "offset": {"x": 0, "y": 0},
      },
       "countylayer": {
           "precedence": 3,
       "src": cb_var+"/"+cb_var+".shp",
       "filter": lambda record: record['STATEFP']==curr_state,
       "sidelayer":"countylayer",
       "specialstyle": "#statelayer[COUNTYFP=<SPECIAL_FIPS>]\n{\n\tfill: #e4744f;\n}\n"+
       "#countylayer[NAME=HighlightThePlace]\n{\n\tstroke-width: 1px;\n\t"+
       "stroke: red;\n}\n"# special style styles a special feature in each layer
       },
       
       "placelayer": {
           "precedence": 1,
        "src": cb_var_place+"/"+cb_var_place+".shp",
        "filter": lambda record: record['DESIRED_GEOM'],
        "main-filter": lambda record: record['STATEFP']==curr_state and record['PLACEFP']==curr_place,
      "sidelayer": "countylayer", 
      "specialstyle": '#placelayer[PLACEFP='+curr_place+']\n{\n\tfill: red;\n}\n'
      }

   },
   "proj":
   {
       "id": "ll",
       "flip": 1
    },
    "sideproj":
    {
        "id": "ll",
        "flip": 1
    },
   "bounds": 
   {
        "padding": 1.0,
        "padding-dict": {
            "left": 0.00,
            "right": 0.00,
            "top": 0.00,
            "bottom": 0.00
        },
        "data": {
            "layer": "statelayer",
            "sidelayer": "countylayer",
            "auto-side": True,
            "min-area": 0.000
            },
        "scale-sidelayer": "auto",
        "scale-sidelayer-factor": 1.25
   },
     "export":
     {
       "width": 500,
       "height": 500,
       "padding": 10
      }
    }

    if show_sub:
        print('Showing sub')
        # Add subdivision to layers
        cfg["layers"]["countysublayer"]={
            "precedence": 2,
       "src": cb_var_cousub+"/"+cb_var_cousub+".shp",
        "filter": lambda record: record['DESIRED_GEOM'],
       "sidelayer":"countylayer"
       }
    print '** Begin generating {0}'.format(curr_place)
    K.generate(cfg, outfile=None,stylesheet=css,render_format='Moo', curr_place=curr_place, cache_bounds=True, cache_view=True, cache_union = True, verbose=True)
    print '** End generating {0}'.format(curr_place)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Output Place Maps to Upload to Wikipedia for gazetteer file")
    parser.add_argument("-s", "--statefips", help="the state gazetteer file FIPS code to read from",default="17")
    parser.add_argument("-y", "--yearfips", help="the year of state gazetteer file FIPS code to read from",default="2016")
    parser.add_argument("-c", "--cssstyle", help="The css style file", default="style.css")
    parser.add_argument("-n", "--nocdp", action="store_true", help="Whether or not to use CDP")
    parser.add_argument("-m", "--minplace", help="Minimum place FIPS code to add", default="00000")
    parser.add_argument("-x", "--maxplace", help="Maximum place FIPS code to add", default="99999")
    parser.add_argument("-p", "--profiler", help="Output file for profiler", default="profile.out")
    #parser.add_argument("-v", "--viewcache", help="Cache the view", action="store_true")
    parser.add_argument("-d", "--showsub", help="Show county subdivisions in output svgs", action="store_true")
    
    args=parser.parse_args()
    css=open(args.cssstyle).read()
    pr = cProfile.Profile()
    pr.enable()
    with open(args.yearfips+'_gaz_place_'+args.statefips+'.txt', 'r') as f:
       first_flag=False
       K=Kartograph()
       for line in f:
           if not first_flag:
               first_flag=True
           else:
               field_list = re.split('\t',line)
               if not (args.nocdp and field_list[4]=='57') and int(args.minplace)<=int(field_list[1][2:]) and int(args.maxplace)>=int(field_list[1][2:]):
                   make_this_place(K,css,args.yearfips,args.statefips,field_list[1][2:], True, args.showsub)

    pr.disable()
#    s=StringIO.StringIO()
    f = io.FileIO(args.profiler,mode='w');
    sortby='cumulative'
    ps = pstats.Stats(pr, stream=f).sort_stats(sortby)
    ps.print_stats()
    #print f.getvalue()
    f.close();
    
 
        



