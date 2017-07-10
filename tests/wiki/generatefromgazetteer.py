import sys
from kartograph import Kartograph
import argparse
import re
import cProfile, pstats, StringIO
import io

def make_this_place(K, css, curr_state, curr_place):
    cb_var='cb_2016_us_county_500k'
    
    cb_var_place='cb_2016_'+curr_state+'_place_500k'
    cfg={
    "layers": 
    {

  
     "statelayer": {
         "src":  cb_var+"/"+cb_var+".shp",
         "filter": lambda record: record['STATEFP'] == curr_state,
         "scale": 1,
         "offset": {"x": 0, "y": 0},
      },
       "countylayer": {
       "src": cb_var+"/"+cb_var+".shp",
       "filter": lambda record: record['STATEFP']==curr_state,
       "offset": {"x": -0., "y": 0},
       "scale": 10.0,
       "sidelayer":"countylayer",
       "specialstyle": "#statelayer[COUNTYFP=<SPECIAL_FIPS>]\n{\n\tfill: #e4744f;\n}\n"  # special style styles a special feature in each layer
       },
       
       "placelayer": {
        "src": cb_var_place+"/"+cb_var_place+".shp",
        "filter": lambda record: record['DESIRED_GEOM'],
        "main-filter": lambda record: record['STATEFP']==curr_state and record['PLACEFP']==curr_place,
       "offset": { "x": 0., "y": 0},
       "scale": 10.0,
      "sidelayer": "countylayer",
      "specialstyle": '#placelayer[PLACEFP='+curr_place+']\n{\n\tfill: red;\n}\n'
      }

   }, 
   "bounds": 
   {
        "padding": 1.0,
        "padding-dict": {
            "left": 0.05,
            "right": 0.05,
            "top": 0.05,
            "bottom": 0.05
        },
        "data": {
            "layer": "statelayer",
            "sidelayer": "countylayer",
            "auto-side": True
            },
        "scale-sidelayer": "auto",
        "scale-sidelayer-factor": 1.5
   },
     "export":
     {
       "width": 500,
       "height": 500,
       "padding": 15
      }
    }
    print '** Begin generating {0}'.format(curr_place)
    #generate_place(cfg, None, css, curr_state, curr_place)
    K.generate(cfg, outfile=None,stylesheet=css,render_format='Moo', curr_place=curr_place, cache_bounds=True)
    print '** End generating {0}'.format(curr_place)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Output Place Maps to Upload to Wikipedia for gazetteer file")
    parser.add_argument("-s", "--statefips", help="the state gazetteer file FIPS code to read from",default="17")
    parser.add_argument("-y", "--yearfips", help="the year of state gazetteer file FIPS code to read from",default="2016")
    parser.add_argument("-c", "--cssstyle", help="The css style file", default="style.css")
    parser.add_argument("-n", "--nocdp", action="store_true", help="Whether or not to use CDP")
    parser.add_argument("-m", "--minplace", help="Minimum place FIPS code to add", default="00000");
    parser.add_argument("-x", "--maxplace", help="Maximum place FIPS code to add", default="99999");
    parser.add_argument("-p", "--profiler", help="Output file for profiler", default="profile.out");
    
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
               if not (args.nocdp and int(field_list[4])==57) and int(args.minplace)<=int(field_list[1][2:]) and int(args.maxplace)>=int(field_list[1][2:]):
                   make_this_place(K,css,args.statefips,field_list[1][2:])

    pr.disable()
#    s=StringIO.StringIO()
    f = io.FileIO(args.profiler,mode='w');
    sortby='cumulative'
    ps = pstats.Stats(pr, stream=f).sort_stats(sortby)
    ps.print_stats()
    #print f.getvalue()
    f.close();
    
 
        



