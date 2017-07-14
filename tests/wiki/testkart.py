import sys
from kartograph import Kartograph
import argparse
import cProfile, pstats, StringIO

# Testing some edge cases


    #if auto-side is set to True, we will hopefully be able to find the county layer automatically (this will make it much easier to automate creation of these things)

def my_filter_place(record):    
    return record['DESIRED_GEOM']
def my_filter_main_place(record):
    return record['STATEFP']==curr_state and record['PLACEFP']==curr_place

def generate_place(cfg, the_file, css, curr_state, curr_place):
    
# configurations for Kartography

    K = Kartograph()
    K.generate(cfg, outfile=the_file,stylesheet=css,render_format='Moo', curr_place=curr_place)

if __name__ == "__main__":
    pr = cProfile.Profile()
    pr.enable()
    parser = argparse.ArgumentParser(description="Output Place Maps to Upload to Wikipedia")
    parser.add_argument("-f", "--filename", help="the filename to output to",default="None")
    parser.add_argument("-s", "--statefips", help="the 2-digit FIPS code for the state the place is located in", default="17")
    parser.add_argument("-p", "--placefips", help="the 5-digit FIPS code for the place being drawn", default="00000")
    parser.add_argument("-c", "--cssstyle", help="The css style file", default="style.css")
    
    args=parser.parse_args()
    css=open(args.cssstyle).read()
    
    curr_state=args.statefips
    #curr_county='131'
    curr_place=args.placefips
    if args.filename=='None':
        the_file=None
    else:
        the_file=args.filename

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
            "auto-side": True,
            },
        "scale-sidelayer": "auto",
        "scale-sidelayer-factor": 1.9

   },
   "export":
   {
       "width": 500,
       "height": 500,
       "padding": 20
    },
    
    }
    generate_place(cfg, the_file, css, curr_state, curr_place)
    pr.disable()
    s=StringIO.StringIO()
    sortby='cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    #ps.print_stats()
    #print s.getvalue()
