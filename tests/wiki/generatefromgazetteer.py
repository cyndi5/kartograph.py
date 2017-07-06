import sys
from kartograph import Kartograph
import argparse
from kartfips import generate_place
import re

def make_this_place(css, curr_state, curr_place):
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
        "scale-sidelayer-factor": 1.2
   },
     "export":
     {
       "width": 500,
       "height": 500,
       "padding": 10
      }
    }
    generate_place(cfg, None, css, curr_state, curr_place)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Output Place Maps to Upload to Wikipedia for gazetteer file")
    parser.add_argument("-s", "--statefips", help="the state gazetteer file FIPS code to read from",default="17")
    parser.add_argument("-y", "--yearfips", help="the year of state gazetteer file FIPS code to read from",default="2016")
    parser.add_argument("-c", "--cssstyle", help="The css style file", default="style.css")
    parser.add_argument("-n", "--nocdp", action="store_true", help="Whether or not to use CDP")
    
    args=parser.parse_args()
    css=open(args.cssstyle).read()
    with open(args.yearfips+'_gaz_place_'+args.statefips+'.txt', 'r') as f:
       first_flag=False
       for line in f:
           if not first_flag:
               first_flag=True
           else:
               field_list = re.split('\t',line)
               if not (args.nocdp and int(field_list[4])==57):
                   make_this_place(css,args.statefips,field_list[1][2:])
    
    
 
        



