import sys
from kartograph import Kartograph
import argparse


parser = argparse.ArgumentParser(description="Output python )

css=open('style.css').read()

cb_var='cb_2016_us_county_500k'


curr_state='17'
curr_county='131'
curr_place=0
the_file='mymap.svg'


#geoid_list=['1711462','1731368','1765481','1781256','1723022','1735203','1738869','1764811','1782088']
geoid_list=[]

if len(sys.argv)>1:
    the_file=sys.argv[1]
if len(sys.argv)>4:
    curr_state=sys.argv[2]
    curr_county=sys.argv[3]
    curr_place=sys.argv[4]
    # Let's make this argparse worthy and we also have to figure out CSS
    # so current county works 
   # css+='#statelayer[COUNTYFP='+curr_county+']\n{\n\tfill: #e4744f;\n}\n';
    #css+='#placelayer[PLACEFP='+curr_place+']\n{\n\tfill: red;\n}\n';


    
cb_var_place='cb_2016_'+curr_state+'_place_500k'
#print('{0}'.format(css))

def my_filter_state(record):
    #if record['STATEFP']==curr_state and record['COUNTYFP']==curr_county:
    return record['STATEFP']==curr_state
def my_filter_county(record):
    #if record['STATEFP']==curr_state and record['COUNTYFP']==curr_county:
    return record['STATEFP']==curr_state#and record['COUNTYFP']==curr_county
def my_filter_place(record):
#    print '{0}: {1}'.format(record['NAME'], record['GEOID'])
    #for i in record:
    #    print '{0}'.format(i)
    #print '\n\n'

    
    return (record['STATEFP']==curr_state and record['GEOID'] in geoid_list) or record['DESIRED_GEOM']

def my_filter_main_place(record):
    return record['STATEFP']==curr_state and record['PLACEFP']==curr_place

cfg={
   "layers": 
   {

  
     "statelayer": {
         "src":  cb_var+"/"+cb_var+".shp",
         "filter": my_filter_state,
         "scale": 1,
         "offset": {"x": 0, "y": 0},
      },
       "countylayer": {
       "src": cb_var+"/"+cb_var+".shp",
       "filter": my_filter_county,
       "offset": {"x": -0., "y": 0},
       "scale": 10.0,
       "sidelayer":"countylayer",
       "specialstyle": "#statelayer[COUNTYFP=<SPECIAL_FIPS>]\n{\n\tfill: #e4744f;\n}\n"  # special style styles a special feature in each layer
       },
       
       "placelayer": {
        "src": cb_var_place+"/"+cb_var_place+".shp",
        "filter": my_filter_place,
        "main-filter": my_filter_main_place,
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
        "scale-sidelayer-factor": 1.
   },
   "export":
   {
       "width": 500,
       "height": 500,
       "padding": 10
    }
}

    #if auto-side is set to True, we will hopefully be able to find the county layer automatically (this will make it much easier to automate creation of these things)



K = Kartograph()

K.generate(cfg, outfile=the_file,stylesheet=css,render_format='Moo', curr_place=curr_place)
