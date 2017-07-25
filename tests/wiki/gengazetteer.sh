#!/bin/sh

args=`getopt s:y:m:x:p:nd $@`
if [ $? != 0 ]
then
    echo 'Usage: ...'
    exit 2
fi
eval set -- $args
# You cannot use the set command with a backquoted getopt directly,
# since the exit code from getopt would be shadowed by those of set,
# which is zero by definition.
for i
do
    echo "i=" $i
    case "$1"
        in
        -n|-b)
            sflags="$1 ${sflags}";
            shift;;
	-d)
	    sflags="$1 ${sflags}";
	    showsub=1
	    shift;;
        -s)
	    echo 's, ' $1 $2
            sarg="$2"; shift; 
            shift;;
	-y)
	    yarg="$2"; shift;
	    shift;;
        --)
            shift; break;;
    esac
done

echo "showsub=${showsub}"

if [ -d "cb_${yarg}_us_county_500k" ]
then
    echo "cb_${yarg}_us_county_500k exists"
else
    url2=http://www2.census.gov/geo/tiger/GENZ${yarg}/shp/cb_${yarg}_us_county_500k.zip
    wget ${url2}
    unzip cb_${yarg}_us_county_500k.zip -d cb_${yarg}_us_county_500k
fi

if [ -d "cb_${yarg}_${sarg}_place_500k" ]
then
    echo "cb_${yarg}_${sarg}_place_500k exists"
else
    url=http://www2.census.gov/geo/tiger/GENZ${yarg}/shp/cb_${yarg}_${sarg}_place_500k.zip
    echo ${url}
    wget ${url}
    unzip cb_${yarg}_${sarg}_place_500k.zip -d cb_${yarg}_${sarg}_place_500k
fi

if [ -e "${yarg}_gaz_place_${sarg}.txt" ]
then
    echo "${yarg}_gas_place_${sarg}.txt exists"
else
    wget https://www2.census.gov/geo/docs/maps-data/data/gazetteer/${yarg}_Gazetteer/${yarg}_gaz_place_${sarg}.txt
fi

if [ ${showsub} != 1 ]
then
    echo "Not showing subdivisions"
elif [  -d "cb_${yarg}_${sarg}_cousub_500k" ]
then
    echo "cb_${yarg}_${sarg}_cousub_500k exists"
else
    wget http://www2.census.gov/geo/tiger/GENZ${yarg}/shp/cb_${yarg}_${sarg}_cousub_500k.zip
    unzip cb_${yarg}_${sarg}_cousub_500k.zip -d cb_${yarg}_${sarg}_cousub_500k
fi

echo ${sarg} ${yarg}
echo python generatefromgazetteer.py -s ${sarg} -y ${yarg} $sflags
python generatefromgazetteer.py -s ${sarg} -y ${yarg} $sflags
#echo single-char flags: "'"$sflags"'"
#echo oarg is "'"$oarg"'"
