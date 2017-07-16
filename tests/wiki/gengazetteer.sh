#!/bin/sh

args=`getopt s:y:n $@`
# you should not use `getopt abo: "$@"` since that would parse
# the arguments differently from what the set command below does.
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

url=http://www2.census.gov/geo/tiger/GENZ${yarg}/shp/cb_${yarg}_${sarg}_place_500k.zip
echo ${url}
wget ${url}

unzip cb_${yarg}_${sarg}_place_500k.zip -d cb_${yarg}_${sarg}_place_500k
wget https://www2.census.gov/geo/docs/maps-data/data/gazetteer/${yarg}_Gazetteer/${yarg}_gaz_place_${sarg}.txt

echo ${sarg} ${yarg}
echo python generatefromgazetteer.py -s ${sarg} -y ${yarg} $sflags
python generatefromgazetteer.py -s ${sarg} -y ${yarg} $sflags
#echo single-char flags: "'"$sflags"'"
#echo oarg is "'"$oarg"'"
