cutting the maps the size of LAtAm

```bash
echo {1992..2018} | xargs -n 1 -P 21 bash -c $’
CLASS=$1
gdalwarp -ot Byte -dstnodata 0 -co COMPRESS=DEFLATE -co ZLEVEL=9 -cutline
/data/brunner/Masterarbeit/maps/Americas/L_C.shp
-crop_to_cutline /data/brunner/Masterarbeit/maps/ESA_LC/ESALC_${CLASS}.tif
/data/brunner/Masterarbeit/maps/ESA_LC/LCLA${CLASS}.tif
’ _
```

grouping of the 11 categories for the first comparison

```bash
pkreclass -i LCLA_2005.tif -o LCLA_2005_reclass.tif -c 10 -r 1 -c 11 -r 1 -c
12 -r 1 -c 20 -r 1 -c 30 -r 1 -c 40 -r 2 -c 50 -r 3 -c 60 -r 3 -c 61 -r 3
-c 62 -r 3 -c 70 -r 3 -c 80 -r 3 -c 90 -r 3 -c 100 -r 4 -c 110 -r 4 -c 120
-r 4 -c 121 -r 4 -c 122 -r 4 -c 130 -r 5 -c 150 -r 6 -c 153 -r 6 -c 160 -r
7 -c 170 -r 7 -c 180 -r 7 -c 190 -r 8 -c 200 -r 9 -c 210 -r 10 -c 220 -r 11
```

calculating the area for the landcover map of 2005 for comparison with
the SSP data

```bash
grass78 -f -text --tmp-location -c /home/jaime/data/annika/OUTPUT_2.tif <<EOF
r.external.out directory=/home/jaime/data/annika/ format="GTiff"
option="COMPRESS=DEFLATE,ZLEVEL=9" r.cell.area output=SAC_area.tif
units=m2 --o r.external.out -r -p EOF
echo {1..11} | xargs -n 1 -P 11 bash -c $’
CLASS=$1
pksetmask -i SAC_area.tif -m LCLA_2005_reclass.tif -co COMPRESS=DEFLATE -co
ZLEVEL=9
--msknodata ${CLASS} --nodata 0 --operator \’!\’ -o
SAC_area_class_${CLASS}.tif pkstat -hist -i SAC_area_class_${CLASS}.tif >
hist_SAC_c${CLASS}.txt awk -i inplace \’$2 > 0\’ hist_SAC_c${CLASS}.txt
sed -i 1d hist_SAC_c${CLASS}.txt awk \’{ printf "%.8f\\n", $1 * $2 }\’
hist_SAC_c${CLASS}.txt > mult_c${CLASS}.txt awk \’{SUM+=$1}END{print
SUM}\’ mult_c${CLASS}.txt > sum_class_${CLASS}.txt
’ _
i## clean a bit
rm hist_SAC_c*.txt mult_c*.txt SAC_area_class_*.tif
## Put together the outputs
for CLASS in {1..11}; do echo $CLASS $( cat sum_class_${CLASS}.txt ) >>
sum_classes.txt; done
## Calculate in hectares
awk ’{ printf "%i %.1f\n", $1, $2 / 10000 }’ sum_classes.txt >
sum_classes_ha.txt
```

reclassing maps 1992-2018 with 6 LC classes

```bash
echo {1992..2018} | xargs -n 1 -P 21 bash -c $’
CLASS=$1
pkreclass -i LCLA_${CLASS}.tif -o
/data/brunner/Masterarbeit/maps/reclass/LCLA_${CLASS}r.tif -co
COMPRESS=DEFLATE -co ZLEVEL=9 --code reclass_tb_4.txt
’ _
```

resampling the maps to 1km resolution

```bash
echo {1992..2018} | xargs -n 1 -P 21 bash -c $’
CLASS=$1
gdalwarp -tr 0.008333333333333 0.008333333333333 -co COMPRESS=DEFLATE -co
ZLEVEL=9 /data/brunner/Masterarbeit/maps/reclass/LCLA_${CLASS}r.tif
LCLA_${CLASS}_1km.tif
’ _
```

cutting the maps the size of the basin map

```bash
echo {1992..2018} | xargs -n 1 -P 27 bash -c $’
CLASS=$1
pkcrop -i /data/brunner/Masterarbeit/maps/1km/LCLA_${CLASS}_1km.tif $(pkinfo
-i /data/brunner/Masterarbeit/maps/Colombia/basin_nodata.tif -bb) -o
/data/brunner/Masterarbeit/maps/Colombia/LCCO_${CLASS}_crop.tif -ot Byte
-nodata 0 -co COMPRESS=LZW -co ZLEVEL=9
’ _
echo {1992..2018} | xargs -n 1 -P 27 bash -c $’
CLASS=$1
iigdal_edit.py -a_nodata 0
/data/brunner/Masterarbeit/maps/Colombia/LCCO_${CLASS}_crop.tif
’ _
echo {1992..2018} | xargs -n 1 -P 27 bash -c $’
CLASS=$1
pksetmask -i
/data/brunner/Masterarbeit/maps/Colombia/LCCO_${CLASS}_crop.tif -m
/data/brunner/Masterarbeit/maps/Colombia/basin_nodata.tif --msknodata 0
--nodata 0 --operator ’<’ -o
/data/brunner/Masterarbeit/maps/Colombia/LCCO_${CLASS}.tif -co
COMPRESS=DEFLATE -co ZLEVEL=9
’ _
```

cutting the maps the size of Colombia

```bash
echo {1992..2018} | xargs -n 1 -P 27 bash -c $’
CLASS=$1
gdalwarp
-ot Byte -dstnodata 0 -co COMPRESS=DEFLATE -co ZLEVEL=9 -cutline
/data/brunner/Masterarbeit/maps/shp_CO/gadm36_COL_0.shp -crop_to_cutline
/data/brunner/Masterarbeit/maps/Colombia/LCCO_${CLASS}.tif
/data/brunner/Masterarbeit/maps/Colombia/country/CO_${CLASS}.tif
’ _
```

calculate area map

```bash
grass78 -f -text --tmp-location -c
/data/brunner/Masterarbeit/maps/ESA_LC/LCLA_2005.tif <<EOF
r.external.out directory=/data/brunner/Masterarbeit/maps/area format="GTiff"
option="COMPRESS=DEFLATE,ZLEVEL=9"
r.cell.area output=SAC_area.tif units=m2 --o
r.external.out -r -p
EOF
```

calculate area for each land cover class

```bash
Echo {1..6} | xargs -n 1 -P 6 bash -c $’
CLASS=$1
pksetmask -i SAC_area_1km.tif -m
/data/brunner/Masterarbeit/maps/1km/LCLA_1996_1km.tif -co
COMPRESS=DEFLATE -co ZLEVEL=9 --msknodata ${CLASS} --nodata 0 --operator
\’!\’ -o SAC_area_class_${CLASS}.tif
iiipkstat -hist -i SAC_area_class_${CLASS}.tif > hist_SAC_c${CLASS}.txt
awk -i inplace \’$2 > 0\’ hist_SAC_c${CLASS}.txt
sed -i 1d hist_SAC_c${CLASS}.txt
awk \’{ printf "%.8f\\n", $1 * $2 }\’ hist_SAC_c${CLASS}.txt >
mult_c${CLASS}.txt
awk \’{SUM+=$1}END{print SUM}\’ mult_c${CLASS}.txt >
sum_class_${CLASS}.txt
’ _
```

