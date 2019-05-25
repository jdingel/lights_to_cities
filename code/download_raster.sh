file_name=$(wget -nv -O- 'https://ngdc.noaa.gov/eog/data/web_data/v4composites/' | grep -o F1[0-9]$1.v4.tar | uniq | tail -1)
wget "https://ngdc.noaa.gov/eog/data/web_data/v4composites/$file_name" -O ../input/$file_name
tar -C ../input -xvf ../input/$file_name
gunzip ../input/F*$1*_web.stable_lights.avg_vis.tif.gz
