# Definitions
geo_name := $(shell grep 'geo_name:*' params.yaml | sed 's/geo_name://' | tr -d '[:space:]' | tr -d '"')
country_name := $(shell grep 'country_name:*' params.yaml | sed 's/country_name://' | tr -d '[:space:]' | tr -d '"')
brightness := $(shell grep "brightness:*" params.yaml | sed "s/brightness://" | tr -d '[:space:]' | tr -d '"')
year := $(shell grep "year:*" params.yaml | sed "s/year://" | tr -d '[:space:]' | tr -d '"')
geo_shapefile := $(shell grep "geo_shapefile:*" params.yaml | sed "s/geo_shapefile://" | tr -d '[:space:]' | tr -d '"' | sed 's/.shp//')
tar_file := $(shell wget -q -O- 'https://ngdc.noaa.gov/eog/data/web_data/v4composites/' | grep -o F[0-9][0-9]$(year).v4.tar | uniq | tail -1)
letter := $(shell if [ $(year) -lt 2010 ];then echo "b";elif [ $(year) -eq 2010 ];then echo "d";else echo "c";fi)
raster_stub := $(shell echo $(tar_file) | sed 's/.tar/$(letter)_web.stable_lights.avg_vis/')

folders = slurmlogs ../input ../output
shp_file_ext := shp shx dbf prj

inputs=\
	$(foreach ext,$(shp_file_ext),$(geo_shapefile).$(ext)) \
	../input/$(raster_stub).tif

outputs=\
	../output/$(raster_stub)_$(country_name).tif \
	$(foreach ext,$(shp_file_ext),../output/nightlight_contour_$(year)_$(brightness).$(ext)) \
	$(foreach ext,$(shp_file_ext),../output/mapping_$(geo_name)_$(year)_NTL$(brightness).$(ext)) \
	../output/$(geo_name)_$(year)_NTL$(brightness).csv \
	../output/nightlight_raster_$(geo_name)_$(year)_map.png \
	../output/radar_map_$(geo_name)_$(year)_$(brightness).png \
	../output/assignment_map_$(geo_name)_$(year)_metros_NTL$(brightness).png

# Recipes
all: $(folders) $(inputs) $(outputs)

$(folders):
	mkdir $@

$(foreach ext,$(shp_file_ext),$(geo_shapefile).$(ext)): | $(geo_shapefile).zip
	unzip $(geo_shapefile).zip -d ../input
../input/$(raster_stub).tif:
	wget "https://ngdc.noaa.gov/eog/data/web_data/v4composites/$(tar_file)" -O ../input/$(tar_file)
	tar -C ../input/ -xvf ../input/$(tar_file)
	gunzip ../input/$(raster_stub).tif.gz
	sed -i.bak 's/raster_download: TRUE/raster_download: FALSE/' params.yaml && rm params.yaml.bak

$(outputs): $(inputs) params.yaml calls.R functions.R run.sbatch | $(folders)
	if command -v sbatch > /dev/null ; then sbatch -W run.sbatch; else Rscript calls.R ; fi
