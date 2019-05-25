# Definitions
folders= slurmlogs ../input ../output

shp_file_ext := shp shx dbf prj

inputs=\
	$(foreach ext,$(shp_file_ext),../input/townships_2000_epsg4326.$(ext)) \
	../input/F152000.v4b_web.stable_lights.avg_vis.tif

outputs=\
	../output/F152000.v4b_web.stable_lights.avg_vis_China.tif \
	$(foreach ext,$(shp_file_ext),../output/nightlight_contour_2000_30.$(ext)) \
	$(foreach ext,$(shp_file_ext),../output/mapping_townships_2000_NTL30.$(ext)) \
	../output/townships_2000_NTL30.csv \
	../output/nightlight_raster_townships_2000_map.png \
	../output/radar_map_townships_2000_30.png \
	../output/assignment_map_townships_2000_metros_NTL30.png

# Recipes
all: $(folders) $(inputs) $(outputs)

$(folders):
	mkdir $@

$(foreach ext,$(shp_file_ext),../input/townships_2000_epsg4326.$(ext)): | ../input/townships_2000_epsg4326.zip
	unzip ../input/townships_2000_epsg4326.zip -d ../input
../input/F152000.v4b_web.stable_lights.avg_vis.tif:
	wget "https://ngdc.noaa.gov/eog/data/web_data/v4composites/F152000.v4.tar" -O ../input/F152000.v4.tar
	tar -C ../input/ -xvf ../input/F152000.v4.tar
	gunzip ../input/F152000.v4b_web.stable_lights.avg_vis.tif.gz
	sed -i.bak 's/raster_download: TRUE/raster_download: FALSE/' params.yaml && rm params.yaml.bak

$(outputs): $(inputs) params.yaml calls.R functions.R run.sbatch | $(folders)
	if command -v sbatch > /dev/null ; then sbatch -W run.sbatch; else Rscript calls.R ; fi
