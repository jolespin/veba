#### Patches

Complete reinstalls of *VEBA* environments and databases is time consuming so we've detailed how to do specific patches **for advanced users**. If you don't feel comfortable running these commands, then just do a fresh install if you would like to update. 


##### 1. Update environment variables of newly installed environments with an existing database

* Run the install/update_environment_variables.sh script. 

```
bash update_environment_variables.sh /path/to/veba_database
```

##### 2. Update human reference genome build (GRCh38 to CHM13v2)

* Remove the old build and download the new build

	```
	# Specify your database directory
	DATABASE_DIRECTORY=/path/to/veba_database # You need to edit this line
	
	# Remove the old build
	rm -rfv ${DATABASE_DIRECTORY}/Contamination/grch38
	
	# Download the new build
	wget -v -P ${DATABASE_DIRECTORY} https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip
	
	# Decompress the archive
	unzip -d ${DATABASE_DIRECTORY}/Contamination/ ${DATABASE_DIRECTORY}/chm13v2.0.zip

	# Remove the archive
	rm -rf ${DATABASE_DIRECTORY}/chm13v2.0.zip

	```

##### 3. Update GTDB-Tk v1.x to v2.x

* First, remove the old *GTDB-Tk* database and download the new one.
	   
	   
	```
	# Specify your database directory
	DATABASE_DIRECTORY=/path/to/veba_database # You need to edit this line
	
	# Remove the old GTDB-Tk (R202 for v1.x)
	rm -rfv ${DATABASE_DIRECTORY}/Classify/GTDBTk
	
	# Download the new GTDB-Tk (R207_v2 for 2.x) archive
	wget -v -P ${DATABASE_DIRECTORY} https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz
	
	# Decompress files
	tar xvzf ${DATABASE_DIRECTORY}/gtdbtk_r207_v2_data.tar.gz -C ${DATABASE_DIRECTORY}
	
	# Rename/move directory
	mv ${DATABASE_DIRECTORY}/release207_v2 ${DATABASE_DIRECTORY}/Classify/GTDBTk
	
	# Remove archive
	rm -rf ${DATABASE_DIRECTORY}/gtdbtk_r207_v2_data.tar.gz	
	```
	
	* Second, remove and add the new Conda environment.

	```
	# Download the v1.0.2a version and decompress .tar.gz
	wget https://github.com/jolespin/veba/archive/refs/tags/v1.0.2a.tar.gz
	
	tar -xvf v1.0.2a.tar.gz && mv veba-1.0.2a veba
	
	# Remove VEBA-binning-prokaryotic_env and VEBA-classify_env
	conda env remove -n VEBA-binning-prokaryotic_env
	conda env remove -n VEBA-classify_env
	
	# Recreate environments
	conda env create -n VEBA-binning-prokaryotic_env -f veba/install/environments/VEBA-binning-prokaryotic_env.yml
	
	conda env create -n VEBA-classify_env -f veba/install/environments/VEBA-classify_env.yml
	
	```
	
	* Third, add scripts to environment bin.

	```
	# Update scripts in all environments
	CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")

	# Environemnts
	for ENV_PREFIX in ${CONDA_BASE}/envs/VEBA-*_env; do
	    # Get environment name
	    ENV_NAME=$(basename ${ENV_PREFIX})
	
	    # Copy over files to environment bin/
	    echo -e "Copying VEBA into ${ENV_NAME} environment"
	    cp -r veba/src/*.py ${CONDA_BASE}/envs/${ENV_NAME}/bin/
	    cp -r veba/src/scripts/ ${CONDA_BASE}/envs/${ENV_NAME}/bin/
	
	    # Symlink the accessory scripts to bin/
	    echo -e "Symlinking VEBA scripts into ${ENV_NAME} environment path"
	    ln -sf ${CONDA_BASE}/envs/${ENV_NAME}/bin/scripts/* ${CONDA_BASE}/envs/${ENV_NAME}/bin/
	
	    done
	```
	* Last, you need to add the environment variables to each conda environment. 

	
	```
	# Specify your database directory
	DATABASE_DIRECTORY=/path/to/veba_database # You need to edit this line
	REALPATH_DATABASE_DIRECTORY=$(realpath $DATABASE_DIRECTORY)
	
	# Get conda location
	CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")
	
	# Add environment variables to activate.d/deactivate.d
	for ENV_NAME in VEBA-binning-prokaryotic_env VEBA-classify_env; do 
    ENV_PREFIX=${CONDA_BASE}/envs/${ENV_NAME}
    
    echo $ENV_PREFIX;
    mkdir -v -p ${ENV_PREFIX}/etc/conda/activate.d/
    mkdir -v -p ${ENV_PREFIX}/etc/conda/deactivate.d/
    
    # VEBA_DATABASE
    echo "export VEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}" > ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset VEBA_DATABASE" > ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
 
    # GTDB-Tk
    echo "export GTDBTK_DATA_PATH=${REALPATH_DATABASE_DIRECTORY}/Classify/GTDBTk/" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset GTDBTK_DATA_PATH" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    
    # CheckM
    echo "export CHECKM_DATA_PATH=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckM/" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset CHECKM_DATA_PATH" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh    
    
    done 

	```
	
