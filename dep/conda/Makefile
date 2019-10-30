# Build conda environments

# Define
MINICONDA_INSTALLER = Miniconda3-4.7.12-Linux-x86_64.sh
SHELL := /bin/bash

## All ###
.PHONY: all
all: build/install_flags/env_python3_install build/install_flags/env_python2_install build/install_flags/env_pacbio_install build/install_flags/env_tools_install


### Environment: tools - Command-line tools other than snakemake and PacBio ###

# Install packages
build/install_flags/env_tools_install: build/install_flags/env_tools_deploy
	@set -e; \
	cd build; \
	. bin/activate tools; \
	source $(abspath install_tools.sh); \
	conda list > PACKAGES_TOOLS; \
	. bin/deactivate
	date > $@

# Deploy environment
build/install_flags/env_tools_deploy: build/bin/conda
	@set -e; \
	cd build; \
	bin/conda create --yes --name tools python=3.6.2
	date > $@


### Environment: pacbio - PacBio conda toolkit ###

# Install packages
build/install_flags/env_pacbio_install: build/install_flags/env_pacbio_deploy
	@set -e; \
	cd build; \
	. bin/activate pacbio; \
	source $(abspath install_pacbio.sh); \
	conda list > PACKAGES_PACBIO; \
	. bin/deactivate
	date > $@

# Deploy environment
build/install_flags/env_pacbio_deploy: build/bin/conda
	@set -e; \
	cd build; \
	bin/conda create --yes --name pacbio python=2.7.15
	date > $@


### Environment: python2 - Base environment for python2 and libraries ###

# Install packages
build/install_flags/env_python2_install: build/install_flags/env_python2_deploy
	@set -e; \
	cd build; \
	. bin/activate python2; \
	source $(abspath install_python2.sh); \
	conda list > PACKAGES_PYTHON2; \
	. bin/deactivate
	date > $@

# Deploy environment
build/install_flags/env_python2_deploy: build/bin/conda
	@set -e; \
	cd build; \
	bin/conda create --yes --name python2 python=2.7.15
	date > $@


### Environment: python3 - Base environment for python3, libraries, and snakemake ###

# Install packages
build/install_flags/env_python3_install: build/install_flags/env_python3_deploy
	@set -e; \
	cd build; \
	. bin/activate python3; \
	source $(abspath install_python3.sh); \
	conda list > PACKAGES_PYTHON3; \
	. bin/deactivate
	date > $@

# Deploy environment
build/install_flags/env_python3_deploy: build/bin/conda
	@set -e; \
	cd build; \
	bin/conda create --yes --name python3 python=3.6.2
	date > $@

### Deploy Miniconda ###

# Deploy and setup directories for deploying environments
build/bin/conda: build/${MINICONDA_INSTALLER}
	cd build; bash ${MINICONDA_INSTALLER} -f -b -p $(abspath build)
	mkdir -p build/install_flags
	touch $@

# Download
build/${MINICONDA_INSTALLER}:
	mkdir -p build
	wget https://repo.continuum.io/miniconda/${MINICONDA_INSTALLER} -O $@


### Clean ###
.PHONY: clean
clean:
	rm -rf build
