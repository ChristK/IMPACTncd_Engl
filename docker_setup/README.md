# IMPACTncd prerequisite dockerfile

This repository contains the Dockerfile used to build the prerequisite container for the IMPACTncd model. The container is based on Ubuntu and includes R version 4.4.3 with package versions frozen as of 31/03/2025, using the [RStudio Package Manager](https://packagemanager.posit.co/client/#/). You can update the R packages by editing the file `r-packages.txt` and running the build script. Similarly, you can update the operating system libraries by editing the file `apt-packages.txt` and running the build script. You can check the system prerequisites for Ubuntu 24.04 (Noble)  for any additional package you may want to add `r-packages.txt` [here](https://packagemanager.posit.co/client/#/repos/cran/packages/overview?search=) and then if necessary add the requires system library to `apt-packages.txt`.To find the version of the system library, for full reproducibility, you can use the command `apt-cache policy <package_name>` from a terminal within the container.

This Docker container is required for the branch 'master' of the IMPACTncd model.

## Build and Push Using Provided Scripts

You can use the provided scripts:

You can use the provided scripts:

- **Prerequisite Container:**
  - On Linux or macOS: `./build_and_push_prerequisite.sh`
    - Use the `--push` argument to push the Docker image to Docker Hub after building.
  - On Windows: `build_and_push_prerequisite.ps1`
    - Use the `-Push` argument to push the Docker image to Docker Hub after building.

- **IMPACTncd Container:**
  - On Linux or macOS: `./build_and_push_IMPACTncd.sh`
    - Use the `--push` argument to push the Docker image to Docker Hub after building.
  - On Windows: `build_and_push_IMPACTncd.ps1`
    - Use the `-Push` argument to push the Docker image to Docker Hub after building.

