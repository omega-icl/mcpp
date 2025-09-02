# Installation instructions for MC++/CRONOS/CANON/MAGNUS (last updated 2025-09-02)

## Windows: install WSL2 with default distribution (Ubuntu 24.04 LTS)
- from Powershell or Command Prompt run `wsl --install`, then follow the instructions to set up a user account (a system restart might be required)
- all subsequent commands are to be run from the **Ubuntu terminal** (use it via Windows Terminal)
- `sudo apt update && sudo apt full-upgrade`
- `explorer.exe .` lets you browse Ubuntu directories using the Windows file explorer

## Required compilers and build tools
- `sudo apt install build-essential gfortran cmake`

## Set up Git
- `git config --global user.name "<NAME>"`
- `git config --global user.email "<EMAIL>"`
- _nano_ is the default text editor; to change it to _vim_, run `git config --global core.editor "vim"`

## SSH key
- `cd`
- generate a new keypair: `ssh-keygen -t ed25519 -m PEM`, then spam Enter to keep the default options
- print the public key with `cat ~/.ssh/id_ed25519.pub`, copy and add it to the Github account settings

---
## MC++ dependencies
- `sudo apt install liblapack-dev libblas-dev libboost-all-dev` 
- get the HSL libraries [MC13](https://www.hsl.rl.ac.uk/catalogue/mc13.html), [MC21](https://www.hsl.rl.ac.uk/catalogue/mc21.html), [MC33](https://www.hsl.rl.ac.uk/catalogue/mc13.html), and install them: `tar -xzf mc13-1.0.0.tar.gz`, `cd mc13-1.0.0` to enter the source folder, then `./configure` and `sudo make install` (repeat for the other HSL libraries)
- Armadillo dependencies: `sudo apt install libopenblas-dev libarpack2-dev libsuperlu-dev`
- Armadillo: download the [latest release](https://sourceforge.net/projects/arma/files/armadillo-14.4.6.tar.xz), extract it with `tar -xJf armadillo-14.4.6.tar.xz`, `cd armadillo-14.4.6` then `./configure` and `sudo make install`
- Boost: `sudo apt update && sudo apt install libboost-all-dev`, then verify that the boost libraries have been installed correctly with `cat /usr/include/boost/version.hpp | grep "BOOST_LIB_VERSION"`

## MC++
- `cd`, then `mkdir -p Programs/bitbucket && cd Programs/bitbucket`
- clone the repository `git clone --recurse-submodules git@github.com:omega-icl/mcpp.git mcpp40`
- `cd mcpp40`, then checkout the relevant branch with `git checkout 4.0`
- `git submodule init && git submodule update`
- `cd src`, then `make install` (make sure that the correct Python path is used in _makeoptions.mk_)

## Environment variables
- `cd`, then append the following lines to .bashrc using a text editor
```
export PYTHONPATH="${PYTHONPATH}:${HOME}/Programs/bitbucket/mcpp40/src/pymc"
```
- restart the terminal for changes to take effect, or run `source .bashrc`

---
## CRONOS dependencies
- `cd`
- SuiteSparse: `sudo apt install libsuitesparse-dev`
### SUNDIALS
- `mkdir sundials && cd sundials`, download the [latest release](https://github.com/LLNL/sundials/releases/download/v7.4.0/sundials-7.4.0.tar.gz) and save it in the folder, then extract it with `tar -xzf sundials-7.4.0.tar.gz`
- create a build directory with `mkdir build`
- `cd build`, then run `cmake -DCMAKE_INSTALL_PREFIX=/opt/sundials-7.4.0 -DEXAMPLES_ENABLE_CXX=ON -DEXAMPLES_INSTALL_PATH=/opt/sundials-7.4.0/examples -DENABLE_LAPACK=ON -DENABLE_PTHREAD=ON -DENABLE_KLU=ON -DKLU_INCLUDE_DIR=/usr/include/suitesparse -DKLU_LIBRARY_DIR=/lib/x86_64-linux-gnu ../sundials-7.4.0`, and install the library with `sudo make install`
### Alternative: CVODES only
- download the [latest release](https://github.com/LLNL/sundials/releases/download/v7.4.0/cvodes-7.4.0.tar.gz), extract it with `tar -xzf cvodes-7.4.0.tar.gz`
- `mkdir build && cd build`, then run `cmake -DCMAKE_INSTALL_PREFIX=/opt/sundials-7.4.0 -DENABLE_LAPACK=ON -DENABLE_PTHREAD=ON -DENABLE_KLU=ON -DKLU_INCLUDE_DIR=/usr/include/suitesparse -DKLU_LIBRARY_DIR=/lib/x86_64-linux-gnu ../cvodes-7.4.0`, and install the library with `sudo make install`

## CRONOS 
- `cd && cd Programs/bitbucket`, clone the repository: `git clone git@github.com:omega-icl/cronos.git cronos40`
- `cd cronos40`, then checkout the relevant branch with `git checkout 4.0`
- `cd src`, then do `make all`

## Environment variables
- `cd`, then append the following lines to .bashrc using a text editor
```
export SUNDIALS_HOME="/opt/sundials-7.4.0"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/lib/x86_64-linux-gnu:${SUNDIALS_HOME}/lib"
export PYTHONPATH="${PYTHONPATH}:${HOME}/Programs/bitbucket/cronos40/src/interface"
```
- restart the terminal for changes to take effect, or run `source .bashrc`

---
## CANON dependencies: SNOPT
- obtain [SNOPT7](https://ccom.ucsd.edu/~optimizers/solvers/snopt/)
- `sudo apt install unzip`
- unzip the file with `unzip snopt7.7.zip`, `cd snopt7`, `./configure --prefix=/opt/snopt77 --with-c --with-cpp` and finally `sudo make install`

## CANON dependencies: GAMS
- download the [installer](https://www.gams.com/download/) for GAMS (v50.4.0)
- `sudo mkdir /opt/gams/ && cd /opt/gams/` 
- move the installer into this directory `sudo mv <FILEPATH> ./` (<FILEPATH> can be a Windows path, for example _/mnt/c/Users/<username>/Downloads/linux_x64_64_sfx.exe_)
- make sure the file can be executed with `sudo chmod u+x linux_x64_64_sfx.exe` and run it `sudo ./linux_x64_64_sfx.exe`

## CANON dependencies: Gurobi
- get the [installer](https://www.gurobi.com/downloads/gurobi-software/) for Gurobi Optimizer (v12.0.3) - you will need to create an account if you do not have one already
- `cd /opt/`
- move the installer `sudo mv <FILEPATH> ./`
- `sudo tar -xzf gurobi12.0.3_linux64.tar.gz`
### Fix for Gurobi under WSL2
- `cd`
- append the following lines to _.bashrc_ to assign a fixed MAC address, otherwise the license does not work after WSL2 restarts 
```
# assign a persistent MAC address for adapter bond0 and rename it eth1
mac=1a:2b:3c:4d:5e:6f
if ! ip link show | grep -q $mac; then
  sudo ip link add bond0 type bond  # only if no bond0 adapter present already
  sudo ip link set dev bond0 down
  sudo ip link set dev bond0 address $mac
  sudo ip link set dev bond0 name eth1
  sudo ip link set dev eth1 up
fi
```
### Gurobi license
- generate a new **Named-User Academic** license on the website, run `grbgetkey <LICENSE_KEY>` and click enter to save it at the default location
- `cd`, add the following line to _.bashrc_ and restart the terminal
```
export GRB_LICENSE_FILE="${HOME}/gurobi.lic"
```

## CANON
- `cd && cd Programs/bitbucket`
- clone the repository: `git clone git@github.com:omega-icl/canon.git canon40`
- `cd canon40`, checkout the relevant branch `git checkout 4.0`
- `cd src`, update the GAMS version in _makeoptions.mk_ to 50.4 using a text editor, then `make install`

## Environment variables
- `cd`
- append the following lines to _.bashrc_
```
export GUROBI_HOME="/opt/gurobi1203/linux64"
export SNOPT_HOME="/opt/snopt77"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib:${SNOPT_HOME}/lib"
export PYTHONPATH="${PYTHONPATH}:${HOME}/Programs/bitbucket/canon40/src/interface"
```
- restart the terminal for changes to take effect, or run `source .bashrc`

---
## MAGNUS
- `cd && cd Programs/bitbucket`
- clone the repository: `git clone git@github.com:omega-icl/magnus.git`
- `cd magnus`, checkout the relevant branch `git checkout 1.0`
- `cd src`, then `make install`

## Environment variables
- `cd`
- append the following lines to _.bashrc_
```
export PYTHONPATH="${PYTHONPATH}:${HOME}/Programs/bitbucket/magnus/src/interface"
```
- restart the terminal for changes to take effect, or run `source .bashrc`

