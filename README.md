
# JJ-s-E45-Geant4-simulation
Geant4 simulation for J-PARC E45 Veto counter

k18geant4


K1.8 geant4 simulation tool.

**Note that this README needs to be modified.**

## Platform

This tool is developed on the platform of KEKCC, CentOS 7.9.2009.
- g++ (GCC) 8.3.0
- ROOT 6.22/08
- Geant4 11.0.2


## Anaconda setting

To use Python,
it is necessary to build the Anaconda local environment once using the `conda` command as follows.
Note that it is recommended to use `conda install` instead of `pip install` in the anaconda environment.

```sh
$ conda create -n py37 python=3.7 # py37 is an example name
$ conda activate py37
$ conda install numpy psutil pyyaml rich
```

Add the following line in .bashrc to activate your environment.

```sh
conda activate py37
```

If the prompt header of conda is annoying, add the following line in .condarc.

```yaml
changeps1: False
```


## How to install

Set environment variables.

```shell
. /opt/python-3.7/etc/profile.d/conda.sh
. /group/had/sks/software/root/6.22.08/bin/thisroot.sh
. /sw/packages/geant4/11.0.2/bin/geant4.sh
. /sw/packages/geant4/11.0.2/share/Geant4-11.0.2/geant4make/geant4make.sh
export PATH=$PATH:/group/had/sks/software/unpacker/s2s/bin
export MAKEFLAGS=-j40
conda activate py37
module load gcc/830
module load git/2260
```

then

```shell
git clone ssh://sks@www-online.kek.jp:8022/~/public_html/git/k18geant4.git
cd k18geant4
git checkout e72
./build.sh
```



## How to use

Arguments of ConfFile and OutputName are necessary.
G4Macro is an optional argument.

```shell
./bin/G4HypTPC [ConfFile] [OutputName] (G4Macro)
./bin/G4HypTPC param/conf/default.conf foo.root
./bin/G4HypTPC param/conf/default.conf foo.root bar.mac
```



## Parameters

Some parameter files that are out of the git control should be linked.

```shell
ln -s /group/had/sks/E72/software/param/BEAM/* param/BEAM/
ln -s /group/had/sks/E72/software/fieldmap .
```

