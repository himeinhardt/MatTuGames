# *MATLAB* TOOLBOX: *MatTuGames* Version 1.8.0

```
Contents
1. Introduction
2. Requirements
3. Installation
3.1 UNIX/Linux
3.2 Windows
3.3 Mac/OS X
4. Documentation
5. Troubleshooting
```

## 1. INTRODUCTION


The game theoretical *Matlab* toolbox *MatTuGames* provides about 400 functions
for modeling, and calculating some solutions as well as properties of cooperative
games with transferable utilities. In contrast to existing Matlab toolboxes to
investigate TU-games, which are written in a C/C++ programming style with the consequence
that these functions are executed relatively slowly, we heavily relied on vectorized
constructs in order to write more efficient Matlab functions. In particular, the toolbox
provides functions to compute the (pre-)kernel, (pre-)nucleolus, and anti (pre-)kernel
as well as game values like the Banzhaf, Myerson, Owen, position, Shapley, solidarity,
and coalition solidarity value and much more. In addition, we will discuss how one can
use *Matlab's Parallel Computing Toolbox* in connection with this toolbox to benefit
from a gain in performance by launching supplementary Matlab workers. Some information
are provided how to call our *Mathematica* package *TuGames* within a running Matlab session.


## 2. REQUIREMENTS

This release of *MatTuGames* was developed and tested using *Matlab
R2019b* and earlier releases. A set of functions use the *Optimization Toolbox*
and the *cdd-library* by *Komei Fukuda*, which can be found at the URL:

* [CDD](http://www.inf.ethz.ch/personal/fukudak/cdd_home/)

as well as the Matlab interface to the cdd solver *CDDMEX*: 

* [CDDMEX](http://people.ee.ethz.ch/~cohysys/cdd.html)

Alternatively, in order to get even full scope of operation of the graphical features, one can also install the *MPT3* toolbox that can be downloaded from 

* [MPT3](https://www.mpt3.org/Main/Installation)

which ships with *CDDMEX*.  We strongly recommend the user to apply
the *MPT3 toolbox*, in particular of using the graphical features of
our toolbox.

For the computation of the pre-kernel and related solutions the *SuiteSparse* for *Matlab* is 
recommend that can be got from the URL

* [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse)

If you do not want to use *SuiteSparse*, then replace the function `qr_dec` by `pinv` in all functions
for the pre-kernel and related solution. The same argument applies for the function `qrginv`. It should be noted that this may cause accuracy issues with the consequence that the result is incorrect. 

To run the toolbox even in parallel mode, *Matlab's Parallel Computing Toolbox* is needed.

For connecting the *Mathematica* Package *TuGames*, the *Mathematica Symbolic Toolbox* is required, which can be found under
the URL:

* [Mathematica Symbolic Toolbox](http://www.mathworks.com/matlabcentral/fileexchange/6044-mathematica-symbolic-toolbox-for-matlab-version-2-0)

whereas *TuGames* Version 2.5.4 can be downloaded from the URL:

* [TuGames](https://github.com/himeinhardt/TuGames)

We recommend a custom installation with paclet, which can be found at

* [Paclet](https://github.com/himeinhardt/TuGames/releases)

The *MatTuGames* toolbox should work with all platforms.

Moreover, the toolbox works also with the game theory toolbox written by *Jean Derks*, which can be requested from:

* [Derks](https://www.maastrichtuniversity.nl/jean.derks/research)

We added some adjusted files that fix a problems with closed loops under certain game classes. 

This toolbox can be used to compute the pre-nucleolus up to 10-persons, if one has no license of Matlab's optimization toolbox. 

Finally, the toolbox offers interfaces to access the solvers of CVX, CPLEX, GLPK, GUROBI, HSL, IPOPT, MOSEK, and OASES. The CPLEX interfaces are compatible with version 12.10.  

To summarize, apart of the mentioned software, the toolbox requires the following MATLAB toolboxes:

*MATLAB Parallel Server*,
*Optimization Toolbox*,
*Parallel Computing Toolbox*,
*Signal Processing Toolbox*,
*Statistics and Machine Learning Toolbox*,
*Symbolic Math Toolbox*

to get full functionality in serial as well as in parallel. 
 

## 3. INSTALLATION
### 3.0 Custom Installation

To install the toolbox, we recommend a custom installation. Having downloaded the .mltbx file, navigate to it within the Matlab file explorer, double click on the mltbx file `mat_tugV1d8.mltbx` and click "install". Alternatively, right click on the .mltbx, and click "Install."

Additional instructions can be found at the URL: 

* [mltbx](https://mathworks.com/matlabcentral/answers/242430-how-do-i-install-a-mltbx-file-from-the-filesharing-site-into-r2015a)

The mltbx file `mat_tugV1d8.mltbx` is provided at

* [mltbx-file](https://github.com/himeinhardt/MatTuGames/releases)


### 3.1 UNIX/Linux (Manual Installation)
#### 3.1.1 INSTALLING FILES


Change in your ` $HOME ` directory to your `MATLAB` sub-directory, and unzip 
there the zip-file ` mat_tugV1d8.zip `. For instance,

```
cd matlab
unzip mat_tugV1d8.zip
```

or in case that you want first check out the contents of the zip file, type

```
unzip -v mat_tugV1d8.zip
```

on the command line.

The first operation above will create a folder named ` mat_tugV1d8 `, 
where all the m-files and documentary files will be copied. In the next
step rename the folder name mat_tugV1d8 to ` mat_tug `. In case of an update
make an backup of your old directory ` mat_tug `.

#### 3.1.2 SETTING ENVIRONMENT VARIABLES


Now edit the `startup.m` file or use the Matlab front-end to make the 
new directories known to your Matlab session. For instance, insert
at the end of your startup.m file the following lines

```
addpath('~/matlab/mat_tug/mama');
addpath('~/matlab/mat_tug/mat_tugames/doc');
addpath('~/matlab/mat_tug/mat_tugames');
addpath('~/matlab/mat_tug/matlink');
addpath('~/matlab/mat_tug/pct_tugames');
addpath('~/matlab/mat_tug/mattug_aux');
```

or add the paths by selecting the appropriate menu of Matlab Command Window.


#### 3.1.3. INSTALLING AUXILIARY FILES


#### 3.1.3.a SHELL-SCRIPT

Getting the functions `CoreVertices()` and `CorePlot()` to work, one has 
to install the files located in the sub-directories `bin`, and `tools` 
in the folder `mat_tug`. These are some auxiliary files that perform 
some reading/writing operations on your hard-disk, and which call the 
cdd-library. Hence, you have to install, the ` lcdd ` and ` lcdd_gmp `
binaries properly on your system, so that these programs can be found by the 
shell-script `corevert`.

Copy the shell script in the directory `bin` to a `bin directory` that is 
known by your environment variable ` $PATH `, that is, for example:

```
cp -v -i mat_tug/bin/corevert $HOME/bin/corevert 
```

#### 3.1.3.b SED-File

Furthermore, create a directory named `tools` in your $HOME directory, 
and copy the sed-file in this new created directory. Hence, invoke

```
mkdir -v $HOME/tools
```
and

```
cp -v -i mat_tug/tools/sed_core $HOME/tools/sed_core
```

This file is needed to convert the game information, which are saved into 
a temporary ASCII-file, into a format that the cdd-library 
can understand. 


#### 3.1.3.c CDD-LIBRARY

The cdd-library must be compiled by following the instructions below. 
We suppose that all compiler tools are installed on your system like a 
c/c++ compiler, binutils, make, etc.

Create first a directory, let us say, "src" somewhere in your $HOME 
directory. For doing so, invoke

```
mkdir src
```

now change to this new directory, and unpack there the source code of 
the cdd-library, hence 

```
cd src
tar xvzf cddlib-094f.tar.gz
```

This creates a sub-directory called "cddlib-094f", change in this directory by

```
cd cddlib-094f
```

and now call consecutively the following four commands or follow 
the instructions given by the cdd-library `README` file. 

```
./configure --prefix=$HOME
make 
make check
make install 
```

In case that one has write permission in the directory ` /usr/local `, 
then the ` --prefix ` option can be omitted. Hence, type consecutively:

```
./configure
make 
make check
```

and finally as a root type: 

```
sudo make install
```

On some systems, the following procedure is required to install 
the cdd-library. First type

```
sudo su
```

then type in the requested root password and finish the installation with 

```
make install
``` 


#### 3.1.4. FINAL COMMENTS

Now, everything should be installed properly. Start a new Matlab session. 
The new Matlab toolbox should now be available.


### 3.2 WINDOWS

To install the *MatTuGames* Toolbox, unzip the zip-file mat_tugV1d8.zip, 
and place the folder containing the functions on a local hard drive or 
a network drive accessible to your computer. In the next step rename 
the folder `mat_tugV1d8` to mat_tug before including the folder location 
in the MATLAB path. To set the *MATLAB* path, start *MATLAB* and then
select the File/Set Path menu item. Then select Add Folder. Use the 
navigation window to select the folder containing the functions. Click 
OK and then click Save. The functions will then be ready for use within
MATLAB.

### 3.3 MAC/OS X


See, the Windows section.

## 4. DOCUMENTATION

See the manual file "manual_mat_tugames.pdf" in the "doc" sub-directory.


## 5. TROUBLESHOOTING

In case that you encounter some problems with the installation or that 
you notice some bugs, please don't hesitate to contact us. The author is
reachable under the e-mail address mentioned in the address field. Of course,
any comments and suggestion of improvement are highly appreciated.

Address:
Holger I. Meinhardt 
Institute of Operations Research
Karlsruhe Institute of Technology (KIT)
Englerstr. 11, Building: 11.40 
D-76128 Karlsruhe 
E-mail: Holger.Meinhardt ät wiwi.uni-karlsruhe.de
