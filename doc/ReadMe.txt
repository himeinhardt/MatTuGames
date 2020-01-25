MATLAB TOOLBOX: MatTuGames

Contents
1. Introduction
2. Requirements
3. Installation
3.1 UNIX/Linux
3.2 Windows
3.3 Mac/OS X
4. Documentation
5. Troubleshooting


1. INTRODUCTION
#################

The Matlab toolbox MatTuGames provides about 130 functions 
for modeling and calculating solutions and properties of 
cooperative games with transferable utilities. This toolbox 
is a partial port of the Mathematica package TuGames. 
In contrast to existing Matlab toolboxes to investigate 
TU games, which are written in a C/C++ programming style 
with the consequence that these functions are executed 
relatively slowly, we heavily relied on vectorized constructs 
in order to write more efficient Matlab functions.  


2. REQUIREMENTS
#################

Although it is not required to install the toolbox written 
by Jean Derks, we recommend it to do so in order to compute 
the pre-nucleolus without the optimization toolbox. However, 
we experienced for certain game classes closed loops. 
This toolbox can be downloaded from: 

http://www.personeel.unimaas.nl/Jean-Derks/downlds/CGinMatlab20050731.zip

In order to get full operation of the consistency functions we
provide a set of adjusted and optimized files of the Derks
toolbox. This set of files have been adjusted by us to conduct a 
proper and more rapid consistency investigation. It also fixes a 
problem with closed loops under certain game classes. This set of
files can be made available upon request. Alternatively, one can 
rely on the PreNucl() function, but this requires a license of 
Matlab's optimization toolbox.

Under UNIX/Linux it is also possible to enable vertex enumerations 
of the core as well as to draw the core of a game in connection 
with some cooperative solution concepts. Getting these features operative, 
one has to install the cdd-library by Komei Fukuda, which can be 
found at the URL:

http://www.cs.mcgill.ca/~fukuda/download/cdd

There are also some binaries for Windows available on this website. 
Nevertheless, it is not sufficient to install them on your computer
to get the features above operative. In addition, one has also to 
port our shell-script "corevert" to Windows.

To have a direct access to this library, the Matlab interface for the 
CDD solver -- Cddmex -- must be installed on your system that can be found
at the URL:

http://control.ee.ethz.ch/~hybrid/cdd.php

To access the Mathematica package TuGames, the Mathematica Symbolic Toolbox
must be installed from

http://www.mathworks.com/matlabcentral/fileexchange/6044-mathematica-symbolic-toolbox-for-matlab-version-2-0

The Mathematica package TuGames Version 1.1 is available at the URL:

http://library.wolfram.com/infocenter/MathSource/5709/

The most recent version can be made available by the author of this 
toolbox upon request. 

3. INSTALLATION
#################
3.1 UNIX/Linux
-----------------
3.1.1 INSTALLING FILES
......................

Change in your $HOME directory to your "matlab" sub-directory, and unzip 
there the zip-file "mat_tugV0d4.zip". For instance,

cd matlab
unzip mat_tugV0d4.zip

or in case that you want first check out the contents of the zip file, type

unzip -v mat_tugV0d4.zip
 
on the command line.

The first operation above will create a folder named "mat_tugV0d4", 
where all the m-files and documentary files will be copied. In the next
step rename the folder name mat_tugV0d4 to mat_tug. In case of an update
move the old directory mat_tug to mat_tugV0d3.

3.1.2 SETTING ENVIRONMENT VARIABLES
...................................

Now edit the "startup.m" file or use the Matlab front-end to make the 
new directory known to your Matlab session.

3.1.3. INSTALLING AUXILIARY FILES
.................................

3.1.3.a SHELL-SCRIPT

Getting the functions "CoreVertices()" and "CorePlot()" to work, one has 
to install the files located in the sub-directories "bin", and "tools" 
in the folder "mat_tug". These are some auxiliary files that perform 
some reading/writing operations on your hard-disk, and which call the 
cdd-library. Hence, you have to install, the "lcdd" and "lcdd_gmp" binaries 
properly on your system, so that these programs can be found by the 
shell-script "corevert".

Copy the shell script in the directory "bin" to a "bin" directory that is 
known by your environment variable $PATH, that is, for example:

cp -v -i mat_tug/bin/corevert $HOME/bin/corevert 


3.1.3.b SED-File

Furthermore, create a directory named "tools" in your $HOME directory, 
and copy the sed-file in this new created directory. Hence, invoke

mkdir -v $HOME/tools

and

cp -v -i mat_tug/tools/sed_core $HOME/tools/sed_core

This file is needed to convert the game information, which are saved into 
a temporary ASCII-file, into a format that the cdd-library 
can understand. 


3.1.3.c CDD-LIBRARY

The cdd-library must be compiled by following the instructions below. 
We suppose that all compiler tools are installed on your system like a 
c/c++ compiler, binutils, make, etc.

Create first a directory, let us say, "src" somewhere in your $HOME 
directory. For doing so, invoke

mkdir src

now change to this new directory, and unpack there the source code of 
the cdd-library, hence 

cd src
tar xvzf cddlib-094f.tar.gz

This creates a sub-directory called "cddlib-094f", change in this directory by

cd cddlib-094f

and now call consecutively the following four commands or follow 
the instructions given by the cdd-library README file. 

./configure --prefix=$HOME
make 
make check
make install 

In case that one has write permission in the directory "/usr/local", 
then the "--prefix" option can be omitted. Hence, type consecutively:

./configure
make 
make check

and finally as a root type: 

sudo make install

On some systems, the following procedure is required to install 
the cdd-library. First type

su

then type in the requested root password and finish the installation with 
 
make install
 


3.1.4. FINAL COMMENTS
......................

Now, everything should be installed properly. Start a new Matlab session. 
The new Matlab toolbox should now be available.


3.2 WINDOWS
-------------

To install the MatTuGames Toolbox, unzip the zip-file mat_tugV0d4.zip, 
and place the folder containing the functions on a local hard drive or 
a network drive accessible to your computer. In the next step rename 
the folder mat_tugV0d4 to mat_tug before including the folder location 
in the MATLAB path. To set the MATLAB path, start MATLAB and then
select the File/Set Path menu item. Then select Add Folder. Use the 
navigation window to select the folder containing the functions. Click 
OK and then click Save. The functions will then be ready for use within
MATLAB.

3.3 MAC/OS X
-------------

See, the Windows section.


4. DOCUMENTATION
################## 

See the manual file "manual_mat_tugames.pdf" in the "doc" sub-directory.



5. TROUBLESHOOTING
###################

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
E-mail: Holger.Meinhardt@wiwi.uni-karlsruhe.de
