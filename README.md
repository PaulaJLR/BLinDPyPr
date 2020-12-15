# BLinDPyPr
**Perform probe-guided blind docking with FTMap and DOCK6**  

BLinDPyPr --Blind Ligand Docking through Preliminary Probing-- is a Python pipeline that associates automation and conversion scripts with well established programs such as FTMap and DOCK6 in order to introduce a novel approach to blind docking.  
FTMap docked probe clusters are converted into DOCK6 spheres for determining binding regions. Because these spheres are solely derived from FTMap probes, their locations are contained in and specific to multiple potential binding pockets, which become the regions that are simultaneously probed and chosen by the search algorithm based on the properties of each candidate ligand.  
For more information, please access our [paper]().  
## Instalation
BLinDPyPr consists of a single python script which requires some packages and programs installed and on path.
### Python packages  
**For the main BLinDPyPr docking module**
1. Pandas
2. Numpy

**For the FTMap online access module**  
BLinDPyPr can use Selenium and Chromedriver to automatically submit jobs to FTMap and download them to the working directory after they are finished. It is possible to use the program without this option and its requirements, in this case the user must manually submit the target to FTMap and provide the results to BLinDPyPr.
1. Selenium
2. Chromedriver

### Programs  
1. PyMOL
2. UCSF Chimera
3. DOCK6 (versions 6.8 or later)

## Usage
```
~~~~~~~~~~~~~~~~~~ BLinDPyPr ~~~~~~~~~~~~~~~~~~~~~~
 Perform probe-guided docking with FTMap and DOCK6
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  By Paula Jofily
  Please cite:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

optional arguments:
-h, --help           show this help message and exit

Receptor and ftmap arguments:
--runftmap           Automatically submit a PDB file to ftmap server.
                     Requires Selenium package and a webdriver
--ftconf CONF        Configuration file for ftmap submission
--fteg               Print an example configuration file and quit blindpypr
--ftsimple PDB/ID    Submit to ftmap without the need for a configuration
                     file. Specify "--ftsimple [PDBID/file.pdb]". All other
                     config options will be set to None/false.
--ftfile PDB         If automatic submission to FTMap is not required,
                     provide the FTMap PDB result (e.g. fftmap.12345.pdb)
--quit               Quit blindpypr after running ftmap (useful if
                     crossclusters need to be chosen before docking).

Docking required arguments:
--lig_file MOL2      Ligand(s) file. Format: .mol2
--dock_path PATH     DOCK6 installation path where bin and parameters
                     folders are located. Ex.: /home/username/dock6.9

Docking optional arguments:
--cross [000,001]    FTmap crosscluster numbers to include as spheres. Ex.:
                     000,001,003,006,011. Default: all.
--receptor MOL2      User input mol2 receptor file, ready for docking. If
                     not provided, receptor.mol2 will be created from
                     protein in the ftmap pdb file using chimera dockprep.
--boxmarg angstroms  Margin (angstroms) to create dock6 receptor box around
                     ftmap probes. Default: 5
--input input.in     Input file for dock6. If not provided, a default one
                     will be written.
--nochem             Do not use chemical matching.
--grid               Use this flag if grid.bmp and grid.nrg already exist in
                     the working directory (and are compatible with receptor
                     and box).
--sphgen [0,1,2]     Create classic dock6 spheres using sphgen. Options: 0
                     (use all spheres); [X] (use sphere cluster [X]). If
                     --cross is used, sphere_selector will be used to select
                     spheres around specified crossclusters.
```
