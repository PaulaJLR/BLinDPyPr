# imports
import argparse
import numpy as np
import pandas as pd
import datetime
import os
import re
import time
from random import randint


### parser
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description="""
  ~~~~~~~~~~~~~~~~~~ BLinDPyPr ~~~~~~~~~~~~~~~~~~~~~~
   Perform probe-guided docking with FTMap and DOCK6
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    By Paula Jofily
    Please cite:
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~""")
# args
## positional
#parser.add_argument('lig_file', metavar='ligand', type=str, help='Ligand(s) file. Format: .mol2')
#parser.add_argument('dock_path', metavar='dock6_home', type=str, help='DOCK6 installation path where bin and parameters folders are located. Ex.: /home/username/dock6.9')

#receptor and ftmap
ft = parser.add_argument_group('Receptor and ftmap arguments')
ft.add_argument('--runftmap', action='store_const', const=True, help='Automatically submit a PDB file to ftmap server. Requires Selenium package and a webdriver')
ft.add_argument('--ftconf', metavar='CONF', type=str, help='Configuration file for ftmap submission')
ft.add_argument('--fteg', action='store_const', const=True, help='Print an example configuration file and quit blindpypr')
ft.add_argument('--ftsimple', metavar='PDB/ID', type=str, help='Submit to ftmap without the need for a configuration file. Specify "--ftsimple [PDBID/file.pdb]". All other config options will be set to None/false.')
ft.add_argument('--ftfile', metavar='PDB', type=str, help='If automatic submission to FTMap is not required, provide the FTMap PDB result (e.g. fftmap.12345.pdb)')
ft.add_argument('--quit', action='store_const', const=True, help='Quit blindpypr after running ftmap (useful if crossclusters need to be chosen before docking).')

##docking
dock_req = parser.add_argument_group('Docking required arguments')
dock_req.add_argument('--lig_file', metavar='MOL2', type=str, help='Ligand(s) file. Format: .mol2')
dock_req.add_argument('--dock_path', metavar='PATH', type=str, help='DOCK6 installation path where bin and parameters folders are located. Ex.: /home/username/dock6.9')
dock = parser.add_argument_group('Docking optional arguments')
dock.add_argument('--cross', metavar='[000,001]', type=str, help='FTmap crosscluster numbers to include as spheres. Ex.: 000,001,003,006,011. Default: all.')
dock.add_argument('--receptor', metavar='MOL2', type=str, help='User input mol2 receptor file, ready for docking. If not provided, receptor.mol2 will be created from protein in the ftmap pdb file using chimera dockprep.')
dock.add_argument('--boxmarg', metavar='angstroms', type=str, help='Margin (angstroms) to create dock6 receptor box around ftmap probes. Default: 5')
dock.add_argument('--input', metavar='input.in', type=str, help='Input file for dock6. If not provided, a default one will be written.')
dock.add_argument('--nochem', action='store_const', const=True, help='Do not use chemical matching.')
dock.add_argument('--grid', action='store_const', const=True, help='Use this flag if grid.bmp and grid.nrg already exist in the working directory (and are compatible with receptor and box).')
dock.add_argument('--sphgen', metavar='[0,1,2]', type=str, help='Create classic dock6 spheres using sphgen. Options: 0 (use all spheres); [X] (use sphere cluster [X]). If --cross is used, sphere_selector will be used to select spheres around specified crossclusters.')


args = parser.parse_args();

### while testing
#os.chdir('/mnt/sda/paula/programs/autochemical')
#args = parser.parse_args(['fftmap.78236.pdb', 'crystal_ligand.mol2', '/home/paula/dock.6.8_source/dock6', '--sphgen', '--cross', '000,001,003,006,011'])
##
workingdir = os.getcwd()
os.chdir(workingdir)


### logfile
logfile = open('blindpypr.log', 'a+')
logfile.write("""
  ~~~~~~~~~~~~~~~~~~ BLinDPyPr ~~~~~~~~~~~~~~~~~~~~~~
   Perform probe-guided docking with FTMap and DOCK6
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    By Paula Jofily
    Please cite:
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n""")
logfile.write("  -- Section created on "+datetime.datetime.now().strftime("%d %b %Y %H:%M"+" --\n\n"))

def logit(string): #write to log and print stuff
    print(string)
    logfile.write(string+'\n')

##############################################
### AUTO ftmap                             ###
##############################################
def autoftmap(args, workingdir, logfile):
    #try:
    import selenium
    from selenium import webdriver
    if args.ftsimple != None: # if ftsimple on, se if it is a pdbid or a pdbfile and set user variables
        if args.ftsimple[-4:] == '.pdb':
            pdbpath = workingdir+'/'+args.ftsimple; pdbid = None
        else:
            pdbpath = None
            pdbid = args.ftsimple
        chains = None; mask = None; ppi_mode = False; has_nucleic = False
    else: # get user variables from config file
        try:
            infile = open(args.ftconf,'r').readlines()
        except:
            print('config file for FTMap not provided or not found.')
            quit()
        for line in infile:
            if 'pdbpath' in line:
                pdbpath = line.split()[2]
                if pdbpath == 'None': pdbpath = None
            elif 'chains' in line:
                chains = line.split()[2]
                if chains == 'None': chains = None
            elif 'pdbid' in line:
                pdbid = line.split()[2]
                if pdbid == 'None': pdbid = None
            elif 'mask' in line:
                mask = line.split()[2]
                if mask == 'None': mask = None
            elif 'ppi_mode' in line:
                ppi_mode = line.split()[2]
                if ppi_mode == 'True': ppi_mode = True
                else: ppi_mode = False
            elif 'has_nucleic' in line:
                has_nucleic = line.split()[2]
                if has_nucleic == 'True': has_nucleic = True
                else: has_nucleic = False

    ### start webdriver
    jobnm_inp = 'job'+str(randint(10000, 99999))
    options = webdriver.ChromeOptions()
    prefs = {'download.default_directory' : workingdir}
    options.add_experimental_option('prefs', prefs)
    options.add_argument('headless')
    driver = webdriver.Chrome(options=options)
    driver.get('https://ftmap.bu.edu/login.php')
    time.sleep(10)
    ### use ftmap without an account
    no_account = driver.find_element_by_partial_link_text('Use FTMap without your own account')
    no_account.click()
    time.sleep(5)

    ### job inputs
    ## job name
    jobname = driver.find_element_by_name("jobname")
    jobname.send_keys(jobnm_inp)
    ## protein
    if pdbid != None: # use a pdb code
        pdbcode = driver.find_element_by_id("protpdb")
        pdbcode.send_keys(pdbid)
        logit('      > using pdb code '+pdbid)
    else:             # upload pdb
        upload = driver.find_element_by_id("showprotfile")
        upload.click()
        choose_file = driver.find_element_by_id("prot")
        choose_file.send_keys(pdbpath)
        logit('      > uploading pdb file: '+pdbpath)
    ## specific chains?
    if chains != None:
        chain_inp = driver.find_element_by_id("protchains")
        chain_inp.send_keys(chains)
        logit('      > specified chain(s): '+chains)
    else:
        logit('      > no chains specified')
    ## advanced options
    if mask != None or ppi_mode or has_nucleic:
        advanced = driver.find_element_by_id("advancedtoggle")
        advanced.click()
        advanced_options = True
        logit('      > advanced options used')
    else:
        advanced_options = False
        logit('      > no advanced options defined')
    if advanced_options:
        if mask != None:
            mask_inp = driver.find_element_by_id("protmask")
            mask_inp.send_keys(mask)
            logit('      > using mask defined in: '+mask)
        else:
            logit('      > no protein mask defined')
        if ppi_mode:
            ppi = driver.find_element_by_id("ppimode")
            ppi.click()
            logit('      > ppi mode ON')
        else:
            logit('      > ppi mode OFF')
        if has_nucleic:
            nucl = driver.find_element_by_id("nucleic_acid")
            nucl.click()
            logit('      > has nucleic acid ON')
        else:
            logit('      > has nucleic acid OFF')

    ## run!
    mapit = driver.find_element_by_name("action")
    mapit.click()
    logit(datetime.datetime.now().strftime("%H:%M")+' > Job submitted')

    ### go to queue
    queue = driver.find_element_by_id('tabQueue')
    queue.click()
    time.sleep(5)

    ## find job ID assigned by FTMap
    page = driver.page_source.encode("utf-8")
    page = str(page)
    find = re.compile('href="jobdetail(.*)?>'+jobnm_inp); found = re.search(find, page)
    job_id = found.group().split('?')[1].split('"')[0].replace('job=', '')

    ## tell the user the details
    logit('      > job name: '+jobnm_inp)
    logit('      > job ID: '+job_id)

    ## wait until the job is done; log what is going on in the queue
    finished = False
    old_found = None
    while finished == False:
        queue = driver.find_element_by_id('tabQueue')
        queue.click()
        time.sleep(5)
        page = driver.page_source.encode("utf-8")
        page = str(page)
        try:
            find = re.compile(jobnm_inp+'(.*?</tr>)'); found = re.search(find, page)
            found = found.group().split('>')[-3].split('<')[0]
            if found != old_found:
                logit(datetime.datetime.now().strftime("%H:%M")+' > '+found)
                old_found = found
            else: continue
        except:
            finished = True
            logit(datetime.datetime.now().strftime("%H:%M")+' > job finished')

    ### go to results

    results = driver.find_element_by_id('tabResults')
    results.click()
    time.sleep(5)

    ## find and report job status
    page = driver.page_source.encode("utf-8")
    page = str(page)
    find = re.compile(jobnm_inp+'(.*?</tr>)'); found = re.search(find, page)
    #2v7a</td><td>error on local system</td></tr>
    found = found.group().split('>')[-3].split('<')[0]
    if found != 'finished':
        success = False
        logit('      > Job status: '+found+'. Will not try to download results.')
        return(success)
    else:
        success = True
        logit('      > Job status: '+found+'. Will download results to working dir.')

    try: os.system('notify-send "BLinDPyPr - FTMap Finished"')
    except: None
    ## download results if job was successful
    if success:
        job_page = driver.find_element_by_partial_link_text(job_id); job_page.click()
        time.sleep(10)

        dl_pdb = driver.find_element_by_partial_link_text('PDB')
        dl_pdb.click()
        dl_pse = driver.find_element_by_partial_link_text('PyMol session')
        dl_pse.click()
        dl_nb = driver.find_element_by_id('nblink')
        dl_nb.click()
        dl_hb = driver.find_element_by_id('hblink')
        dl_hb.click()
        dl_sum = driver.find_element_by_id('summarylink')
        dl_sum.click()
        return(job_id)
    #except:
        #return(False)


if args.fteg is True:
    print('Example configuration file for automatic FTMap submission. Do not delete lines.')
    print('If a parameter will not be not used, replace with None.')
    print('For advanced options ppi_mode and has_nucleic, use "True" or "False".')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('''pdbpath = /path/to/pdbfile/file.pdb
chains = None
pdbid = None
# advanced mode options:
mask = /path/to/maskfile/mask.pdb
ppi_mode = False
has_nucleic = False''')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    quit()

if args.runftmap is True:
    logit('Automatically running FTMap through online access:')
    autoftmap = autoftmap(args, workingdir, logfile)
    if autoftmap == False:
        print('The automatic submission to FTMap failed. Please check your internet connection, the availability of the server and if all dependencies are installed.')
        quit()
    else: ftmap_file = 'fftmap.'+autoftmap+'.pdb'
elif args.runftmap is None and args.ftfile is not None:
    ftmap_file = args.ftfile
    logit('> Using provided FTMap PDB file.')

if args.quit is True:
    logit("> Quitting BLinDPyPr after FTMap run.")
    quit()

##############################################
### ftmap probes pharmacophore preparation ###
##############################################

# functions
def separate_probes(probes, cross):
    '''separate probes'''
    # get all probes, separate
    all_probes = []
    for i in probes:
        if i[:6]=='HEADER':
            try:
                all_probes.append((probenum, probe))
            except: None
            probe = ''; probenum=i[20:23]
            probe+=(i+"\n")
        elif i == probes[-1]:
            probe+=(i+"\n"); all_probes.append((probenum, probe))
        else:
            probe+=(i+"\n")
    # filter (or not) by user chosen crossclusters
    if cross == ['all']:
        sep_probes = all_probes
    elif cross != ['all']:
        sep_probes = []
        for i in all_probes:
            if i[0] in cross:
                sep_probes.append(i)
            else: continue
    return(sep_probes)

def make_farm_files(probe_numbers, dock_path):
    farm_input = '''conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           receptor.sph
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       grid
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                yes
fms_score_use_ref_mol2                                       yes
fms_score_ref_mol2_filename
fms_score_write_reference_pharmacophore_mol2                 yes
fms_score_write_reference_pharmacophore_txt                  yes
fms_score_reference_output_mol2_filename
fms_score_reference_output_txt_filename
fms_score_write_candidate_pharmacophore                      no
fms_score_write_matched_pharmacophore                        no
fms_score_compare_type                                       overlap
fms_score_full_match                                         yes
fms_score_match_rate_weight                                  5.0
fms_score_match_dist_cutoff                                  1.0
fms_score_match_proj_cutoff                                  0.7071
fms_score_max_score                                          20
minimize_ligand                                              yes
simplex_max_iterations                                       1000
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_secondary_minimize_pose                              yes
use_advanced_secondary_simplex_parameters                    no
simplex_secondary_max_iterations                             100
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                path/vdw_AMBER_parm99.defn
flex_defn_file                                               path/flex.defn
flex_drive_file                                              path/flex_drive.tbl
pharmacophore_defn_file                                      path/ph4.defn
ligand_outfile_prefix                                        output
write_orientations                                           no
num_primary_scored_conformers_rescored                       1
num_secondary_scored_conformers                              1
rank_primary_ligands                                         no
rank_secondary_ligands                                       no'''

    os.system('mkdir pharm_farm'); os.chdir('pharm_farm')
    os.system('touch receptor.sph'); os.system('touch grid.bmp'); os.system('touch grid.nrg')
    dock_commands = []
    for number in probe_numbers:
        tmp_farm_input = farm_input
        tmp_farm_input = tmp_farm_input.replace('ligand_atom_file', 'ligand_atom_file                                             ../probe_clusters/cluster'+number+'.mol2')
        tmp_farm_input = tmp_farm_input.replace('fms_score_ref_mol2_filename', 'fms_score_ref_mol2_filename                                  ../probe_clusters/cluster'+number+'.mol2')
        tmp_farm_input = tmp_farm_input.replace('fms_score_reference_output_mol2_filename', 'fms_score_reference_output_mol2_filename                     ../probe_clusters/cluster'+number+'pharm.mol2')
        tmp_farm_input = tmp_farm_input.replace('fms_score_reference_output_txt_filename', 'fms_score_reference_output_txt_filename                      ../probe_clusters/cluster'+number+'pharm.txt')
        tmp_farm_input = tmp_farm_input.replace('path', dock_path+'/parameters')
        inputfile = open('mkpharm'+number+'.in','w+')
        inputfile.write(tmp_farm_input); inputfile.close()
        dock_commands.append('dock6 -i mkpharm'+number+'.in -o mkpharm'+number+'.out')
    return(dock_commands)

# transform cross into a list; open ftmapfile and turn it into a list
if args.cross is None: cross = ['all']
else: cross = args.cross.split(',')

ftmap = open(ftmap_file, 'r').read().splitlines();
logit('\nBegin docking steps:')
mssg = '> FTmap file used: '+ftmap_file+'\n'; logfile.write(mssg)
mssg = '> FTMap cross clusters used: ' + ','.join(map(str, cross))+'\n'
print(mssg); logfile.write(mssg)
logfile.write('> Ligand file: '+args.lig_file+'\n')
# get index where every structure starts (ptn ou crosscl):
ind = np.where(np.asarray([i[:6] for i in ftmap])=='HEADER')[0]

### separate protein if mol2 receptor not given by user:
if args.receptor is None:
    print('Writing receptor.mol2')
    receptor = ftmap[:ind[1]]; receptor = '\n'.join(map(str, receptor))
    recfile = open('receptor.pdb','w+'); recfile.write(receptor); recfile.close() # save receptor.pdb for conversion
    chimera_py = '''import chimera
from DockPrep import prep
models = chimera.openModels.list(modelTypes=[chimera.Molecule])
prep(models)
from WriteMol2 import writeMol2
writeMol2(models, "receptor.mol2")'''
    # ^ python script provided to chimera to run dockprep
    chimera_in = open('dockprep.py', 'w+'); chimera_in.write(chimera_py); chimera_in.close() # write the script
    os.system('chimera --nogui receptor.pdb dockprep.py') # chimera command
    os.system('rm dockprep.py dockprep.pyc') # delete the script and unecessary file created by chimera
    rec_file = 'receptor.mol2'
    logfile.write('> Receptor generated by chimera dockprep and saved as receptor.mol2\n' )
# save rec_file as generated receptor or save rec_file as generated receptor or
else:
    rec_file = args.receptor
    logfile.write('> Receptor provided by user: '+args.receptor+'\n')

### separate probes
probes = ftmap[ind[1]:];
sep_probes = separate_probes(probes, cross)
probe_numbers= [number for number, cluster in sep_probes]

### convert separate probe clusters to mol2 with pymol
os.system('mkdir probe_clusters'); os.chdir('probe_clusters')
print('Writing probes in pdb and mol2 format in ./probe_clusters')
for probenum, probe in sep_probes:
    probemol = open('cluster'+probenum+'.pdb', 'w+'); probemol.write(probe); probemol.close()
    pymol_command = "pymol -c -d 'load cluster"+probenum+".pdb'"+" -d 'save cluster"+probenum+".mol2'"
    os.system(pymol_command)
os.chdir('..')
logfile.write("> Wrote probe clusters' pdb file in ./probe_clusters and converted to mol2 with pymol\n")

### pharmacophore farm
dock_commands = make_farm_files(probe_numbers, args.dock_path)
# run dock to get bundle pharmacophores
print('Running dock6 to get probe pharmacophores. Folder ./pharm_farm will be created then deleted')
for command in dock_commands:
    os.system(command) # pharmacophore files were created in ../probe_clusters/
# done with pharm_farm, can delete:
os.chdir('..'); os.system('rm -r pharm_farm')
logfile.write("> Used dock6 to get probe pharmacophores. Folder ./pharm_farm was created then deleted\n")


##################################################
### convert pharmacophore txt files to spheres ###
##################################################

# functions
def readpharms(pharm_files):
    pharm = pd.read_csv(pharm_files[0], skipinitialspace=True, sep=' ', skiprows=[0,1,2])
    for i in pharm_files[1:]:
        pharm=pharm.append(pd.read_csv(i, skipinitialspace=True, sep=' ', skiprows=[0,1,2]))
    return(pharm)

def formatstr(i):
    string='{:>5}{:>10}{:>10}{:>10}{:>8}{:>5}{:>2}{:>3}'.format(i[4],i[0],i[1],i[2],i[3],i[4],i[6],i[7])
    return(string)

def make_sphgen():
    chimera_in='''from chimera import runCommand, openModels, MSMSModel
# generate surface using 'surf' command
runCommand("surf")
# get the surf object
surf = openModels.list(modelTypes=[MSMSModel])[0]
# write DMS
from WriteDMS import writeDMS
writeDMS(surf, "receptor.ms")'''
    chimeradms = open('ms.py', 'w+'); chimeradms.write(chimera_in); chimeradms.close()
    os.system('chimera --nogui receptor.mol2 ms.py')
    os.system('rm ms.py ms.pyc')
    sphgen_in='''receptor.ms
R
X
0.0
4.0
1.4
spheres.sph'''
    sphgen = open('INSPH', 'w+'); sphgen.write(sphgen_in); sphgen.close()
    logfile.write('> Flag --sphgen used.\n> Created receptor.dms with chimera\n> Wrote INSPH:\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    logfile.write(sphgen_in+'\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')

def chimera_sph_ref():
    chim_in='open probe_clusters/cluster*.pdb; combine #* close true modelId 0; write format mol2 #0 sph_ref.mol2'
    a = open('chim.in', 'w+'); a.write(chim_in); a.close()
    os.system('chimera --nogui < chim.in')
    os.system('rm chim.in')

def make_showsphere(sphere_file, args):
    if args.cross is None and args.sphgen is not None:
        sphcl = args.sphgen
    else:
        sphcl = '1'
    showsph_in=sphere_file+"\n"+sphcl+"\nn\n"+sphere_file[:-3]+"pdb"
    showsph_file = open('showsphere.in','w+'); showsph_file.write(showsph_in); showsph_file.close()

if args.sphgen is not None:
    args.nochem = True
    make_sphgen(); os.system('sphgen')
    logfile.write('> Ran sphgen\n')
    if cross == ['all']:
        sphere_file = 'spheres.sph'
        logfile.write('> No crossclusters specified as reference for sphere_selector.\n')
        if args.sphgen == '0':
            logfile.write('> Argument sphgen set to 0. Will use all spheres for blind docking\n')
        elif args.sphgen == '1':
            logfile.write('> Argument sphgen set to 1. Will use sphere cluster 1 for docking.\n')
    else:
        chimera_sph_ref()
        os.system('sphere_selector spheres.sph sph_ref.mol2 3')
        logfile.write('> Ran sphere_selector with specified crossclusters as reference. Command: sphere_selector spheres.sph sph_ref.mol2 3')
        sphere_file = 'selected_spheres.sph'
else:
    ## read all pharmacophores into a dataframe
    pharm_files = ['./probe_clusters/cluster'+i+'pharm.txt' for i in probe_numbers]
    pharm = readpharms(pharm_files)
    ## rearrange dataframe into sphere format
    sph=pharm.copy(); sph=sph.reset_index(drop=True)
    naX=  ['1' for i in range(len(sph))] # atom number columns
    crit= [0 for i in range(len(sph))] # critical column
    sph['coox'] = [str(i)+'0'*(8-len(str(i))) for i in sph['coox']] # make columns x, y, z and
    sph['cooy'] = [str(i)+'0'*(8-len(str(i))) for i in sph['cooy']] # format them into
    sph['cooz'] = [str(i)+'0'*(8-len(str(i))) for i in sph['cooz']] # 8 algarisms
    sph['radius'] = sph['radius'].round(3)
    sph['radius'] = [str(i)+'0'*(5-len(str(i))) for i in sph['radius']]
    # transformar o ph4name nos numeros das esferas do chemical matching:
    # transform ph4name into chemical matching sphere numbers:
    mapping = {'PHO': 1, 'HBD': 2, 'HBA':3, 'ARO':4, 'RNG':5, 'NEG':6,'POS':7}
    chem=sph['ph4name'].replace(mapping)
    sph=sph.drop(['ph4id','vecx','vecy','vecz','ph4name'], axis=1) # remove unnecessary columns
    # insert created columns
    sph['na1']=naX; sph['na2']=naX; sph['crit']=crit; sph['chem']=chem;
    ## create sphere file
    sph_header= """DOCK 3.5 receptor_spheres
color hydrophobic 1
color donor 2
color acceptor 3
color aromatic 4
color aroAcc 5
color negative 6
color positive 7
cluster     1   number of spheres in cluster   """+str(len(sph))+'\n' # in order: PHO, HBD, HBA, ARO, RNG, NEG, POS
    # open sphere file and write header
    sphfile = open('spheres_chem.sph', 'a+')
    sphfile.write(sph_header)
    sphere_file = 'spheres_chem.sph'
    # format sphere each line
    print("Converting probes' pharmacophores into dock6 sphere file: spheres_chem.sph")
    for line in range(len(sph)):
        i=[line for line in sph.loc[line]]
        sphfile.write(formatstr(i)+'\n')
    sphfile.close()
    logfile.write("> Wrote spheres_chem.sph with chemical information from probe pharmacophores\n")

### using showsphere to convert spheres_chem.sph to spheres_chem.pdb
make_showsphere(sphere_file, args)
os.system('showsphere < showsphere.in')
logfile.write("> Used dock tool showsphere to convert spheres_chem.sph to pdb to allow visualization of the spheres\n")

##################################################
###                 docking                    ###
##################################################

# functions
def make_params():
    chem_custom = '''______________________________________________________________________
name		null
definition	*
______________________________________________________________________
name		hydrophobic

definition	C. [ O. ] [ N. ] [ S. ] [ F ] [ P ] ( * )
definition      C. ( N.pl3 ( 2 C. ) ) ( * )
definition	N.pl3 ( 3 C. )
______________________________________________________________________
name		donor

definition      H ( O. )
definition      H ( N. )
definition      H ( S. )
definition      H ( F )
______________________________________________________________________
name		acceptor

definition	O. ( * )
definition      N.1 ( 1 * )
definition      N.2 [ 3 * ]
definition      N.3 ( 3 * )
definition      N.pl3 ( 2 * ) [ H ]
definition      S.2 [ O. ] [ N. ]
definition      S.3 ( 2 * )
definition      F  ( * )
definition      Cl ( * )
______________________________________________________________________
name            aromatic


definition      C.ar
definition      N.ar
______________________________________________________________________
name            aroAcc

definition      N.ar [ H ] [ 3 * ] ( * )
______________________________________________________________________
name            negative

definition      C. ( 2 O.co2 )
definition      C.2 ( O.2 ) ( O.3 [ * ] )
definition      P. ( 4 O. ) ( O.3 [ * ] )
definition      S. ( 3 O. ) ( O.3 [ * ] )
definition      S. ( 4 O. ) ( O.3 [ * ] )
definition      F [ * ]
definition      Cl [ * ]
______________________________________________________________________
name            positive

definition      C.cat ( * )
definition      N.4 ( * )
definition      N.3 ( 4 * )
definition      N.2 ( 3 * )
definition      Zn [ * ]
definition      Mg [ * ]
definition      Ca [ * ]
definition      Mn [ * ]
definition      K  [ * ]
definition      Fe [ * ]
______________________________________________________________________
'''
    chem_match='''label	null
label	hydrophobic
label	donor
label	acceptor
label   aromatic
label   aroAcc
label	negative
label   positive

table
1
1	1
1	0	1
1	0	0	1
1	0	0	0	1
1	0	0	0	0	1
1	0	0	0	0	0	1
1	0	0	0	0	0	0	1'''
    os.system('mkdir parameters'); os.chdir('parameters')
    a = open('chem_custom.defn','w+'); a.write(chem_custom); a.close()
    b = open('chem_match_custom.tbl', 'w+'); b.write(chem_match); b.close()
    os.chdir('..')

def make_box(args):
    if args.boxmarg is None: margin = '5'
    else: margin = args.boxmarg
    if args.cross is None and args.sphgen is not None: clnum = args.sphgen
    else: clnum = '1'
    mssg = '> Box margin: '+margin+' angstroms\n'; print(mssg); logfile.write(mssg)
    box_in = 'y\n'+margin+'\n'+sphere_file+'\n'+clnum+'\nbox.pdb'; open('box.in','w+').write(box_in)
    os.system('showbox < box.in')

def make_grid_input(args):
    if args.receptor is None: receptor = 'receptor.mol2'
    else: receptor = args.receptor
    grid_in = '''compute_grids                  yes
grid_spacing                   0.3
output_molecule                no
contact_score                  no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            6
repulsive_exponent             12
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  '''+receptor+'''
box_file                       box.pdb
vdw_definition_file            '''+args.dock_path+'''/parameters/vdw_AMBER_parm99.defn
score_grid_prefix              grid'''
    a = open('grid.in', 'w+'); a.write(grid_in); a.close()
    return(grid_in)

def make_dock_input():
    dock_in = '''conformer_search_type                                        flex
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              5
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               100.0
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            yes
clash_overlap                                                0.5
write_growth_tree                                            no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             X
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           X
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            X
chem_match_tbl                                               X
use_ligand_spheres                                           no
bump_filter                                                  yes
bump_grid_prefix                                             grid
max_bumps_anchor                                             12
max_bumps_growth                                             12
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         yes
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       grid
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  500
simplex_grow_tors_premin_iterations                          0
simplex_secondary_minimize_pose                              yes
use_advanced_secondary_simplex_parameters                    no
simplex_secondary_max_iterations                             100
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                X
flex_defn_file                                               X
flex_drive_file                                              X
chem_defn_file                                               X
ligand_outfile_prefix                                        X
write_orientations                                           no
num_primary_scored_conformers_rescored                       10000
write_primary_conformations                                  no
cluster_primary_conformations                                yes
cluster_rmsd_threshold                                       2.0
num_clusterheads_for_rescore                                 5
num_secondary_scored_conformers                              10000
write_secondary_conformations                                no
rank_primary_ligands                                         yes
max_primary_ranked                                           10000
rank_secondary_ligands                                       yes
max_secondary_ranked                                         10000'''
    return(dock_in)


# make chemical matching parameter folder
print('Creating ./parameter folder and parameter files for chemical_matching')
make_params()
logfile.write("> Created ./parameter folder and parameter files for docking with chemical matching\n")

### make grid
print('Creating box.in and running showbox to create box around spheres')
make_box(args) # make box
logfile.write('> Created box.\n')
if args.grid == True:
    print('Using preexisting grid files')
    logfile.write("> Did not run grid. User provided preexisting grid files.\n")
else:
    grid_in = make_grid_input(args) # make grid.in
    # run dock6 grid
    print('Running dock6 grid')
    start = datetime.datetime.now()
    os.system('grid -i grid.in -o grid.out')
    end = datetime.datetime.now(); elaps = end - start
    logfile.write("> Ran dock tool grid. Elapsed time: "+str(round(elaps.seconds/60, 2))+" minutes\n")
    logfile.write("> grid.in parameters:\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"+grid_in+"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")

### dock
# create input var from default or user provided input file
if args.input is None:
    dock_in = make_dock_input(); dock_in = dock_in.splitlines()
    logfile.write("> Created default dock.in input file.\n")
else:
    dock_in = open(args.input, 'r').read().splitlines()
    logfile.write("> User provided dock input file. Did not generate default.\n")

if args.nochem == True: match = 'no'
else: match = 'yes'
logfile.write('> Docking to be performed with chemical_matching: '+match+'\n')

flags = [('ligand_atom_file' ,      'ligand_atom_file                                             '+args.lig_file),
         ('receptor_site_file' ,    'receptor_site_file                                           '+sphere_file),
         ('chemical_matching' ,     'chemical_matching                                            '+match),
         ('chem_match_tbl' ,        'chem_match_tbl                                               ./parameters/chem_match_custom.tbl'),
         ('vdw_defn_file' ,         'vdw_defn_file                                                '+args.dock_path+'/parameters/vdw_AMBER_parm99.defn'),
         ('flex_defn_file' ,        'flex_defn_file                                               '+args.dock_path+'/parameters/flex.defn'),
         ('flex_drive_file' ,       'flex_drive_file                                              '+args.dock_path+'/parameters/flex_drive.tbl'),
         ('chem_defn_file' ,        'chem_defn_file                                               ./parameters/chem_custom.defn'),
         ('ligand_outfile_prefix' , 'ligand_outfile_prefix                                        ligand')]

# make dock.in
print('> Creating dock input file dock.in')
for flag, value in flags:
    exists = False
    for i in range(len(dock_in)):
        line = dock_in[i]
        if line[:len(flag)]==flag:
            exists = True
            dock_in[i]=value
    if exists == False:
        dock_in.append(value)
dock_in = '\n'.join(map(str, dock_in))
dock_infile = open('dock.in', 'w+'); dock_infile.write(dock_in); dock_infile.close()
logfile.write("> Updated dock input parameters and created dock.in:\n\n")
logfile.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
logfile.write(dock_in+'\n')
logfile.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")

# run dock6
print('> Running dock6')
start = datetime.datetime.now()
os.system('dock6 -i dock.in -o dock.out')
end = datetime.datetime.now(); elaps = end - start
logfile.write("Ran dock6. Command: dock6 -i dock.in -o dock.out\n")
if elaps.seconds > 60: timed = str(round(elaps.seconds/60, 2))+" minutes\n"
else: timed = str(elaps.seconds)+" seconds\n"
logfile.write("Elapsed time for docking: "+timed)
logfile.write("Done :D\n\n\n"); logfile.close()
print('Done!')














### area de teste
