import argparse
import os
import pdb
import shutil
import sys
import intermol.unit as units
import subprocess

#--------DESMOND energy evaluation methods---------#
def parse_args():
    parser = argparse.ArgumentParser(description = 'testing file conversions')
    parser.add_argument('-i', '--in', dest='infile', help="input struture/topology .cms file")
    parser.add_argument('--cfg', dest='cfgfile', default='Inputs/Desmond/onepoint.cfg', help="input .cfg file (for DESMOND)")
    # to do: add arguments for other file types
    args = parser.parse_args()
    return args

def get_desmond_energy_from_file(energy_file):
    ''' 
    parses the desmond energy file
    for now, returns just the total energy
    '''
    with open(energy_file, 'r') as f:
        for line in f:
            if 'Total' in line:
                tot_energy = line.split()[-1]
                break
    return tot_energy

def desmond_energies(cms, cfg, despath='/opt/schrodinger2013'):
    """
    Evalutes energies of DESMOND files
    Args:
        cms = cms file
        cfg = cfg file
        despath = path to DESMOND binaries

    """
    cms = os.path.abspath(cms)
    cfg = os.path.abspath(cfg)
    direc, cms_filename = os.path.split(cms)
    cwd = os.getcwd()
    name = 'system'
    energy_file = '%s/%s.enegrp.dat' % (direc, name)

    # first see if the file already exists
    if os.path.exists(energy_file):
        print '%s already exists, not running DESMOND' % energy_file
        tot_energy = get_desmond_energy_from_file(energy_file)
        return tot_energy, energy_file

    # use DESMOND To evaluate energy
    #    cd to directory of cms file so that files generated by desmond
    #    don't clog the working directory
    os.chdir(direc)   
    if os.path.exists('trj'):
        shutil.rmtree('trj')
    cmd = '{despath}/desmond -WAIT -P 1 -in {cms} -JOBNAME {name} -c {cfg}'.format(despath = despath, name = name, cms = cms, cfg = cfg)
    print 'Running DESMOND with command'
    print cmd
    exit = os.system(cmd)
    if exit: # exit status not 0
        print 'Failed evaluating energy of {0}'.format(cms)
        os.chdir(cwd)
        sys.exit(1)
    
    # parse desmond energy file
    os.chdir(cwd)
    tot_energy = get_desmond_energy_from_file(energy_file)
    return tot_energy, energy_file

#--------GROMACS energy evaluation methods---------#
# to do: clean up

def gromacs_energies(name, top=None, gro=None, in_out='in', gropath='',grosuff=''):
    """

    gropath = path to gromacs binaries
    grosuff = suffix of gromacs binaries, usually '' or '_d'

    """
    mdp = 'Inputs/Gromacs/grompp.mdp'

    if in_out == 'in':
        base = 'Inputs/Gromacs'
        if top == None:
            top = os.path.join(base, name, 'topol.top')
        if gro == None:
            gro = os.path.join(base, name, 'conf.gro')
    elif in_out == 'GtoG':
        base = 'Outputs/GromacsToGromacs'
        if top == None:
            base = os.path.join(base, name, 'topol.top')
        if gro == None:
            base = os.path.join(base, name, 'conf.gro')
    elif in_out == 'LtoG':
        base = 'Outputs/LammpsToGromacs'
        if top == None:
            base = os.path.join(base, name, 'topol.top')
        if gro == None:
            base = os.path.join(base, name, 'conf.gro')
    else:
        raise Exception("Unknown flag: {0}".format(in_out))

    tpr  = os.path.join(base, name, 'topol.tpr')
    ener  = os.path.join(base, name, 'ener.edr')
    ener_xvg  = os.path.join(base, name, 'energy.xvg')
    conf  = os.path.join(base, name, 'confout.gro')
    mdout = os.path.join(base, name, 'mdout.mdp')
    state  = os.path.join(base, name, 'state.cpt')
    traj  = os.path.join(base, name, 'traj.trr')
    log  = os.path.join(base, name, 'md.log')

    grompp_bin = os.path.join(gropath, 'grompp' + grosuff)
    mdrun_bin = os.path.join(gropath, 'mdrun' + grosuff)
    genergy_bin = os.path.join(gropath, 'g_energy' + grosuff)

    # grompp'n it up
    os.system(grompp_bin + " -f {mdp} -c {gro} -p {top} -o {tpr} -po {mdout} -maxwarn 1".format(mdp=mdp,
            top=top, gro=gro, tpr=tpr, mdout=mdout))

    # mdrunin'
    os.system(mdrun_bin + " -s {tpr} -o {traj} -cpo {state} -c {conf} -e {ener} -g {log}".format(tpr=tpr,
            traj=traj, state=state, conf=conf, ener=ener, log=log))

    # energizin'
    select = " ".join(map(str, range(1, 15))) + " 0 "
    os.system("echo {select} | ".format(select=select) + genergy_bin + " -f {ener} -o {ener_xvg} -dp".format(ener=ener,
            ener_xvg=ener_xvg))

    # extract g_energy output and parse initial energies
    with open(ener_xvg) as f:
        all_lines = f.readlines()

    # take last line
    sec_last = all_lines[-1].split()[1:]
    data = map(float, sec_last)

    # give everything units
    temp = data[-1] * units.kelvin
    data = [value * units.kilojoules_per_mole for value in data[:-1]]
    data.append(temp)

    # pack it all up in a dictionary
    types = ['Bond', 'Angle', 'Proper Dih.', 'Ryckaert-Bell.', 'LJ-14', 'Coulomb-14',
            'LJ (SR)', 'Disper. corr.', 'Coulomb (SR)', 'Coul. recip.', 'Potential',
            'Kinetic En.', 'Total Energy', 'Temperature']
    e_out = dict(zip(types, data))
    return e_out

#--------LAMMPS energy evaluation methods---------#
#to do: clean up

def lammps_energies(name, in_out='in', lmppath='', lmpbin='lmp_openmpi'):
    """Evaluate energies of LAMMPS files

    Args:
        lmppath = path to LAMMPS binaries
        lmpbin = name of LAMMPS binary
    """

    if in_out == 'in':
        base = 'Inputs/Lammps'
    elif in_out == 'GtoL':
        base = 'Outputs/GromacsToLammps'
    elif in_out == 'LtoL':
        base = 'Outputs/LammpsToLammps'
    else:
        raise Exception("Unknown flag: {0}".format(in_out))

    lmpbin = os.path.join(lmppath, lmpbin)
    sim_dir = os.path.join(base, name)
    log = os.path.join(sim_dir, 'log.lammps')

    # mdrunin'
    saved_path = os.getcwd()
    os.chdir(sim_dir)
    run_lammps = "{lmpbin} < data.input".format(lmpbin=lmpbin)
    #run_lammps = "{lmpbin} < input_file.out".format(lmpbin=lmpbin)
    os.system(run_lammps)
    os.chdir(saved_path)

    # energizin'
    proc = subprocess.Popen(["awk '/E_bond/{getline; print}' %s" % (log)],
            stdout=subprocess.PIPE, shell=True)
    (energies, err) = proc.communicate()

    data = map(float, energies.split())

    # give everything units
    #temp = data[-1] * units.kelvin
    data = [value * units.kilocalories_per_mole for value in data]
    #data.append(temp)

    # pack it all up in a dictionary
    types = ['Bond', 'Angle', 'Proper Dih.', 'Improper', 'Pairs', 'vdW', 'Coulomb', 'Potential']

    e_out = dict(zip(types, data))
    return e_out

def main():
    args = parse_args()
    energy, energy_file = desmond_energies(args.infile, args.cfgfile)
    print 'Total energy from %s:' % energy_file
    print energy
    
if __name__ == '__main__':
    main()