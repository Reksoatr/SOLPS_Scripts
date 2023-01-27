import os
#need to fix the script later
def batch_extended(lib, start=1, end=25, alg='T', path = '/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_'):
    '''Function that creates a batch script running multiple runs on an HPC.
    Requires a baserun file and a /base 500 step run with a transport coefficient guess and corresponding 
    params.txt file. A g-file will also be required.
    (find examples from different machines in the literature)
    
    
    Parameters
    ----------
    lib : library number where you are running this attempt
    
    path: path name found using pwd in the optimization
    /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_ should be replaced with the rest of your directory

    start: starting attempt (set to 1 unless continuing from previous iterations)

    end: last attempt (set to 25 for a 24 hour job)
    
    alg: algorithm being used for fitting the transport profile. The single Gaussian was determined
    to be the best, so it is recommended for users. The double gauss is provided for users as an option as well,
    and users can edit train_test2 with other algorithms.
    
    Returns batch_use_extended file and runs it using qsub
    -------
    '''
    f = open('batch_use_extended', 'w')
    f.writelines(['#!/bin/tcsh','\n#PBS -l nodes=1:hima:ppn=1','\n#PBS -l walltime=24:00:00',f'\n#PBS -N LARGE_{lib}','\n#PBS -j oe','\n','\nenv','\n',f'\ncd {path}{lib}/'])
    f.writelines([f'\npython3 -m {alg} <<Inputblock','\n1','\nn',f'\n{lib}', '\nInputblock'])
    f.writelines([f'\ncd {path}{lib}/Attempt_1\n','\nb2run b2mn > run.log'])
    f.writelines([f'\ncd {path}{lib}'])
    for i in range(start,end):
        f.writelines([f'\npython3 -m {alg} <<Inputblock',f'\n{i}','\ny',f'\n{lib}', '\nInputblock'])
        f.writelines([f'\npython3 -m {alg} <<Inputblock',f'\n{i+1}','\nn',f'\n{lib}', '\nInputblock'])
        f.writelines(['\n',f'cd {path}{lib}/Attempt_{i+1}\n','\nb2run b2mn > run.log'])
        f.writelines([f'\ncd {path}{lib}'])
    f.close()
if __name__ == '__main__':
    blip = int(input('Which Directory Number?'))
    #other algorithms can be added here
    functtt = input('T (Trainer) or D (DoubleGauss)')
    if functtt == 'T':
        batch_extended(blip, 1,25,alg='steepest_descent')        
    elif functtt == 'D':
        batch_extended(blip, 1,25,alg='steepest_descent_dblgauss')
    elif functtt == 'T2':
        batch_extended(blip, 1,25,alg='train_test2')        
    else:
        print('Invalid input')
    os.system('qsub batch_use_extended')