import os
#need to fix the script later
def batch_extended(dest, start, end):
    f = open('batch_use_extended', 'w')
    f.writelines(['#!/bin/tcsh','\n#PBS -l nodes=1:hima:ppn=1','\n#PBS -l walltime=24:00:00','\n#PBS -N LARGE','\n#PBS -j oe','\n','\nenv','\n','\n'])
    for i in range(start,end):
        f.writelines(['\n python3 -m steepest_descent_analysis <<Inputblock',f'\n{i}','\ny', '\nInputblock'])
        f.writelines(['\n python3 -m steepest_descent_analysis <<Inputblock',f'\n{i}','\nn', '\nInputblock'])
        f.writelines(['\n',dest,f'{i}\n','\nb2run b2mn > run.log'])
    f.close()
if __name__ == '__main__':
    blip = 'cd /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_3/Attempt_'
    batch_extended(blip, 2,15)