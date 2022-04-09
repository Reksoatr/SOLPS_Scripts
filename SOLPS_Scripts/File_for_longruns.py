import os
#need to fix the script later
def batch_extended(lib, start, end):
    f = open('batch_use_extended', 'w')
    f.writelines(['#!/bin/tcsh','\n#PBS -l nodes=1:hima:ppn=1','\n#PBS -l walltime=24:00:00','\n#PBS -N LARGE','\n#PBS -j oe','\n','\nenv','\n','\n'])
    f.writelines(['\n python3 -m steepest_descent_analysis <<Inputblock','\n1','\nn', '\nInputblock'])
    f.writelines([f'\ncd /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_1\n','\nb2run b2mn > run.log'])
    f.writelines(['\ncd /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}'])
    for i in range(start,end):
        f.writelines(['\n python3 -m steepest_descent_analysis <<Inputblock',f'\n{i}','\ny', '\nInputblock'])
        f.writelines(['\n python3 -m steepest_descent_analysis <<Inputblock',f'\n{i}','\nn', '\nInputblock'])
        f.writelines(['\n',f'cd /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{i}\n','\nb2run b2mn > run.log'])
        f.writelines(['\ncd /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}'])
    f.close()
if __name__ == '__main__':
    blip = 3
    batch_extended(blip, 2,15)
    os.system('qsub batch_use_extended')