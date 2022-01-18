'''
run poar setup in order to create scenario files
'''
import sys
from epos.par import setup_params as sepa

if __name__ =='__main__':
    if len(sys.argv) <2:
        print('Using default super_parameters: par_v01.json')
        filename_super_parameters = 'par_v01.json' # Should be user input
    else:
        print('Using super_parameters: ', str(sys.argv[-1]))
        fl_nm = sys.argv[-1]
        if '.json' not in fl_nm:
            fl_nm = fl_nm+'.json'
        filename_super_parameters = fl_nm 

    sepa.ScenarioSetup(filename_super_parameters)
