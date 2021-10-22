'''
run poar setup in order to create scenario files
'''
import sys
from par import setup_params as sepa

if __name__ =='__main__':
    if len(sys.argv) <2:
        print('Using default super_parameters: par_v002.json')
        filename_super_parameters = 'par_v002.json' # Should be user input
    else:
        print('Using super_parameters: ', str(sys.argv[-1]))
        filename_super_parameters = sys.argv[-1]

    sepa.ScenarioSetup(filename_super_parameters)
