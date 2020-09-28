'''
run poar setup in order to create scenario files
'''
from par import setup_params as sepa

if __name__ =='__main__':
    filename_super_parameters = 'par_v001.json'
    sepa.ScenarioSetup(filename_super_parameters)
