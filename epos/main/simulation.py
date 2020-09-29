'''
contains basic simualtion class
'''
import os
import datetime
import logging
import uuid

import time

import aux.readingfiles as rf


class ElSim():

    def __init__(self, parameter_filename):
        ### auxilliary parameters
        #logging.basicConfig(filename='example_df.log',level=logging.INFO)

        # date
        self.tdd            = datetime.datetime.now()
        self.today_ymd      = str(self.tdd.year) +str(self.tdd.strftime("%m")) + str(self.tdd.day)
        self.today_ymdhs    = self.today_ymd + str(self.tdd.hour) + str(self.tdd.minute)
        self.cwd = os.getcwd()
        ### Parameters
        self.prms = rf.read_json_file(filename=parameter_filename) # Parameters as dict
        self.name = self.prms['scen_name'].replace('Scen','Sim')
        self.tag = uuid.uuid1()

        ### input data
        self.data_sig, self.metadata_sig = rf.read_in_signal_dataset(self, filename=self.prms['relpth_sig_data'])
        # check properties of df ?
        #logging.info
        print('Finished __init__() of Simulation: ', self.name)


    def check_properties_sig_data():
        '''
        compare properties of data with parameters
        resolve conflicts
        '''

        # see: hndl_data.py !


        return


    def run(self,):
        #logging.info
        print('Starting Simulation: ', self.name)
        time.sleep(1)
        #logging.info
        print('End Simulation: ', self.name)
        return None
