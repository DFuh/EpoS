'''
contains basic simualtion class
'''


class Simulation():

    def __init__(parameter_file):
        ### auxilliary parameters
        # date
        self.tdd            = datetime.datetime.now()
        self.today_ymd      = str(self.tdd.year) +str(self.tdd.month) + str(stef.tdd.day)
        self.today_ymdhs    = self.today_yms + str(self.tdd.hour) + str(self.tdd.minute)

        ### Parameters
        self.prms = self.read_in_parameter_file(parameter_file) # Parameters as dict

        ### input data
        self.data_sig_input = self.read_in_signal_dataset()
        # check properties of df ?







    def check_properties_sig_data():
        '''
        compare properties of data with parameters
        resolve conflicts
        '''

        # see: hndl_data.py !


        return


    def run(self,):

        return
