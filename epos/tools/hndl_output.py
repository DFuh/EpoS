'''
handle output data
'''
from dataclasses import dataclass

@dataclass
class Data():
    '''
    Container for calculation values
    '''
    t_abs:      float = # Absolute time simulated
    t_d:        float = # Time difference of calc.

    T_st:       float = # Stack temperature     # // in K
    m_c:        float = # Mass floww of coolant // in kg/s

    i:          float = # Current density
    u_cell:     float = # in V
    u_anode:    float =     # in V
    u_cathode:  float = # in V
    u_dgr:      float = # in V

    P_i:        float =
    P_act:      float =
    P_st:       float =
    P_rect:     float =
    P_aux:      float =       # in W

    p_anode:    float = #
    p_cathode:  float = #                           # in Pa ->

    n_H2_an:    float =
    n_H2_ca:    float =
    n_O2_an:    float =
    n_O2_ca:    float = #
    n_H2_perm:  float = #
    n_O2_perm:  float = #
    n_H2O:      float = #

    V_ely, d_mem, dummy1,dummy2,dummy_cH2,dummy_cO2,dummy_pH2,dummy_pO2,
    x_H2inO2, x_O2inH2)

    def ini_data_output():

        return


'''
def ini_data_output(self):


        - make output-directory
        - prepare csv-files (one per year)

        Returns
        -------
        None.



        ### specify output pth (directory)
        output_pth_basic = self.basepath + '/data/out'
        if
        subdir0 =
        subdir = '/' + str(self.tdd.year) +'/'+ str(self.tdd.month) +'/'+ str(self.tdd.day)

        path_out = output_pth_basic + subdir
        if not os.path.exists(path_out):
                os.makedirs(path_out)
        out_pth, out_flnm = None
        ###

        return
'''
