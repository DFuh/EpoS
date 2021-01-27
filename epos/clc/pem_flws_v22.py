'''
materialbalance PEM
'''
print(__name__ + ' imported...')

from collections import namedtuple
import epos.aux.faux as fx

xflws = fx.dyn_aux_import(__file__, __name__)


def materialbalance(obj, T, i, m_H2O_in_an, p_an, p_ca, c_in, n_in):#T, i, c_in, n_in): #sns=False):
    '''
    materialbalance PEM
    flow balance in cell

    // in mol/s
    '''
    # - water diffusion through membrane disabled

    # TODO: water-balance
    # TODO: check density and viscosity calc. !
    # TODO: check permeation (especially O2!)
    # TODO: rho_H2O
    # TODO: D_eff_
    # TODO: check concentrations of permeated gases at membrane !
    # TODO: what about pp_H2O (gaseous)
    A_cell = obj.pcll.active_cell_area

    xflws.clc_flws_auxpars(obj, T )

    ### production // in mol/s
    n_H2_prd, n_O2_prd, n_H2O_cns = xflws.clc_mlr_flws_prod(obj, obj.pec, i, A_cell)

    ### permeation
    #n_H2O_eo, n_H2O_dd, n_H2O_pe =xflws.clc_crssflw_membrane_H2O(obj,obj.pec, i, A_cell)
    n_H2O_prm =xflws.clc_crssflw_membrane_H2O_chandesris(obj,obj.pec, i, T, A_cell)

    n_H2_prm = xflws.clc_hydrogen_permeation(obj, obj.pec, i)

    n_O2_prm = xflws.clc_oxygen_permeation(obj, )

    ### H2O
    # conc in channels
    c_H2O_ch = obj.av.rho_H2O / obj.pec.M_H2O
    #n_H2O_prm  = n_H2O_eo + n_H2O_dd + n_H2O_pe

    n_H2O_ch_in_an  =  m_H2O_in_an / obj.pec.M_H2O # Convert massflow (kg/s) to molar flow (mol/s)
    n_H2O_ch_out_an = n_H2O_ch_in_an - (n_H2O_prm + n_H2O_cns)
    n_H2O_ch_in_ca  = 0
    n_H2O_ch_out_ca = n_H2O_ch_in_ca + n_H2O_prm

    # H2
    n_H2_ch_out_an    = n_H2_prm
    #n_H2_ch_in_ca = 0
    n_H2_ch_out_ca    = n_H2_prd - n_H2_prm #

    # O2
    n_O2_ch_out_an    = n_O2_prd - n_O2_prm
    #n_H2_ch_in_ca = 0
    n_O2_ch_out_ca    = n_O2_prm #

    # CHEck !!!!
    n_H2O_ca = (n_H2O_prm + n_H2O_cns)
    n_H2O_an = n_H2O_ch_out_an

    ### molar fractions
    '''
    calc. molar fractions of species in outlet-channels of cell
    Marangio_2009 eq.: 37
    '''
    x_H2_ch_ca = xflws.clc_mlr_frc(n_H2_ch_out_ca,
                            [n_H2_ch_out_ca, n_O2_ch_out_ca, n_H2O_ch_out_ca])
    #n_H2_ch_ca / (n_H2_ch_ca + n_O2_ch_ca + n_H2O_ch_ca)

    x_H2_ch_an = xflws.clc_mlr_frc(n_H2_ch_out_an,
                            [n_H2_ch_out_an, n_O2_ch_out_an, n_H2O_ch_out_an])
    #n_H2_ch_an / (n_H2_ch_an + n_O2_ch_an + n_H2O_ch_an)

    x_O2_ch_ca = xflws.clc_mlr_frc(n_O2_ch_out_ca,
                            [n_H2_ch_out_ca, n_O2_ch_out_ca, n_H2O_ch_out_ca])
    #n_O2_ch_ca / (n_H2_ch_ca + n_O2_ch_ca + n_H2O_ch_ca)

    x_O2_ch_an = xflws.clc_mlr_frc(n_O2_ch_out_an,
                            [n_H2_ch_out_an, n_O2_ch_out_an, n_H2O_ch_out_an])
    #n_O2_ch_an / (n_H2_ch_an + n_O2_ch_an + n_H2O_ch_an)

    # TODO: check n_H2O_an/ca !!!
    x_H2O_ch_ca = n_H2O_ca / (n_H2_ch_out_ca + n_O2_ch_out_ca + n_H2O_ca)
    x_H2O_ch_an = n_H2O_an / (n_H2_ch_out_an + n_O2_ch_out_an + n_H2O_an)

    ### concentrations
    '''
    clc concentration based on pressure (in compartment)
    and molar fraction

    Concentration of species in channels
    Marangio 2009 eq. 38
    '''
    # -> in channels
    c_H2_ch_ca = p_ca * x_H2_ch_ca / (obj.pec.R * T)
    c_H2_ch_an = p_an * x_H2_ch_an / (obj.pec.R * T)

    c_O2_ch_ca = p_ca * x_O2_ch_ca / (obj.pec.R * T)
    c_O2_ch_an = p_an * x_O2_ch_an / (obj.pec.R * T)

    c_H2O_ch_an = c_H2O_ch_ca = obj.av.rho_H2O / obj.pec.M_H2O

    # -> at membrane/ electrode
    '''
    Concentration of species at membrane
    Marangio 2009 eq. 39
    '''
    c_H2_mem_ca = (c_H2_ch_ca
                    + (obj.pec.d_el_ca / obj.av.D_eff_H2 /A_cell) * n_H2_prd)
    #TODO: check H2 at Anode
    c_H2_mem_an = (c_H2_ch_an
                    + (obj.pec.d_el_an / obj.av.D_eff_H2 /A_cell) * n_H2_ch_out_an)
    #TODO: check O2 at Cathode
    c_O2_mem_ca = (c_O2_ch_ca
                    + (obj.pec.d_el_ca / obj.av.D_eff_O2 /A_cell) * n_O2_ch_out_ca)
    c_O2_mem_an = (c_O2_ch_an
                    + (obj.pec.d_el_an / obj.av.D_eff_O2 /A_cell) * n_O2_prd)

    #c_H2O_mem = obj.av.rho_H2O / obj.pec.M_H2O
    # TODO: check D_eff for H2O
    # liquid water?
    c_H2O_mem_ca = (c_H2O_ch_ca
                    + (obj.pec.d_el_ca/obj.av.D_eff_H2 /A_cell) * n_H2O_ca )
    c_H2O_mem_an = (c_H2O_ch_an
                    + (obj.pec.d_el_an/obj.av.D_eff_O2 /A_cell) * n_H2O_an) # Check !

    ### clc molar fractions at membrane
    #TODO: check, if correct !
    x_H2_mem_ca = c_H2_mem_ca * (obj.pec.R * T) / p_ca
    x_H2_mem_an = c_H2_mem_an * (obj.pec.R * T) / p_an

    x_O2_mem_ca = c_O2_mem_ca * (obj.pec.R * T) / p_ca
    x_O2_mem_an = c_O2_mem_an * (obj.pec.R * T) / p_an

    x_H2O_mem_ca = c_H2O_mem_ca / (c_H2O_mem_ca + c_H2_mem_ca)
    x_H2O_mem_an = c_H2O_mem_an / (c_H2O_mem_an + c_O2_mem_an)

    ### clc partial pressures
    pp_H2_mem_ca = x_H2_mem_ca *p_ca
    pp_H2_mem_an = x_H2_mem_an *p_an
    pp_O2_mem_ca = x_O2_mem_ca *p_ca
    pp_O2_mem_an = x_O2_mem_an *p_an

    pp_H2O_mem_an = x_H2O_mem_an *p_an


    #TT = namedtuple('TT', ['test_a', 'test_aa','test_b', 'test_bb'])
    FLWS = namedtuple('FLWS', '''n_H2_out_ca n_H2_out_an n_O2_out_ca n_O2_out_an
                                c_H2_out_ca c_H2_out_an c_O2_out_ca c_O2_out_an
                                x_H2_out_ca x_H2_out_an x_O2_out_ca x_O2_out_an
                                pp_H2_mem_ca pp_H2_mem_an pp_O2_mem_ca pp_O2_mem_an,
                                pp_H2O_mem_an ''')

    flws_out = FLWS(n_H2_ch_out_ca, n_H2_ch_out_an, n_O2_ch_out_ca, n_O2_ch_out_an,
                    c_H2_ch_ca, c_H2_ch_an, c_O2_ch_ca, c_O2_ch_an,
                    x_H2_ch_ca, x_H2_ch_an, x_O2_ch_ca, x_O2_ch_an,
                    pp_H2_mem_ca, pp_H2_mem_an, pp_O2_mem_ca, pp_O2_mem_an,
                    pp_H2O_mem_an )
    #tt = TT(1001,1,2002,2)
    return flws_out #tt#c_out, n_out
