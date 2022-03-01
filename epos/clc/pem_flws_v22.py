'''
materialbalance PEM
'''
print(__name__ + ' imported...')

import numpy as np
from collections import namedtuple
import epos.auxf.faux as fx

xflws = fx.dyn_aux_import(__file__, __name__)


def materialbalance(obj, T, i, m_H2O_in_an, p, c_in, n_in,
                    stf=1, sns=False, ntd=None, m=None):#T, i, c_in, n_in): #sns=False):
    '''
    materialbalance PEM
    flow balance in cell
    level: 1 cell

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

    #A_cell = obj.pcll.active_cell_area
    #print(' ----------- PEM_flws -------------------')
    xflws.clc_flws_auxpars(obj, T )
    # print(f'm_H2O_in_an = {m_H2O_in_an}')
    ### production // in mol/(m² s)
    n_H2_prd, n_O2_prd, n_H2O_cns = xflws.clc_mlr_flws_prod(obj, obj.pec, i)#, A_cell)

    #print('i (flws): ', i)
    #print(f'n_ (flws) -> n_H2: {n_H2_prd}, n_O2: {n_O2_prd}, n_H2O: {n_H2O_cns}')
    ### permeation
    #n_H2O_eo, n_H2O_dd, n_H2O_pe =xflws.clc_crssflw_membrane_H2O(obj,obj.pec, i, A_cell)
    # H2O transport through membrane // in mol/(m² s)
    n_H2O_prm =xflws.clc_crssflw_membrane_H2O_chandesris(obj,obj.pec, T, i)#, A_cell)

    # Hydrogen permeation (ca -> an) // in mol/ (m² s)
    n_H2_prm = xflws.clc_hydrogen_permeation(obj, obj.pec, i)

    # Oxygen permeation // in mol/(m² s)
    n_O2_prm = xflws.clc_oxygen_permeation(obj, T, i)
    obj.av.n_O2_prm = n_O2_prm
    obj.av.n_H2_prm = n_H2_prm
    #print(f'n_ (flws, prm) -> n_H2: {n_H2_prm}, n_O2: {n_O2_prm}, n_H2O: {n_H2O_prm}')
    ### H2O
    # conc in channels
    c_H2O_ch = obj.av.rho_H2O / obj.pec.M_H2O
    #n_H2O_prm  = n_H2O_eo + n_H2O_dd + n_H2O_pe
    #print('n_H2O_cns: ', n_H2O_cns)

    #n_H2O_ch_in_an  =  (m_H2O_in_an / (obj.pec.M_H2O * obj.pcll.active_cell_area)) /stf # Convert massflow (kg/s) to molar flow (mol/m²s)
    n_H2O_ch_in_an  =  (m_H2O_in_an / (obj.pec.M_H2O)) # Convert massflow (kg/(m²s) to molar flow (mol/m²s)
    #print('n_H2O_in: ', n_H2O_ch_in_an)
    n_H2O_ch_out_an = n_H2O_ch_in_an - (n_H2O_prm + n_H2O_cns)
    n_H2O_ch_in_ca  = 0
    n_H2O_ch_out_ca = n_H2O_ch_in_ca + n_H2O_prm

    # H2

    #n_H2_ch_in_ca = 0

    n_H2_ch_out_ca    = n_H2_prd  - n_H2_prm if (n_H2_prm <= n_H2_prd) else 0#
    n_H2_ch_out_an    = n_H2_prm if (n_H2_prm < n_H2_prd) else n_H2_prd

    # O2

    # print(f'flws: n_O2_prd={n_O2_prd} | n_O2_prm={n_O2_prm}')
    # print(f'flws: n_H2_prd={n_H2_prd} |n_H2O_cns={n_H2O_cns}  | n_H2O_in={n_H2O_ch_in_an}')
    #n_H2_ch_in_ca = 0
    n_O2_ch_out_ca    = n_O2_prm if (n_O2_prm < n_O2_prd) else n_O2_prd#
    n_O2_ch_out_an    = n_O2_prd  - n_O2_prm if (n_O2_prm <= n_O2_prd) else 0

    # TODO  CHEck !!!!
    n_H2O_ca = (n_H2O_prm ) #+ n_H2O_cns)
    n_H2O_an = n_H2O_ch_out_an

    if False:
        plst = [n_H2_ch_out_an,n_H2_ch_out_ca,n_O2_ch_out_an,n_O2_ch_out_ca,
                n_H2O_an, n_H2O_ca]
        pnms = ['n_H2_ch_out_an','n_H2_ch_out_ca','n_O2_ch_out_an','n_O2_ch_out_ca',
                    'n_H2O_an', 'n_H2O_ca']
        for j,elm in enumerate(plst):
            if not np.isnan(plst[j]):
                print(pnms[j]+'('+str(j)+') :',plst[j])
    ### molar fractions
    '''
    calc. molar fractions of species in outlet-channels of cell
    Marangio_2009 eq.: 37
    '''
    #if i>0:
    # x_H2_ch_ca = n_H2_ch_out_ca / (n_H2_ch_out_ca + n_O2_ch_out_ca + n_H2O_ch_out_ca)
    # x_H2_ch_an = n_H2_ch_out_an / (n_H2_ch_out_an + n_O2_ch_out_an + n_H2O_ch_out_an)
    #
    # x_O2_ch_ca = n_O2_ch_out_ca / (n_H2_ch_out_ca + n_O2_ch_out_ca + n_H2O_ch_out_ca)
    # x_O2_ch_an = n_O2_ch_out_an / (n_H2_ch_out_an + n_O2_ch_out_an + n_H2O_ch_out_an)
    #
    # x_H2O_ch_ca = n_H2O_ch_out_ca / (n_H2_ch_out_ca + n_O2_ch_out_ca + n_H2O_ch_out_ca)
    # x_H2O_ch_an = n_H2O_ch_out_an / (n_H2_ch_out_an + n_O2_ch_out_an + n_H2O_ch_out_an)

    x_H2_ch_ca = xflws.division(n_H2_ch_out_ca, (n_H2_ch_out_ca + n_O2_ch_out_ca + n_H2O_ch_out_ca))
    x_H2_ch_an = xflws.division(n_H2_ch_out_an, (n_H2_ch_out_an + n_O2_ch_out_an + n_H2O_ch_out_an))

    x_O2_ch_ca = xflws.division(n_O2_ch_out_ca, (n_H2_ch_out_ca + n_O2_ch_out_ca + n_H2O_ch_out_ca))
    x_O2_ch_an = xflws.division(n_O2_ch_out_an, (n_H2_ch_out_an + n_O2_ch_out_an + n_H2O_ch_out_an))
    # print(f'n_H2O_ch_out_an={n_H2O_ch_out_an}')
    # print(f'x_O2_ch_an={x_O2_ch_an}')
    x_H2O_ch_ca = xflws.division(n_H2O_ch_out_ca,  (n_H2_ch_out_ca + n_O2_ch_out_ca + n_H2O_ch_out_ca))
    x_H2O_ch_an = xflws.division(n_H2O_ch_out_an,  (n_H2_ch_out_an + n_O2_ch_out_an + n_H2O_ch_out_an))

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
    print(f'x (ch_out): x_H2_ca= {x_H2_ch_ca} x_H2_an= {x_H2_ch_an} x_O2_ca= {x_O2_ch_ca} x_O2_an= {x_O2_ch_an}')
    # TODO: check n_H2O_an/ca !!!
    x_H2O_ch_ca = n_H2O_ca / (n_H2_ch_out_ca + n_O2_ch_out_ca + n_H2O_ca)
    x_H2O_ch_an = n_H2O_an / (n_H2_ch_out_an + n_O2_ch_out_an + n_H2O_an)
    '''
    obj.av.x_H2O_ch_ca = x_H2O_ch_ca
    obj.av.x_H2O_ch_an = x_H2O_ch_an

    ### concentrations
    '''
    clc concentration based on pressure (in compartment)
    and molar fraction

    Concentration of species in channels
    Marangio 2009 eq. 38
    '''
    #print('x_H2_ch_an: ',x_H2_ch_an)

    # -> in channels
    c_H2_ch_ca = p.cathode * x_H2_ch_ca / (obj.pec.R * T)
    c_H2_ch_an = p.anode * x_H2_ch_an / (obj.pec.R * T)

    c_O2_ch_ca = p.cathode * x_O2_ch_ca / (obj.pec.R * T)
    c_O2_ch_an = p.anode * x_O2_ch_an / (obj.pec.R * T)

    c_H2O_ch_an = c_H2O_ch_ca = obj.av.rho_H2O / obj.pec.M_H2O

    #print('c_H2_ch_ca: ',c_H2_ch_ca)
    #print('c_H2_ch_an: ',c_H2_ch_an)
    # -> at membrane/ electrode
    '''
    Concentration of species at membrane
    Marangio 2009 eq. 39
    '''
    # edit 20210618: remove A_cell -> area specific molar flows (in mol/ m² s)
    #Trinke d_cl = 10...20 micrometer
    #checked D_eff -> ok

    c_H2_mem_ca = (c_H2_ch_ca
                    + (obj.pec.d_el_ca *1e-4/ obj.av.D_eff_H2 ) * n_H2_prd)
    #print('(obj.pec.d_el_ca / obj.av.D_eff_H2 )=', (obj.pec.d_el_ca / obj.av.D_eff_H2 ))
    #TODO: check H2 at Anode
    c_H2_mem_an = (c_H2_ch_an
                    + (obj.pec.d_el_an *1e-4/ obj.av.D_eff_H2 ) * n_H2_ch_out_an)
    #TODO: check O2 at Cathode
    c_O2_mem_ca = (c_O2_ch_ca
                    + (obj.pec.d_el_ca *1e-4/ obj.av.D_eff_O2 ) * n_O2_ch_out_ca)
    c_O2_mem_an = (c_O2_ch_an
                    + (obj.pec.d_el_an *1e-4/ obj.av.D_eff_O2 ) * n_O2_prd)
    if c_O2_mem_an<0:
        print('c_O2_mem_an =', c_O2_mem_an)
        print('(obj.pec.d_el_an / obj.av.D_eff_O2 )=', (obj.pec.d_el_an / obj.av.D_eff_O2 ))
        c_O2_mem_an = 0

    if (c_H2_mem_an<0) and False:
        c_H2_mem_an = 0

    # print('c_H2_mem_ca: ',c_H2_mem_ca)
    # print('c_H2_mem_an: ',c_H2_mem_an)
    #c_H2O_mem = obj.av.rho_H2O / obj.pec.M_H2O
    # TODO: check D_eff for H2O
    # liquid water?
    c_H2O_mem_ca = (c_H2O_ch_ca
                    + (obj.pec.d_el_ca/obj.av.D_eff_H2 ) * n_H2O_ca ) # // in mol/m³
    c_H2O_mem_an = (c_H2O_ch_an
                    + (obj.pec.d_el_an/obj.av.D_eff_O2 ) * n_H2O_an) # Check !

    ### clc molar fractions at membrane
    #TODO: check, if correct !
    x_H2_mem_ca = c_H2_mem_ca * (obj.pec.R * T) / p.cathode
    x_H2_mem_an = c_H2_mem_an * (obj.pec.R * T) / p.anode

    x_O2_mem_ca = c_O2_mem_ca * (obj.pec.R * T) / p.cathode
    x_O2_mem_an = c_O2_mem_an * (obj.pec.R * T) / p.anode

    x_H2O_mem_ca = c_H2O_mem_ca / (c_H2O_mem_ca + c_H2_mem_ca)
    x_H2O_mem_an = c_H2O_mem_an / (c_H2O_mem_an + c_O2_mem_an)

    ### for plr clc
    obj.av.cnc_O2_mem_an = c_O2_mem_an
    obj.av.cnc_H2_mem_ca = c_H2_mem_ca
    ### clc partial pressures
    p.pp_H2_mem_ca = x_H2_mem_ca *p.cathode
    p.pp_H2_mem_an = x_H2_mem_an *p.anode
    p.pp_O2_mem_ca = x_O2_mem_ca *p.cathode
    p.pp_O2_mem_an = x_O2_mem_an *p.anode

    pp_H2O_mem_an = x_H2O_mem_an *p.anode

    if (sns==False) and (ntd!=None):

        ntd.n_H2_ca[m], ntd.n_H2_an[m] = n_H2_ch_out_ca*stf, n_H2_ch_out_an*stf
        ntd.n_O2_ca[m], ntd.n_O2_an[m] = n_O2_ch_out_ca*stf, n_O2_ch_out_an*stf
        ntd.c_H2_ca[m], ntd.c_H2_an[m] = c_H2_ch_ca, c_H2_ch_an
        ntd.c_O2_ca[m], ntd.c_O2_an[m] = c_O2_ch_ca, c_O2_ch_an
        ntd.x_H2_ca[m], ntd.x_H2_an[m] = x_H2_ch_ca, x_H2_ch_an
        ntd.x_O2_ca[m], ntd.x_O2_an[m] = x_O2_ch_ca, x_O2_ch_an
        ntd.pp_H2_ca[m], ntd.pp_H2_an[m] = p.pp_H2_mem_ca, p.pp_H2_mem_an
        ntd.pp_O2_ca[m], ntd.pp_O2_an[m] = p.pp_O2_mem_ca, p.pp_O2_mem_an
        #pp_H2O_mem_an
        ntd.n_H2O_cns[m] = n_H2O_cns*stf
        ntd.pp_H2O_an[m] = pp_H2O_mem_an
        if n_H2_ch_out_an>0:
            ntd.X_H2inO2[m] = n_H2_ch_out_an/(n_H2_ch_out_an+n_O2_ch_out_an)
        else:
            ntd.X_H2inO2[m] = 0
        return
    else:
        #TT = namedtuple('TT', ['test_a', 'test_aa','test_b', 'test_bb'])
        FLWS = namedtuple('FLWS', '''n_H2_out_ca n_H2_out_an n_O2_out_ca n_O2_out_an
                                    c_H2_out_ca c_H2_out_an c_O2_out_ca c_O2_out_an
                                    x_H2_out_ca x_H2_out_an x_O2_out_ca x_O2_out_an
                                    x_H2_mem_ca x_H2_mem_an x_O2_mem_ca x_O2_mem_an
                                    pp_H2_mem_ca pp_H2_mem_an pp_O2_mem_ca pp_O2_mem_an,
                                    pp_H2O_mem_an ''')

        flws_out = FLWS(n_H2_ch_out_ca, n_H2_ch_out_an, n_O2_ch_out_ca, n_O2_ch_out_an,
                        c_H2_ch_ca, c_H2_ch_an, c_O2_ch_ca, c_O2_ch_an,
                        x_H2_mem_ca, x_H2_mem_an, x_O2_mem_ca, x_O2_mem_an,
                        x_H2_ch_ca, x_H2_ch_an, x_O2_ch_ca, x_O2_ch_an,
                        p.pp_H2_mem_ca, p.pp_H2_mem_an, p.pp_O2_mem_ca, p.pp_O2_mem_an,
                        pp_H2O_mem_an )
        #tt = TT(1001,1,2002,2)
        return flws_out #tt#c_out, n_out

def partial_pressure(obj, pec, T, p):
    ''' partial pressure of product gases dependend on water-vapor-pressure'''
    #T in K
    #p_ca, p_an = p_in

    # ====== CEHCK below !
    #pp_H2O = xflws.clc_pp_H2O(obj, pec, T, )
    pp_H2O    = 1e5 *(6.1078 * 1e-3 * np.exp( 17.2694 * (T - 273.15) / (T - 34.85) )) # espinoza-lopez // in Pa

    pp_H2_ca = p.cathode - pp_H2O

    pp_O2_an = p.anode - pp_H2O

    return pp_H2_ca, pp_O2_an, pp_H2O #av.plr_ppr[0][1:3]
