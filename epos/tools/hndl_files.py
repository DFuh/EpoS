'''
File handling
'''
import os
import datetime
import pandas as pd


class handleFiles():

    def __init__(self,):
        pass
    '''
    class InputF():
    '''
    '''
        handling/ arranging of input files
    '''
    '''

        def __init():
            pass


    class OutputF():

        def __init__(self,):
            pass

        def ini_output_files():

            for file in fllst:
                pass

            return
    '''

    def ini_output_file(self,):
        '''
        initialize output file:
        --- metadata ---
        specs of simu
        --- end metadata ---

        --- data ---
        header of data
        data of simu

        '''

        return

    def mk_output_dir(full_path, now=None):
        '''
        make new output directory (if not existing)
        if now is not none -> create subdir with todays date
        '''

        if now:
            full_path += '/'+now
        if not os.path.exists(full_path):
            print('+++ make new dir: ', full_path)
            os.makedirs(full_path)
        return



    def ini_logfile(simu_obj):
        '''
        create logfile for specific simulation

        ++ redundant code lines!
        # TODO: compact/ simplify
        '''
        s_pre = 'test1'
        s_suf = 'test2'
        filename = simu_obj.log_filename
        print(' ... creating logfile ... -> ', filename )
        log_now = str(datetime.datetime.now())
        specs_dict =    {
                        'name': 'logfile',
                        'ini_date': log_now,
                        'time_duration_simu': 0,


                        }# specs of logfile

        log_dict = simu_obj.__dict__.copy() #specs of data
        del log_dict['tec_parameters'] # avoid duplicating tec_parameters
        s_par_dict = simu_obj.s_parameters._asdict()
        #s_par_df = pd.DataFrame(simu_obj.s_parameters) #superior simu parameters
        #par_df = pd.DataFrame(simu_obj.tec_parameters)# tec-specific parameters
        par_dict = simu_obj.tec_parameters._asdict()
        with open(simu_obj.path_data_out + '/eposLog_' + filename + '.txt', 'w') as f:  # Just use 'w' mode in 3.x
            head_s = f'--- {s_pre}_log_specs_{s_suf} ---\n'
            l_head = len(head_s)
            f.write(l_head*'-'+'\n')
            f.write(head_s)
            [f.write(f'{key}'+(25-len(key))*' '+f'{value}\n') for key, value in specs_dict.items()]
            f.write('--- end specs ---\n\n')

            head_s = f'--- {s_pre}_log_simuspecs_{s_suf} ---\n'
            l_head = len(head_s)
            f.write(l_head*'-'+'\n')
            f.write(head_s)
            [f.write(f'{key}'+(25-len(key))*' '+f'{value}\n') for key, value in log_dict.items()]

            #TODO: sig-specs

            f.write('--- end dataspecs ---\n\n')

            head_s = f'--- {s_pre}_log_s_par_{s_suf} ---\n'
            l_head = len(head_s)
            f.write(l_head*'-'+'\n')
            f.write(head_s)
            #s_par_df.to_csv(f, index=False)
            [f.write(f'{key}'+(25-len(key))*' '+f'{value}\n') for key, value in s_par_dict.items()]
            f.write('--- end s_par\n\n')

            head_s = f'--- {s_pre}_log_tec_par_{s_suf} ---\n'
            l_head = len(head_s)
            f.write(l_head*'-'+'\n')
            f.write(head_s)
            #par_df.to_csv(f, index=False)
            #[f.write(f'{key}'+(20-len(key))*' '+f'{value}\n') for key, value in par_dict.items()]
            for key,value in par_dict.items():
                if hasattr(value, 'items'):
                    f.write(f'{key}\n')#+(20-len(key)*' ')#+f'{value}\n')
                    for subkey, subvalue in value.items():
                        f.write(5*' '+f'{subkey}'+(30-len(subkey))*' '+f'{subvalue}\n')
                else:
                    f.write(key+(30-len(key))*' '+f'{value}\n')
            f.write('--- end tec_par')
        return


    def pth_mirror(filepath='', bsc_pth='', ref_dir_path='/data/in', tday='', fl_prfx='', mk_newdir=False, out_pth=None):
        '''
        input:
        sig_filepath    | string | full input filepath
        bsc_pth         | string | ?
        ref_dir_pth     | string | reference path (subpath used for mirroring)

        create directory according to input file-location

        sig_filepath: string | relative path from
        bsc_pth: string | (rel) path of data

        #TODO: Code structure // redundant variable assignment
        '''
        #
        # -> sig_filepath

        # full basic filepath for input data # default: ./data/in
        #abs_pth = os.path.abspath(bsc_pth)
        #print('abs_pth: ', abs_pth)

        # split filepath
        filepath_head, filename = os.path.split(filepath)
        #print('pathmirror: filepath: ', filepath)
        #print('pathmirror: filepath: ', filepath_head)
        if (not ref_dir_path) or ('data/in' in ref_dir_path):
            subpath = filepath_head.split('/data/in')[1]
        else:
            subpath = filepath_head.split(ref_dir_path)[1] +'/ext'

        out_path = os.path.abspath('./data/out'+subpath+'/'+tday)


        '''
        if bsc_dir is None:
            bsc_dir = os.path.basename(bsc_pth) # basic input directory

        #out_pth_lst = []
        #for pth in in_pth_lst: # pth list: full pth+filename.suffix
        print('pth_in: ', pth_in)
        pure_pth, flnm = os.path.split(pth_in) #full pth
        last_dir = os.path.basename(pure_pth)
        #mat_pth0 = '/mat'
        pth_0 = bsc_pth
        dir_append = ''
        ### --->>>> use os.path.relpath(pth1, pth2) !!!
        #https://stackoverflow.com/questions/7287996/python-get-relative-path-from-comparing-two-absolute-paths#7288019
        '''
        '''
        k=0
        while (last_dir != bsc_dir) and k<20:
            print('last_dir: ', last_dir)
            print('bsc_dir: ',bsc_dir)
            dir_append = last_dir + '/' + dir_append
            pure_pth = os.path.split(pure_pth)[0] #??
            last_dir = os.path.basename(pure_pth)
            k+=1
        '''


        #add_pth = '/' + dir_append

        #dirnm = self.basepath
        #<<<< need relative path to file >>>>>
        #newpath = os.path.split(abs_pth)[0]+'/out'+ add_pth
        #full_filepath = newpath + fl_prfx + flnm
        #out_pth_lst.append(full_filepath)
        '''
        if out_pth is not None:
            newpath = os.path.abspath(out_pth)
            outpath = newpath
            #full_filepath = newpath + fl_prfx + flnm
        '''
        '''
        elif mk_newdir:
            print('+++ make new dir: ', newpath)
            if not os.path.exists(newpath):
                os.makedirs(newpath)
            outpath = newpath
            #full_filepath = newpath + fl_prfx + flnm
        ###########################################
        if now:
            dat = str(now.year)+str(now.month)+str(now.day)
            add_pth = '/' + dir_append + dat
            new_spath = newpath + dat
            if mk_newdir:
                print('+++ make new subdir: ', new_spath)
                if not os.path.exists(new_spath):
                    os.makedirs(new_spath)
            outpath = new_spath
        full_filepath = outpath + fl_prfx + flnm
        '''
        return out_path#, full_filepath
