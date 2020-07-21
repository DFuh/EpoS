'''
File handling
'''
import os



class handleFiles():

    def __init__():
        pass

    class InputF():
        '''
        handling/ arranging of input files
        '''

        def __init():
            pass


    class OutputF():

        def __init__():
            pass

        def ini_output_files():

            for file in fllst:
                pass

            return

    def pth_mirror(pth_in, bsc_pth='data/in', bsc_dir=None, now=None, fl_prfx='', mk_newdir=False, out_pth=None):
        '''
        create directory according to input file-location

        pth_in: string | relative path from

        #TODO: Code structure // redundant variable assignment
        '''

        abs_pth = os.path.abspath(bsc_pth)
        print('abs_pth: ', abs_pth)
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
        k=0
        while (last_dir != bsc_dir) and k<20:
            print('last_dir: ', last_dir)
            print('bsc_dir: ',bsc_dir)
            dir_append = last_dir + '/' + dir_append
            pure_pth = os.path.split(pure_pth)[0]
            last_dir = os.path.basename(pure_pth)
            k+=1


        add_pth = '/' + dir_append

        #dirnm = self.basepath
        #<<<< need relative path to file >>>>>
        newpath = os.path.split(abs_pth)[0]+'/out'+ add_pth
        full_filepath = newpath + fl_prfx + flnm
        #out_pth_lst.append(full_filepath)
        if out_pth is not None:
            newpath = os.path.abspath(out_pth)
            outpath = newpath
        elif mk_newdir:
            print('+++ make new dir: ', newpath)
            if not os.path.exists(newpath):
                os.makedirs(newpath)
            outpath = newpath
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

        return outpath, full_filepath
