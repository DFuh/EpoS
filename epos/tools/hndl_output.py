'''
'''


def ini_data_output(self):
        '''

        - make output-directory
        - prepare csv-files (one per year)

        Returns
        -------
        None.

        '''

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
