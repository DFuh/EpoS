'''
files for handling signal input
@dafu, 20200818
'''
import csv

def group_by_time(df):
    '''
    group signal-df by time(years)
    ---
    requires:
    df.columns = ['date', ... , ...]
    ---
    returns
    list of dicts with
    d { df: ,
        year: ,
        date_begin: ,
        date_end: ,}
    '''
    per_year = df.date.dt.to_period('Y')
    agg = df.groupby([per_year])

    lod = []
    for num,item in enumerate(agg):
        d = {}
        # full df of year
        #dfi[num] = item[1]
        d['df'] = item[1]
        # year
        #yri[num] = item[0].year
        d['year'] = item[0].year
        # startdate
        #mindate[num] = item[1].date.min()

        d['date_begin'] = item[0].start_time
        # enddate
        #maxdate[num]  = item[1].date.max()
        d['date_end'] = item[0].end_time    # CHECK! (output may cause trouble: Timestamp('2010-12-31 23:59:59.999999999'))

        #print(item[1].date.min())
        lod.append(d)
    return lod


def find_line(pth_to_file, search_key, search_key2=None, num_end=30):
    '''
        get line in csv, with search key
        or any specified search text
    '''
    with open(pth_to_file, 'r') as f:
        for num, line in enumerate(f):#,1):
            if search_key in line:
                return num
            if num > num_end:
                return None

def read_metadata(pth, ln):
    d = {}
    with open(pth, 'r') as f:
        rf =csv.reader(f)
        for num, line in enumerate(rf):
            print('num: ', num, 'ln: ', ln)
            if num < ln:
                if line[0].strip()[0].isalpha():
                    print('--line: ', line)
                    key, value = line[0], line[1]
                    d[key.strip()] = value.strip()
    return d

def read_data():
    return
