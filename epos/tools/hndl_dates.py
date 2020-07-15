import datetime
from collections import namedtuple

def todaysdate():
    now = datetime.datetime.now()
    NT = namedtuple('TODD', 'year month day hour minute')
    todd = NT(now.year, now.month, now.day, now.hour, now.minute)
    return todd
