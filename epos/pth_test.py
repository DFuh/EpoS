import sys
import os



def t_pth():
    home_pth = os.path.expanduser('~')
    curr_dir = os.path.dirname(__file__)
    print('home_pth: ', home_pth)
    print('curr_dir: ', curr_dir)
    return

if __name__ == '__main__':
    t_pth()
