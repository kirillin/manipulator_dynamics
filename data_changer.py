
PATH = 'data_for_identification_Youbot'
# PATH = 'data_for_identification_2dof'

DATA_FILE_NAME = PATH + '/data/data_{:d}{:s}.txt'
BIG_TAUs_FILE_NAME = PATH + '/bigs/{:s}big_tau{:d}.txt'
BIG_XIs_FILE_NAME = PATH + '/bigs/{:s}big_xi{:d}.txt'
EE_FILE_NAME = PATH + '/ee/EE{:d}.txt'

if __name__ == '__main__':
    for i in range(0,1):
        file_input = open(DATA_FILE_NAME.format(i, ''))
        file_output = open()
