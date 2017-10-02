from fileios import *
msg = 'Enter Absolute Path to file: '
f_name = raw_input(msg).strip()

path = file_data_and_path(f_name)
if path != None:
       print 'Path:',path
