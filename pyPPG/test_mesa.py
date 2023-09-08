import numpy as np

from pyPPG.example import*

DB_dir='C:/DataBase/MESA/EDF/'
IDs=[]
output_file=''
with open(DB_dir+'IDs.txt', "r") as file:
    # Iterate through each line in the file
    for line in file:
        # Remove leading and trailing whitespaces and add the ID to the list
        IDs.append(line.strip())


fid_names = ('on', 'sp', 'dn','dp', 'off', 'u', 'v', 'w', 'a', 'b', 'c', 'd', 'e', 'f', 'p1', 'p2')
SQI_mat=pd.DataFrame(columns=fid_names)

for id in IDs:
    print(id)
    data_path = DB_dir+id # the path of the file containing the PPG signal to be analysed
    start_sig = 1000000 # the first sample of the signal to be analysed
    end_sig = start_sig+256*30 # the last sample of the signal to be analysed
    savingfolder = 'temp_dir/MESA'
    output_file = savingfolder+'/SQI.csv'

    # run example code
    fiducials,s=ppg_example(data_path=data_path, start_sig=start_sig, end_sig=end_sig, process_type="fiducials", savingfolder=savingfolder, savefig=True)
    fp_keys=fiducials.keys()
    s_type = ['ppg', 'ppg', 'ppg','ppg', 'vpg', 'vpg', 'vpg', 'apg', 'apg', 'apg', 'apg', 'apg', 'apg', 'jpg','jpg']
    i=0

    for key in fp_keys:
        tmp_sig=s_type[0]

        # Add a single element to a specific cell
        row_index = id  # Index of the row where you want to add the element
        column_name = key  # Name of the column where you want to add the element
        try:
            sqi =  round(np.mean(SQI.get_ppgSQI(s.get_s()[tmp_sig].values, s.fs, fiducials[key])) * 100, 2)
        except:
            sqi = np.NAN

        # Update the DataFrame at the specified cell
        SQI_mat.at[row_index, column_name] = sqi

        SQI_mat.to_csv(output_file)
        i=i+1

    print(id)

print('Program finished')
