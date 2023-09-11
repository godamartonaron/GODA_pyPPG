from pyPPG.example import*
from dotmap import DotMap
import matplotlib.pyplot as plt

# Define the file path
DB_dir='D:/ALL_DATA/Uni/Subjects/ITK_Adjunktus/HAIFA/TECHNION-BME/Research/Australia/Data/'
rec_name='003-20Mar2023 PPGpreintervention .txt'
file_path = DB_dir+rec_name # Replace 'your_file.txt' with the actual file path

# Create a dictionary to store the header information
header_data = {}

# Create a list to store the data lines
data_lines = []

# Flag to indicate when to start capturing data lines
capture_data = False

# Open the file and read its contents line by line
with open(file_path, 'r') as file:
    for line in file:
        # Strip leading and trailing whitespace
        line = line.strip()

        # Check if the line contains an equal sign (=) to split into key and value
        if '=' in line:
            key, value = line.split('=')
            # Store the key-value pair in the dictionary
            header_data[key.strip()] = value.strip()

            # Check if the line contains the marker to start capturing data
            if key.strip() == 'Range':
                capture_data = True
        elif capture_data:
            # If we're capturing data and the line doesn't contain '=', add it to data_lines
            data_lines.append(line)

# # Print the header data
# for key, value in header_data.items():
#     print(f"{key}: {value}")
#
# # Print the data lines
# print("\nData Lines:")
# for line in data_lines:
#     print(line)

channel_titles = header_data.get('ChannelTitle', '').split()

# Create dictionaries to store data for each channel
channel_data = {channel: [] for channel in channel_titles}

# Get the index of each channel in the data based on the header
channel_indices = {channel: channel_titles.index(channel) for channel in channel_titles}

# Iterate through the data lines and split them into values for each channel
for line in data_lines:
    values = line.split()
    for channel in channel_titles:
        channel_index = channel_indices[channel]
        channel_data[channel].append(float(values[channel_index]))

# # Create a dictionary with renamed keys
# old_keys=channel_data.keys()
# new_keys=['Time','ECG','PPG_f','PPG_e']
#
# # Define a mapping of old keys to new keys
# for i in range(0,len(old_keys)):
#     exec ('key_mapping = {"'+old_keys[i]+':"new_keys[i]}')
#
# # Create a new dictionary with renamed keys
# new_dict = {key_mapping.get(key, key): value for key, value in channel_data.items()}


# # Print the data for each channel
# for channel, data in channel_data.items():
#     print(f"{channel} Data:")
#     print(data)

# DB_dir='C:/DataBase/MESA/EDF/'
# IDs=[]
# output_file=''
# with open(DB_dir+'IDs.txt', "r") as file:
#     # Iterate through each line in the file
#     for line in file:
#         # Remove leading and trailing whitespaces and add the ID to the list
#         IDs.append(line.strip())


fid_names = ('on', 'sp', 'dn','dp', 'off', 'u', 'v', 'w', 'a', 'b', 'c', 'd', 'e', 'f', 'p1', 'p2')
SQI_mat=pd.DataFrame(columns=fid_names)

start_sig = 0 # the first sample of the signal to be analysed
end_sig = -1 # the last sample of the signal to be analysed
savingfolder = 'temp_dir/MESA'
output_file = savingfolder+'/SQI.csv'


sig=channel_data['Pulse']
signal = DotMap()
signal.fs=1000

signal.start_sig = start_sig
if start_sig < end_sig:
    signal.end_sig = end_sig
else:
    signal.end_sig = len(sig)

try:
    signal.v = sig[signal.start_sig:signal.end_sig]
except:
    raise ('There is no valid PPG signal!')

signal.name = rec_name
filtering=True

## Preprocessing
signal.ppg, signal.vpg, signal.apg, signal.jpg = Preprocessing(signal, filtering=filtering)

## Create a PPG class
signal.filtering = filtering
signal.correct = True
s = PPG(signal)

# Init the fiducials package
fpex = FP.FpCollection(s)
fiducials = fpex.get_fiducials(s, signal.correct)

fp = Fiducials(fiducials)
ppgSQI = round(np.mean(SQI.get_ppgSQI(s.ppg, s.fs, fp.sp)) * 100, 2)
print('Mean PPG SQI: ', ppgSQI, '%')

# Plot fiducial points
plot_fiducials(s, fp, savingfolder, show_fig=True)

print('Program finished')
