from FiducialPoints import*
from Biomarkers import*
from Statistics import*
from Summary import*

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap
from tkinter import filedialog
import mne
import time

from six.moves import cPickle as pickle

import matplotlib.mlab
from scipy.io import savemat

###########################################################################
####################### Data Acquisition from Files #######################
###########################################################################
if __name__ == '__main__':

    ppg_sig_dir='D:/ALL_DATA/Uni/Subjects/ITK_Adjunktus/HAIFA/TECHNION-BME/Research/PPG/PETE_Matlab'
    ppg_file = '/PPG-BP1.mat'
    annot_path = ppg_sig_dir+'/ANNOTS/MG_PPG-BP_annot/merged'
    sig_path=(ppg_sig_dir + ppg_file)
    input_sig = scipy.io.loadmat(sig_path)

    dn, u, v, a, b, c, d, e = [], [], [], [], [], [], [], []
    dnr, ur, vr, ar, br, cr, dr, er = [], [], [], [], [], [], [], []
    dn_dist, u_dist,v_dist,a_dist,b_dist,c_dist,d_dist,e_dist = [],[],[],[],[],[],[],[]
    dist_error = np.empty([218,8])
    dist_error [:] = np.NaN

    for i in range(0,218):
        hr = input_sig['ppg_data']['filt_sig'][0,i]
        hr =np.squeeze(hr)

        fs = input_sig['ppg_data']['fs'][0,0][0][0]
        fs = np.squeeze(fs)


        drt1 = input_sig['ppg_data']['d1'][0,i]
        drt1 = np.squeeze(drt1)

        drt2 = input_sig['ppg_data']['d2'][0,i]
        drt2 = np.squeeze(drt2)

        name=input_sig['ppg_data']['name'][0, i][0]
        annot_file=annot_path+'/'+name+'.mat'
        annot=scipy.io.loadmat(annot_file)
        peaks_loc = [np.squeeze(np.round(annot['annot']['pk'][0, 0][0][0]['t'] * fs).astype(int))]
        onsets_loc = np.squeeze(np.round(annot['annot']['os'][0, 0][0][0]['t']*fs).astype(int))

        ref_dn = np.squeeze(np.round(annot['annot']['dn'][0, 0][0][0]['t'] * fs).astype(int))
        ref_u = np.squeeze(np.round(annot['annot']['u'][0, 0][0][0]['t'] * fs).astype(int))
        ref_v = np.squeeze(np.round(annot['annot']['v'][0, 0][0][0]['t'] * fs).astype(int))
        ref_a = np.squeeze(np.round(annot['annot']['a'][0, 0][0][0]['t'] * fs).astype(int))
        ref_b = np.squeeze(np.round(annot['annot']['b'][0, 0][0][0]['t'] * fs).astype(int))
        ref_c = np.squeeze(np.round(annot['annot']['c'][0, 0][0][0]['t'] * fs).astype(int))
        ref_d = np.squeeze(np.round(annot['annot']['d'][0, 0][0][0]['t'] * fs).astype(int))
        ref_e = np.squeeze(np.round(annot['annot']['e'][0, 0][0][0]['t'] * fs).astype(int))

        det_dn=getDicroticNotch(hr, fs, peaks_loc, onsets_loc)
        det_u, det_v = getFirstDerivitivePoints(hr, fs, onsets_loc)
        det_a, det_b, det_c, det_d, det_e = getSecondDerivitivePoints(hr,fs, onsets_loc)

        dnr.append(np.squeeze(ref_dn))
        ur.append(np.squeeze(ref_u))
        vr.append(np.squeeze(ref_v))
        ar.append(np.squeeze(ref_a))
        br.append(np.squeeze(ref_b))
        cr.append(np.squeeze(ref_c))
        dr.append(np.squeeze(ref_d))
        er.append(np.squeeze(ref_e))

        dn.append(np.squeeze(det_dn))
        u.append(np.squeeze(det_u))
        v.append(np.squeeze(det_v))
        a.append(np.squeeze(det_a))
        b.append(np.squeeze(det_b))
        c.append(np.squeeze(det_c))
        d.append(np.squeeze(det_d))
        e.append(np.squeeze(det_e))

        dn_dist.append(np.squeeze(ref_dn - det_dn))
        u_dist.append(np.squeeze(ref_u - det_u))
        v_dist.append(np.squeeze(ref_v - det_v))
        a_dist.append(np.squeeze(ref_a - det_a))
        b_dist.append(np.squeeze(ref_b - det_b))
        c_dist.append(np.squeeze(ref_c - det_c))
        d_dist.append(np.squeeze(ref_d - det_d))
        e_dist.append(np.squeeze(ref_e - det_e))

        temp_dist = np.squeeze([ref_dn - det_dn, ref_u - det_u, ref_v - det_v, ref_a - det_a,
                        ref_b - det_b, ref_c - det_c, ref_d - det_d,ref_e - det_e])
        for j in range(0,8):
            if temp_dist[j].size>0:
                dist_error[i][j] = np.squeeze(temp_dist[j])

        print(str(i),' Abs errors of',name,': dn=',dn_dist[i],': u=',u_dist[i],' | v=',v_dist[i],' | a=',a_dist[i],
             ' | b=',b_dist[i],' | c=',v_dist[i],' | d=',d_dist[i],' | e=',e_dist[i])

        plt.plot(hr,'r',label='x')
        plt.plot(drt1,'b',label='dx')
        plt.plot(drt2,'k',label='ddx')

        plt.plot(onsets_loc, hr[onsets_loc], 'ks', label='onset')
        plt.plot(ref_dn, hr[ref_dn], 'ko', label='ref dn')
        plt.plot(ref_u, drt1[ref_u], 'go', label='ref u')
        plt.plot(ref_v, drt1[ref_v], 'mo', label='ref v')
        plt.plot(ref_a, drt2[ref_a], 'bs', label='ref a')
        plt.plot(ref_b, drt2[ref_b], 'cs', label='ref b')
        plt.plot(ref_c, drt2[ref_c], 'ks', label='ref c')
        plt.plot(ref_d, drt2[ref_d], 'bs', label='ref d')
        plt.plot(ref_e, drt2[ref_e], 'gs', label='ref e')

        plt.plot(det_dn, hr[det_dn], 'm*', label='det dn')
        plt.plot(det_u, drt1[det_u], 'r*', label='det u')
        plt.plot(det_v, drt1[det_v], 'b*', label='det v')
        plt.plot(det_a, drt2[det_a], 'mx', label='det a')
        plt.plot(det_b, drt2[det_b], 'kx', label='det b')
        plt.plot(det_c, drt2[det_c], 'rx', label='det c')
        plt.plot(det_d, drt2[det_d], 'cx', label='det d')
        plt.plot(det_e, drt2[det_e], 'r+', label='det e')

        plt.legend(loc=4, prop={'size': 10})
        plt.title(name,  fontsize=20)
        plt.xlabel('Time [ms]',  fontsize=20)
        plt.ylabel('Pulse Wave',  fontsize=20)
        plt.grid(color = 'green', linestyle = '--', linewidth = 0.5)
        # plt.show()
        plt.close('all')

        file_name = 'temp_dir/PPG-BP1_eval02.mat'
        mdic = {"dist_error": dist_error}
        savemat(file_name, mdic)

    MAE_dn = np.round(np.mean(np.absolute(dn_dist)),2)
    MAE_u = np.round(np.mean(np.absolute(u_dist)),2)
    MAE_v = np.round(np.mean(np.absolute(v_dist)),2)
    MAE_a = np.round(np.mean(np.absolute(a_dist)),2)
    MAE_b = np.round(np.mean(np.absolute(b_dist)),2)
    MAE_c = np.round(np.mean(np.absolute(c_dist)),2)
    MAE_d = np.round(np.mean(np.absolute(d_dist)),2)
    MAE_e = np.round(np.mean(np.absolute(e_dist)),2)

    std_dn = np.round(np.std(dn_dist),2)
    std_u = np.round(np.std(u_dist),2)
    std_v = np.round(np.std(v_dist),2)
    std_a = np.round(np.std(a_dist),2)
    std_b = np.round(np.std(b_dist),2)
    std_c = np.round(np.std(c_dist),2)
    std_d = np.round(np.std(d_dist),2)
    std_e = np.round(np.std(e_dist),2)

    bias_dn = np.round(np.mean(dn_dist),2)
    bias_u = np.round(np.mean(u_dist),2)
    bias_v = np.round(np.mean(v_dist),2)
    bias_a = np.round(np.mean(a_dist),2)
    bias_b = np.round(np.mean(b_dist),2)
    bias_c = np.round(np.mean(c_dist),2)
    bias_d = np.round(np.mean(d_dist),2)
    bias_e = np.round(np.mean(e_dist),2)

    print('-------------------------------------------')
    print('Bias',': dn=', bias_dn,' | u=', bias_u, ' | v=', bias_v, ' | a=', bias_a,
          ' | b=', bias_b, ' | c=', bias_c, ' | d=', bias_d, ' | e=', bias_e)
    print('MAE',': dn=', MAE_dn,' |  u=', MAE_u, ' | v=', MAE_v, ' | a=', MAE_a,
          ' | b=', MAE_b, ' | c=', MAE_c, ' | d=', MAE_d, ' | e=', MAE_e)

    print('STD',': dn=', std_dn,' |  u=', std_u, ' | v=', std_v, ' | a=', std_a,
          ' | b=', std_b, ' | c=', std_c, ' | d=', std_d, ' | e=', std_e)


