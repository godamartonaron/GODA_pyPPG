import matplotlib.pyplot as plt
import numpy as np

def plot_fig(bm_vals,ppgSQI,s, fp):
    tmp_keys=list(bm_vals['ppg_sig'].keys())
    figure_width = 1200
    figure_height = 700

    fig = plt.figure(figsize=(figure_width / 100, figure_height / 100))
    fig.subplots_adjust(hspace=.3, wspace=0)

    ax = plt.subplot(211)
    plt.plot(s.ppg,color='black')
    SQI_txt='Mean PPG SQI: '+str(ppgSQI)+ '%'

    plt.text(0.9, 0.95,SQI_txt , transform=plt.gca().transAxes, horizontalalignment='right',verticalalignment='top', fontweight='bold', fontsize=20)

    marker = ['o', 's', 's', 'o']
    color = ['r', 'b', 'g', 'm']
    fid_names = ('sp', 'on', 'dn', 'dp')
    for n in fid_names:
        tmp_pnt = eval("fp." + n + ".values")
        ind = fid_names.index(n)
        plt.scatter(tmp_pnt, s.ppg[tmp_pnt.astype('int')], s=60,linewidth=2, marker = marker[ind],facecolors='none', color=color[ind],label=n)

    plt.xlabel('Time [s]', fontsize=20)
    plt.ylabel('PPG', fontsize=20)

    str_sig = 0
    end_sig = len(s.ppg)
    len_sig = end_sig - str_sig
    step_small = 1
    step_big = step_small * 5

    major_ticks_names = range(0, int(len_sig / s.fs), step_big)
    len_ticks_names = len(major_ticks_names)
    major_diff = len_sig / len_ticks_names
    major_ticks = np.arange(str_sig, end_sig, major_diff)
    ax.set_xticks(major_ticks)
    plt.xticks(major_ticks, major_ticks_names, fontsize=20)
    plt.yticks([])

    leg = plt.legend(loc='upper left', fontsize=20, ncol=2, facecolor="white", frameon=False)
    for text in leg.get_texts():
        text.set_weight('bold')

    for i in range(0,6):
        plt.subplot(212)
        tmp_key=tmp_keys[i]
        plt.plot(bm_vals['ppg_sig'][tmp_key],label=tmp_key)
        leg = plt.legend(loc='upper right', fontsize=20, ncol=3, facecolor="white", frameon=False)
        for text in leg.get_texts():
            text.set_weight('bold')

        plt.xlabel('Index of pulse wave', fontsize=20)
        plt.xticks(fontsize=20)
        plt.ylabel('Time [s]', fontsize=20)
        plt.yticks(fontsize=20)

    plt.show()
