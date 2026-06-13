import numpy as np
import matplotlib.pyplot as plt
from color import ColorsData
from scipy.io import loadmat
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
import os
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


light_green = "#A8E6A3"
peach_orange     = "#FFCC99"   # 稍微暖一点



def apply_axis_style(ax):

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xticks(range(0, 21, 5))
    ax.set_xticks(np.arange(0, 8.4, 2))
    ax.tick_params(axis='both', which='major', labelsize=23)

    for spine in ax.spines.values():
        spine.set_linewidth(1)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))


def plot_error(GC_data, save_dir):

    GC_outsideL = GC_data['outsideL']
    GC_outsidebox = GC_data['outsidebox']
    t = GC_outsideL['t_span'].squeeze()

    
    Lpt        = GC_outsideL['Lpt']
    Lvt        = GC_outsideL['Lvt']
    L          = GC_outsideL['L']

    GCL_ep_list = GC_outsideL['ep_list']
    GCL_ev_list = GC_outsideL['ev_list']
    GCL_lyap_V     = GC_outsideL['lyap_list_V']

    GCB_ep_list = GC_outsidebox['ep_list']
    GCB_ev_list = GC_outsidebox['ev_list']  
    GCB_lyap_V     = GC_outsidebox['lyap_list_V']

    crop_index = len(t)
    legend_loc = (1., 1.0)
    colors = ColorsData().colors
    # separate figure to three plots
    # ----- LV -----
    f1 = plt.figure(figsize=(8,6), tight_layout=True)
    ax1 = f1.add_subplot(111)
    max_V = np.maximum(np.amax(GCL_lyap_V[:, :crop_index], axis=0),
                       np.amax(GCB_lyap_V[:, :crop_index], axis=0))
    # ax1.plot(t[:crop_index], L[:crop_index], color=red, linewidth=2.5)
    # ax1.plot(t[:crop_index], max_V, color='blue', linewidth=2.5)
    for i in range(GCL_lyap_V.shape[0]):
        ax1.plot(t[:crop_index], GCL_lyap_V[i, :crop_index],
                 color=light_green,  alpha=1, linewidth=1)
        ax1.plot(t[:crop_index], GCB_lyap_V[i, :crop_index],
                 color=peach_orange,  alpha=1, linewidth=1)
    ax1.set_ylabel('$V(t)$', fontsize=25)
    ax1.set_xlabel('$t$ [s]', fontsize=25)
    # ax1.set_xlim([0, t[crop_index-1]])
    ax1.set_xlim([0, 8+0.1])
    ax1.set_ylim([0, np.max(max_V) + 0.04])
    # ax1.legend([r'$\tilde{\mathcal{L}}^2$', r'$\max V_i$'], fontsize=18)
    ax1.plot([], [], color=light_green,  alpha=1, linewidth=1.0,
                label=r'$V_i$ (Exceeding the theoretical bound)')
    ax1.plot([], [], color=peach_orange,  alpha=1, linewidth=1.0,
                label=r'$V_i$ (Outside initial safe set)')
    ax1.plot(t[:crop_index], L[:crop_index], color=colors['L'], linewidth=3,
                label=r'$\tilde{\mathcal{L}}^2(\overline{\mathcal{V}}_{1}, \overline{\mathcal{V}}_{2}, t)$')
    apply_axis_style(ax1)
    ax1.legend(bbox_to_anchor=legend_loc, loc='upper right', fontsize=18)
    f1.savefig(os.path.join(save_dir, 'Fig_LV_unsatisfy.pdf'),format='pdf',bbox_inches='tight',pad_inches=0.1)
    f1.savefig(os.path.join(save_dir, 'Fig_LV_unsatisfy.eps'),format='eps',bbox_inches='tight',pad_inches=0.1)




    # ----- ep -----
    f2 = plt.figure(figsize=(8,6), tight_layout=True)
    ax2 = f2.add_subplot(111)
    max_ep = np.maximum(np.amax(GCL_ep_list[:, :crop_index], axis=0),
                       np.amax(GCB_ep_list[:, :crop_index], axis=0))
    for i in range(GCL_ep_list.shape[0]):
        ax2.plot(t[:crop_index], GCL_ep_list[i, :crop_index],
                 color=light_green,  alpha=1, linewidth=1)
        ax2.plot(t[:crop_index], GCB_ep_list[i, :crop_index],
                 color=peach_orange,  alpha=1, linewidth=1)
    # ax2.plot(t[:crop_index], Lpt[:crop_index], color=red, linewidth=3)
    ax2.set_ylabel(r'$\|e_p(t)\|$ [m]', fontsize=25)
    ax2.set_xlabel(r'$t$ [s]', fontsize=25)
    # ax2.set_xlim([0, t[crop_index-1]])
    ax2.set_xlim([0, 8+0.1])
    ax2.set_ylim([0, np.max(max_ep) + 0.04])
    ax2.plot([], [], color=light_green,  alpha=1, linewidth=1.0,
                label=r'$\|e_{i,p}(t)\|$ (Exceeding the theoretical bound)')
    ax2.plot([], [], color=peach_orange,  alpha=1, linewidth=1.0,
                label=r'$\|e_{i,p}(t)\|$ (Outside initial safe set)')

    ax2.plot(t[:crop_index], Lpt[:crop_index], color=colors['L'], linewidth=3,
                label=r'$\tilde{\mathcal{L}}_{p}(\overline{\mathcal{V}}_{1}, \overline{\mathcal{V}}_{2}, t)$')
    # ax2.legend([r'GC trajectories', r'$\tilde{\mathcal{L}}_p$'], fontsize=18)
    ax2.legend(bbox_to_anchor=legend_loc, loc='upper right', fontsize=18)
    apply_axis_style(ax2)
    f2.savefig(os.path.join(save_dir, 'Fig_ep_unsatisfy.pdf'),format='pdf',bbox_inches='tight',pad_inches=0.1)
    f2.savefig(os.path.join(save_dir, 'Fig_ep_unsatisfy.eps'),format='eps',bbox_inches='tight',pad_inches=0.1)

    # ----- ev -----
    f3 = plt.figure(figsize=(8,6), tight_layout=True)
    ax3 = f3.add_subplot(111)
    max_ev = np.maximum(np.amax(GCL_ev_list[:, :crop_index], axis=0),
                       np.amax(GCB_ev_list[:, :crop_index], axis=0))
    for i in range(GCL_ev_list.shape[0]):
        ax3.plot(t[:crop_index], GCL_ev_list[i, :crop_index],
                 color=light_green, alpha=1, linewidth=1)
        ax3.plot(t[:crop_index], GCB_ev_list[i, :crop_index],
                color=peach_orange, alpha=1, linewidth=1)
    ax3.plot(t[:crop_index], Lvt[:crop_index], color=colors['L'], linewidth=3)
    ax3.set_ylabel(r'$\|e_v(t)\|$ [m/s]', fontsize=25)
    ax3.set_xlabel(r'$t$ [s]', fontsize=25)
    # ax3.set_xlim([0, t[crop_index-1]])
    ax3.set_xlim([0, 8+0.1])
    ax3.set_ylim([0, np.max(max_ev) + 0.04])
    ax3.plot([], [], color=light_green,  alpha=1, linewidth=1.0,
                label=r'$\|e_{i,v}(t)\|$ (Exceeding the theoretical bound)')
    ax3.plot([], [], color=peach_orange,  alpha=1, linewidth=1.0,
                label=r'$\|e_{i,v}(t)\|$ (Outside initial safe set)')
    ax3.plot(t[:crop_index], Lvt[:crop_index], color=colors['L'], linewidth=3,
                label=r'$\tilde{\mathcal{L}}_{v}(\overline{\mathcal{V}}_{1}, \overline{\mathcal{V}}_{2}, t)$')
    # ax2.legend([r'GC trajectories', r'$\tilde{\mathcal{L}}_p$'], fontsize=18)
    ax3.legend(bbox_to_anchor=legend_loc, loc='upper right', fontsize=18)
    apply_axis_style(ax3)
    f3.savefig(os.path.join(save_dir, 'Fig_ev_unsatisfy.pdf'),format='pdf',bbox_inches='tight',pad_inches=0.1)
    f3.savefig(os.path.join(save_dir, 'Fig_ev_unsatisfy.eps'),format='eps',bbox_inches='tight',pad_inches=0.1)



def main():
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    save_dir = os.path.join(BASE_DIR, 'figures')
    os.makedirs(save_dir, exist_ok=True)
    
    GC_path = os.path.join(BASE_DIR, '..', 'GC', 'results_submit')
    GC_outsideL = loadmat(os.path.join(GC_path, 'plot_error_outside.mat'))
    Gc_outsidebox = loadmat(os.path.join(GC_path, 'plot_error_outsidebox.mat'))
    GC_date = {
        "outsideL": GC_outsideL,
        "outsidebox": Gc_outsidebox,
    }

    plot_error(GC_date, save_dir)
    plt.show()

if __name__ == "__main__":
    main()
