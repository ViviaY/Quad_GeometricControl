import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
from scipy.io import loadmat
from color import ColorsData 
from matplotlib import gridspec
from matplotlib.patches import ConnectionPatch
import matplotlib.ticker as mtick
import matplotlib.patches as patches

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

import matplotlib.ticker as mtick   # 你已经 import 过就不用重复

def xfmt(x, pos):
    """set the format 0.00 """
    if abs(x) < 1e-9:     # avoid -0.00
        return "0.0"
    return f"{x:.2f}"

def plot_with_break(ax_bottom, ax_top, t, y, threshold, **kwargs):
    """Plot curve into bottom axis (normal region) and top axis (spike region).
    """
    # creat mask
    below_mask = y <= threshold
    above_mask = y > threshold
    
    # expand mask: if adjacent points are on the other side, include them to avoid broken lines
    below_extended = below_mask.copy()
    above_extended = above_mask.copy()
    
    # expand one point forward and backward
    for i in range(len(y)):
        if below_mask[i]:
            if i > 0:
                above_extended[i-1] = True
            if i < len(y) - 1:
                above_extended[i+1] = True
        if above_mask[i]:
            if i > 0:
                below_extended[i-1] = True
            if i < len(y) - 1:
                below_extended[i+1] = True
    
    ax_bottom.plot(t, np.where(below_extended, y, np.nan), **kwargs)
    ax_top.plot(t, np.where(above_extended, y, np.nan), **kwargs)




def draw_figlevel_zoom_indicator(fig,
                                 box_xy=(0.62, 0.43), 
                                 box_w=0.10,
                                 box_h=0.08,
                                 arrow_to=(0.78, 0.20)):



    # Draw dashed rectangle patch at figure level 
    rect = patches.Rectangle(
        box_xy,
        box_w,
        box_h,
        linewidth=1.2,
        linestyle="--",
        edgecolor='black',
        facecolor='none',
        alpha=0.8,
        transform=fig.transFigure,
        zorder=2000
    )
    fig.patches.append(rect)

    # Draw arrow using a hidden Axes
    ax_tmp = fig.add_axes([0, 0, 1, 1], frameon=False)
    ax_tmp.set_axis_off()     

    ax_tmp.annotate(
        "",
        xy=arrow_to,
        xytext=(box_xy[0] + box_w*0.85, box_xy[1]),
        xycoords=fig.transFigure,
        textcoords=fig.transFigure,
        arrowprops=dict(arrowstyle="->", ls="--", lw=1.1, alpha=0.8, color='black'),
        zorder=3000
    )

def split_legend_NMPC(ax_bottom, ax_top, keyword="NMPC"):

    fig = ax_bottom.figure   
    handles, labels = ax_bottom.get_legend_handles_labels()

    left_h, left_l = [], []   
    right_h, right_l = [], []   

    for h, lab in zip(handles, labels):
        if keyword in lab:
            right_h.append(h)
            right_l.append(lab)
        else:
            left_h.append(h)
            left_l.append(lab)

    # note: the position of legend depends on ax_top ---
    legend = fig.legend(
        left_h + right_h,
        left_l + right_l,
        ncol=2,
        columnspacing=1.5,
        labelspacing=0.3,
        fontsize=17,
        frameon=True,
        loc='upper right',
        bbox_to_anchor=(1, 1),         # position anchor uses top subplot
        bbox_transform=ax_top.transAxes
    )

    legend.set_zorder(9999)
    legend.get_frame().set_alpha(0.9)

    return legend



def extract_envelope(data, N=100):
    """
    trajectory_list: list of arrays, each array shape = (T, d)
                     d can be 3 for v, or 1 for f
    return: envelope curve shape (T,)
    """
    mats = []
    for traj in data:
        mats.append(np.max(np.abs(traj), axis=1))  # (T,)

    mats = np.stack(mats, axis=0)  # (N, T)
    env = np.max(mats, axis=0)     # (T,)
    return env



def save_subfigure_v_panel(t,
                           GC_v_max,
                           MPC_v_list,
                           CBF_v_list,
                           vmax_line=2.0, 
                           save_dir=''):
    
    fig = plt.figure(figsize=(8, 6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 5], hspace=0.05)
    ax_top = fig.add_subplot(gs[0])
    ax_bottom = fig.add_subplot(gs[1], sharex=ax_top)
    colors = ColorsData().colors
    ax_bottom.axhline(vmax_line, color=colors["L"],
                      linewidth=2.5, label=r'$|v_{\max}|$')

    threshold = vmax_line + 0.65

    plot_with_break(ax_bottom, ax_top, t, GC_v_max, threshold,
                    color=colors["GC_1"], linewidth=2.5, label='GC')

    H_list = [10, 15, 20, 25]
    NMPC_colors = ["NMPC10", "NMPC15", "NMPC20", "NMPC25"]
    for H, v_curve, ckey in zip(H_list, MPC_v_list, NMPC_colors):
        plot_with_break(ax_bottom, ax_top, t, v_curve, threshold,
                        color=colors[ckey], linestyle='--',
                        linewidth=2.5, label=rf'NMPC $(H={H})$')

    CBF_labels = [r'CBF $(\epsilon_{i,p})$',
                  r'CBF $(\epsilon_{safe})$',
                  r'CBF $(\epsilon_{reach})$']
    CBF_color_keys = ["CBF_epi", "CBF_safety", "CBF_reach"]
    for v_curve, lab, ckey in zip(CBF_v_list, CBF_labels, CBF_color_keys):
        plot_with_break(ax_bottom, ax_top, t, v_curve, threshold,
                        color=colors[ckey],
                        linestyle=(0, (6, 2, 2, 2)),
                        linewidth=2.5, label=lab)

    max_v_vals = [np.max(GC_v_max)] + [np.max(v) for v in MPC_v_list] + \
                 [np.max(v) for v in CBF_v_list]
    ylim_max = np.max(max_v_vals) + 0.02

    ax_top.spines['bottom'].set_visible(False)
    ax_bottom.spines['top'].set_visible(False)
    ax_top.tick_params(labelbottom=False, bottom=False)
    ax_top.xaxis.set_visible(False)

    d = 0.015
    kwargs_diag = dict(transform=ax_top.transAxes, color='k',
                       clip_on=False, linewidth=1)
    ax_top.plot((-d, +d), (-d, +d), **kwargs_diag)
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs_diag)
    kwargs_diag.update(transform=ax_bottom.transAxes)
    ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs_diag)
    ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs_diag)

    ax_top.set_ylim(threshold + 0.5, ylim_max)
    ax_bottom.set_ylim(0, threshold)
    ax_bottom.set_xlim(0, t[-1])

    ax_parent = fig.add_subplot(gs[:, :], frameon=False)
    ax_parent.tick_params(labelcolor='none', bottom=False, left=False)
    ax_parent.set_ylabel(r'$|v(t)|$ [m/s]', fontsize=26, labelpad=25)
    ax_bottom.set_xlabel('$t$ [s]', fontsize=26)

    split_legend_NMPC(ax_bottom, ax_top, keyword="NMPC")
    axins = ax_bottom.inset_axes([0.4, 0.453, 0.55, 0.28])
    axins.plot(t[:10], CBF_v_list[2][:10],   # reach
               color=colors["CBF_reach"],
               linestyle=(0, (6, 2, 2, 2)), linewidth=2.5)
    axins.plot(t[:10], CBF_v_list[1][:10],   # safety
               color=colors["CBF_safety"],
               linestyle=(0, (6, 2, 2, 2)), linewidth=2.5)
    axins.tick_params(axis='both', which='major', labelsize=12)
    axins.set_xlim([0, 0.035])
    axins.set_xticks([0, 0.01, 0.02, 0.03])

    def xfmt(x, pos):
        if abs(x) < 1e-9:
            return "0.0"
        return f"{x:.2f}"

    axins.xaxis.set_major_formatter(mtick.FuncFormatter(xfmt))
    axins.set_ylim([2, np.max(CBF_v_list[1][:10]) + 0.3])
    axins.set_yticks([2.0, 4.0, 6.0, 8.0, 10.0])
    axins.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))


    for ax in [ax_top, ax_bottom]:
        ax.set_xticks(range(0, 21, 5))
        ax.tick_params(axis='both', which='major', labelsize=22)
        for spine in ax.spines.values():
            spine.set_linewidth(1)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    fig.canvas.draw()
    draw_figlevel_zoom_indicator(
        fig,
        box_xy=(0.043, 0.832),  
        box_w=0.02,
        box_h=0.11,
        arrow_to=(0.36, 0.6),
    )


    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.12)
    fig.savefig(os.path.join(save_dir, 'Fig_v_sub.pdf'),format='pdf',bbox_inches='tight',pad_inches=0.1)
    fig.savefig(os.path.join(save_dir, 'Fig_v_sub.eps'),format='eps',bbox_inches='tight',pad_inches=0.1)



def save_subfigure_f_panel(t,
                           GC_f_max,
                           MPC_f_list,
                           CBF_f_list,
                           f_line=50.04,
                           save_dir=''):

    fig = plt.figure(figsize=(8, 6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 5], hspace=0.05)
    ax_top = fig.add_subplot(gs[0])
    ax_bottom = fig.add_subplot(gs[1], sharex=ax_top)
    colors = ColorsData().colors
    ax_bottom.axhline(f_line, color=colors["L"],
                      linewidth=2.5, label=r'$\overline{\mathcal{F}}$')

    threshold = f_line + 20

    plot_with_break(ax_bottom, ax_top, t, GC_f_max, threshold,
                    color=colors["GC_1"], linewidth=2.5, label='GC')

    H_list = [10, 15, 20, 25]
    NMPC_colors = ["NMPC10", "NMPC15", "NMPC20", "NMPC25"]
    for H, f_curve, ckey in zip(H_list, MPC_f_list, NMPC_colors):
        plot_with_break(ax_bottom, ax_top, t[:-1], f_curve, threshold,
                        color=colors[ckey], linestyle='--',
                        linewidth=2.5, label=rf'NMPC $(H={H})$')

    CBF_labels = [r'CBF $(\epsilon_{i,p})$',
                  r'CBF $(\epsilon_{safe})$',
                  r'CBF $(\epsilon_{reach})$']
    CBF_color_keys = ["CBF_epi", "CBF_safety", "CBF_reach"]
    for f_curve, lab, ckey in zip(CBF_f_list, CBF_labels, CBF_color_keys):
        plot_with_break(ax_bottom, ax_top, t, f_curve, threshold,
                        color=colors[ckey],
                        linestyle=(0, (6, 2, 2, 2)),
                        linewidth=2.5, label=lab)

    max_f_vals = [np.max(GC_f_max)] + [np.max(f) for f in MPC_f_list] + \
                 [np.max(f) for f in CBF_f_list]
    ylim_max = np.max(max_f_vals) + 0.02

    ax_top.spines['bottom'].set_visible(False)
    ax_bottom.spines['top'].set_visible(False)
    ax_top.tick_params(labelbottom=False, bottom=False)
    ax_top.xaxis.set_visible(False)

    d = 0.015
    kwargs_diag = dict(transform=ax_top.transAxes, color='k',
                       clip_on=False, linewidth=1)
    ax_top.plot((-d, +d), (-d, +d), **kwargs_diag)
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs_diag)
    kwargs_diag.update(transform=ax_bottom.transAxes)
    ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs_diag)
    ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs_diag)

    ax_top.set_ylim(threshold + 5, ylim_max)
    ax_bottom.set_ylim(0, threshold)
    ax_bottom.set_xlim(0, t[-1])

    ax_parent = fig.add_subplot(gs[:, :], frameon=False)
    ax_parent.tick_params(labelcolor='none', bottom=False, left=False)
    ax_parent.set_ylabel(r'$|f(t)|$ [N]', fontsize=26, labelpad=25)
    ax_bottom.set_xlabel('$t$ [s]', fontsize=26)

    split_legend_NMPC(ax_bottom, ax_top, keyword="NMPC")


    # ======== first zoom-in ========
    # axins = ax_top.inset_axes([0.12, 0.25, 0.55, 0.28])
    axins = fig.add_axes([0.13, 0.66, 0.23, 0.20])
    axins.plot(t[:5], CBF_f_list[2][:5],
               color=colors["CBF_reach"],
               linestyle=(0, (6, 2, 2, 2)), linewidth=2.5)
    axins.plot(t[:5], CBF_f_list[1][:5],
               color=colors["CBF_safety"],
               linestyle=(0, (6, 2, 2, 2)), linewidth=2.5)

    for H, f_curve, ckey in zip(H_list, MPC_f_list, NMPC_colors):
        axins.plot(t[:5], f_curve[:5],
                    color=colors[ckey],
                    linestyle="--", linewidth=2.5)

    axins.tick_params(axis='both', which='major', labelsize=12)
    axins.set_xlim([0, t[4]])
    axins.set_xbound(0, t[4])
    axins.autoscale(enable=False)
    # axins.set_xticks([0, 0.01, 0.02, 0.03])
    axins.xaxis.set_major_formatter(mtick.FuncFormatter(xfmt))
    yticks = np.arange(80, ylim_max, 500)
    axins.set_yticks(yticks)
    axins.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    for ax in [ax_top, ax_bottom]:
        ax.set_xticks(range(0, 21, 5))
        ax.tick_params(axis='both', which='major', labelsize=22)
        for spine in ax.spines.values():
            spine.set_linewidth(1)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    draw_figlevel_zoom_indicator(
        fig,
        box_xy=(0.043, 0.832),  
        box_w=0.02,
        box_h=0.11,
        arrow_to=(0.0938, 0.77),
    )

    # ======== second zoom-in =========
    # axins2 = ax_bottom.inset_axes([0.15, 0.55, 0.55, 0.28])  
    axins2 = fig.add_axes([0.35, 0.28, 0.45, 0.20])

    for H, f_curve, ckey in zip(H_list, MPC_f_list, NMPC_colors):
        axins2.plot(t[:400], f_curve[:400],
                    color=colors[ckey],
                    linestyle="--", linewidth=2.5)

    axins2.plot(t[:400], GC_f_max[:400], color=colors["GC_1"], linewidth=2.5)
    axins2.axhline(f_line, color=colors["L"], linewidth=2.5)

    for f_curve, lab, ckey in zip(CBF_f_list, CBF_labels, CBF_color_keys):
        axins2.plot(t[:400], f_curve[:400],
                    color=colors[ckey],
                    linestyle=(0, (6, 2, 2, 2)), linewidth=2.5)


    axins2.tick_params(axis='both', which='major', labelsize=12)
    axins2.set_xlim([0, t[400]])
    axins2.set_xbound(0, t[400])
    axins2.autoscale(enable=False)
    # axins2.set_xticks([t[2], t[4], t[6], t[8], t[10]])
    yticks = np.arange(41.5, 46.5, 2)
    axins2.set_ybound(41.5, 46.5)
    axins2.set_yticks(yticks)
    axins2.set_ylim(41.5, 46.5)
    axins2.xaxis.set_major_formatter(mtick.FuncFormatter(xfmt))
    axins2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    draw_figlevel_zoom_indicator(
        fig,
        box_xy=(0.043, 0.52),  
        box_w=0.1,
        box_h=0.06,
        arrow_to=(0.3056, 0.4375),
    )
    fig.canvas.draw()

    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.12)
    fig.savefig(os.path.join(save_dir, 'Fig_f_sub.pdf'),format='pdf',bbox_inches='tight',pad_inches=0.1)
    fig.savefig(os.path.join(save_dir, 'Fig_f_sub.eps'),format='eps',bbox_inches='tight',pad_inches=0.1)



def save_subfigure_Fd_panel(t, GC_Fd3, save_dir=''):
    fig, ax = plt.subplots(figsize=(8, 6))

    colordata = ColorsData()
    color_d = colordata.colors_v1

    for i in range(GC_Fd3.shape[0]):
        fi = GC_Fd3[i, :]
        ax.plot(t, fi, color=color_d[i % len(color_d)], linewidth=1.5)

    ax.set_ylabel(r'$F_{d,3}(t)$ [N]', fontsize=26)
    ax.set_xlim([0, t[-1]])
    ax.set_xticks(range(0, 21, 5))
    ax.tick_params(axis='both', which='major', labelsize=22)
    for spine in ax.spines.values():
        spine.set_linewidth(1)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlabel('$t$ [s]', fontsize=26)

    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.12)
    fig.savefig(os.path.join(save_dir, 'Fig_Fd3_sub.pdf'),format='pdf',bbox_inches='tight',pad_inches=0.1)
    fig.savefig(os.path.join(save_dir, 'Fig_Fd3_sub.eps'),format='eps',bbox_inches='tight',pad_inches=0.1)


def plot_vf(t, GC_data, NMPC_data, CBF_data, save_dir):
    crop = 1
    NMPC_H10 = NMPC_data["H10"]
    MPC_H10_v = extract_envelope([NMPC_H10[f"track_v_{i}"]for i in range(1, 101)])
    MPC_H10_f = extract_envelope([NMPC_H10[f"track_u_{i}"][:, [0]] for i in range(1, 101)])

    NMPC_H15 = NMPC_data["H15"]
    MPC_H15_v =  extract_envelope([NMPC_H15[f"track_v_{i}"]for i in range(1, 101)])
    MPC_H15_f = extract_envelope([NMPC_H15[f"track_u_{i}"][:, [0]] for i in range(1, 101)])
    
    NMPC_H20 = NMPC_data["H20"]
    MPC_H20_v =  extract_envelope([NMPC_H20[f"track_v_{i}"]for i in range(1, 101)])
    MPC_H20_f = extract_envelope([NMPC_H20[f"track_u_{i}"][:, [0]] for i in range(1, 101)])

    NMPC_H25 = NMPC_data["H25"]
    MPC_H25_v =  extract_envelope([NMPC_H25[f"track_v_{i}"]for i in range(1, 101)])
    MPC_H25_f = extract_envelope([NMPC_H25[f"track_u_{i}"][:, [0]] for i in range(1, 101)])

    CBF_epi = CBF_data["epi"]
    CBF_epi_v =  extract_envelope([CBF_epi[f"track_v_{i}"]for i in range(1, 101)])
    CBF_epi_f = extract_envelope([CBF_epi[f"track_u_{i}"][:, [0]] for i in range(1, 101)])


    CBF_safety = CBF_data["safety"]
    CBF_safety_v =  extract_envelope([CBF_safety[f"track_v_{i}"]for i in range(1, 101)])
    CBF_safety_f = extract_envelope([CBF_safety[f"track_u_{i}"][:, [0]] for i in range(1, 101)])

    CBF_reach = CBF_data["reach"]
    CBF_reach_v =  extract_envelope([CBF_reach[f"track_v_{i}"]for i in range(1, 101)])
    CBF_reach_f = extract_envelope([CBF_reach[f"track_u_{i}"][:, [0]] for i in range(1, 101)])

    GC_v_list = GC_data["GC_v"]["v_list"]
    GC_v_max = extract_envelope([GC_v_list[i].T for i in range(GC_v_list.shape[0])])
    GC_f = GC_data["GC_u"]["f_list"]
    GC_Fd3 = GC_data["GC_u"]["Fd3_list"]
    GC_f_max   = np.max(np.abs(GC_f), axis=0)


    
    save_subfigure_v_panel(
        t,
        GC_v_max,
        [MPC_H10_v, MPC_H15_v, MPC_H20_v, MPC_H25_v],
        [CBF_epi_v, CBF_safety_v, CBF_reach_v],
        vmax_line=2.0,
        save_dir=save_dir
    )

    save_subfigure_f_panel(
        t,
        GC_f_max,
        [MPC_H10_f, MPC_H15_f, MPC_H20_f, MPC_H25_f],
        [CBF_epi_f, CBF_safety_f, CBF_reach_f],
        f_line=50.04,
        save_dir=save_dir
    )

    save_subfigure_Fd_panel(t, GC_Fd3, save_dir=save_dir)
    plt.show()


def main():
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    save_dir = os.path.join(BASE_DIR, 'figures')
    os.makedirs(save_dir, exist_ok=True)

    GC_path = os.path.join(BASE_DIR, '..', 'GC', 'results_submit')
    GC_v = loadmat(os.path.join(GC_path, 'plot_v_satisfy.mat'))
    GC_u = loadmat(os.path.join(GC_path, 'plot_fFd3_satisfy.mat'))

    GC_data = {
        "GC_v": GC_v,
        "GC_u": GC_u,
    }
    data =  loadmat(os.path.join(GC_path, 'plot_error_satisfy.mat')) 
    t = data['t_span'].squeeze()
    initial_n = data['initial_n'][0][0]

    NMPC_path = os.path.join(BASE_DIR, '..', 'NMPC', 'results_submit')
    NMPC_H10_data = loadmat(os.path.join(NMPC_path, 'NMPC_trajs_100_H10.mat'))
    NMPC_H15_data = loadmat(os.path.join(NMPC_path, 'NMPC_trajs_100_H15.mat'))
    NMPC_H20_data = loadmat(os.path.join(NMPC_path, 'NMPC_trajs_100_H20.mat'))
    NMPC_H25_data = loadmat(os.path.join(NMPC_path, 'NMPC_trajs_100_H25.mat'))
    MPC_data = {
        "H10": NMPC_H10_data,
        "H15": NMPC_H15_data,
        "H20": NMPC_H20_data,
        "H25": NMPC_H25_data,
    }

    CBF_path = os.path.join(BASE_DIR, '..', 'CBFs', 'results_submit')
    CBF_epi_data = loadmat(os.path.join(CBF_path, 'CBFs_tracking_100_initial.mat'))
    CBF_safety_data = loadmat(os.path.join(CBF_path, 'CBFs_tracking_100_safety.mat'))
    CBF_reaching_data = loadmat(os.path.join(CBF_path, 'CBFs_tracking_100_reach.mat'))
    CBF_data = {
        "epi": CBF_epi_data,
        "safety": CBF_safety_data,
        "reach": CBF_reaching_data,
    }

    plot_vf(t, GC_data, MPC_data, CBF_data, save_dir)



if __name__ == "__main__":
    main()
