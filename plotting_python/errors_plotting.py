from matplotlib import patches
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as mtick
import matplotlib.patches as patches
from color import ColorsData 
import os

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


def plot_with_break(ax_bottom, ax_top, t, y, threshold, **kwargs):
    """Plot curve into bottom axis (normal region) and top axis (spike region).
    """
 
    below_mask = y <= threshold
    above_mask = y > threshold
    
    below_extended = below_mask.copy()
    above_extended = above_mask.copy()

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


def save_subfigure_V(t, L, lyap_V, save_dir=''):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    colors = ColorsData().colors

    ax.plot(t, L, color=colors["L"], linewidth=2.5,
            label=r'$\tilde{\mathcal{L}}^2(\overline{\mathcal{V}}_{1}, \overline{\mathcal{V}}_{2}, t)$')

    max_V = np.max(lyap_V, axis=0)
    ax.plot(t, max_V, color=colors["GC_1"], linewidth=2.5,
            label=r'$\max_{i\in[1,100]} V_i(t)$')

    ax.set_ylabel('$V(t)$', fontsize=26)
    ax.set_xlabel('$t$ [s]', fontsize=26)

    ax.set_xlim([0, t[-1]])
    ax.set_ylim([0, np.max(L) + 0.02])

    ax.legend(loc="upper right", fontsize=17)

    ax.tick_params(axis='both', which='major', labelsize=22)
    ax.set_xticks(range(0, 21, 5))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'Fig_V.pdf'),format='pdf',bbox_inches='tight',pad_inches=0.1)
    fig.savefig(os.path.join(save_dir, 'Fig_V.eps'),format='eps',bbox_inches='tight',pad_inches=0.1)

    
    # plt.close(fig)

def save_subfigure_ep(t, GC_ep_list, NMPC_ep_max_list, CBF_ep_max_list, Lpt, save_dir=''):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    colors = ColorsData().colors
    # Lp(t)
    ax.plot(t, Lpt, color=colors["L"], linewidth=2.5,
            label=r'$\tilde{\mathcal{L}}_{p}(\overline{\mathcal{V}}_{1},\overline{\mathcal{V}}_{2},t)$')

    # GC max
    max_ep_gc = np.max(GC_ep_list, axis=0)
    ax.plot(t, max_ep_gc, color=colors["GC_1"], linewidth=2.5, label=r'GC')

    # NMPC
    H_list = [10, 15, 20, 25]
    color_keys = ["NMPC10", "NMPC15", "NMPC20", "NMPC25"]
    for ep, H, ck in zip(NMPC_ep_max_list, H_list, color_keys):
        ax.plot(t, ep, color=colors[ck], linestyle="--", linewidth=2.5,
                label=rf'NMPC $(H={H})$')

    # CBF
    cbf_labels = [r'CBF $(\epsilon_{i,p})$', r'CBF $(\epsilon_{safe})$', r'CBF $(\epsilon_{reach})$']
    cbf_colors = ["CBF_epi", "CBF_safety", "CBF_reach"]
    for ep, lbl, ck in zip(CBF_ep_max_list, cbf_labels, cbf_colors):
        ax.plot(t, ep, color=colors[ck], linestyle=(0, (6, 2, 2, 2)), linewidth=2.5, label=lbl)

    ax.set_ylabel(r'$\|e_p(t)\|$ [m]', fontsize=26)
    ax.set_xlabel(r'$t$ [s]', fontsize=26)

    ax.set_xlim([0, t[-1]])
    ax.set_ylim([0, 1.05 * max(
        [np.max(Lpt), np.max(max_ep_gc)] +
        [np.max(ep) for ep in NMPC_ep_max_list] +
        [np.max(ep) for ep in CBF_ep_max_list]
    )])
    ax.set_xticks(range(0, 21, 5))

    ax.tick_params(axis='both', labelsize=22)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    split_legend_NMPC(ax)

    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'Fig_ep.pdf'),format='pdf',bbox_inches='tight',pad_inches=0.1)
    fig.savefig(os.path.join(save_dir, 'Fig_ep.eps'),format='eps',bbox_inches='tight',pad_inches=0.1)


    # plt.close(fig)

def save_subfigure_ev(t, GC_ev_list, NMPC_ev_max, CBF_ev_max, Lvt, save_dir=''):

    threshold = np.max(Lvt) + 0.02

    fig = plt.figure(figsize=(8, 6))
    gs2 = fig.add_gridspec(2, 1, height_ratios=[1, 5], hspace=0.05)
    ax2_top    = fig.add_subplot(gs2[0])
    ax2_bottom = fig.add_subplot(gs2[1], sharex=ax2_top)

    ax2_bottom.margins(x=0)
    ax2_top.margins(x=0)
    colors = ColorsData().colors

    plot_with_break(
        ax2_bottom, ax2_top, t, Lvt, threshold,
        color=colors["L"], linewidth=2.5,
        label=r'$\tilde{\mathcal{L}}_{v}(\overline{\mathcal{V}}_{1}, \overline{\mathcal{V}}_{2}, t)$'
    )

    # GC
    max_ev_GC = np.amax(GC_ev_list, axis=0)
    plot_with_break(
        ax2_bottom, ax2_top, t, max_ev_GC,threshold,
        color=colors["GC_1"], linewidth=2.5, label='GC'
    )

    # NMPC 
    H_list = [10, 15, 20, 25]
    NMPC_colors = [colors["NMPC10"], colors["NMPC15"],
                   colors["NMPC20"], colors["NMPC25"]]
    for ev_max, H, c in zip(NMPC_ev_max, H_list, NMPC_colors):
        plot_with_break(
            ax2_bottom, ax2_top, t, ev_max, threshold,
            color=c, linestyle='--', linewidth=2.5,
            label=rf'NMPC $(H={H})$'
        )

    # CBF 
    CBF_labels = [r'CBF $(\epsilon_{i,p})$',
                  r'CBF $(\epsilon_{safe})$',
                  r'CBF $(\epsilon_{reach})$']
    CBF_colors = [colors["CBF_epi"], colors["CBF_safety"], colors["CBF_reach"]]
    for ev_max, lab, c in zip(CBF_ev_max, CBF_labels, CBF_colors):
        plot_with_break(
            ax2_bottom, ax2_top, t, ev_max, threshold,
            color=c, linestyle=(0, (6, 2, 2, 2)), linewidth=2.5,
            label=lab
        )

    all_max = [
        np.max(Lvt), np.max(max_ev_GC),
        *[np.max(ev) for ev in NMPC_ev_max],
        *[np.max(ev) for ev in CBF_ev_max]
    ]
    ylim_max = np.max(all_max) + 0.02

    ax2_bottom.set_ylim(0, threshold)
    ax2_top.set_ylim(2, ylim_max)  

    ax2_top.spines['bottom'].set_visible(False)
    ax2_bottom.spines['top'].set_visible(False)
    ax2_top.tick_params(labelbottom=False, bottom=False)
    ax2_top.xaxis.set_visible(False)

    d = 0.015
    kwargs_diag = dict(color='k', clip_on=False, linewidth=1)
    ax2_top.plot((-d, +d), (-d, +d), transform=ax2_top.transAxes, **kwargs_diag)
    ax2_top.plot((1 - d, 1 + d), (-d, +d), transform=ax2_top.transAxes, **kwargs_diag)
    ax2_bottom.plot((-d, +d), (1 - d, 1 + d), transform=ax2_bottom.transAxes, **kwargs_diag)
    ax2_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), transform=ax2_bottom.transAxes, **kwargs_diag)

    ax2_bottom.set_xlim(0, t[-1] + 0.1)
    ax2_bottom.set_xlabel('$t$ [s]', fontsize=26)

    ax2_parent = fig.add_subplot(gs2[:, :], frameon=False)
    ax2_parent.tick_params(labelcolor='none', bottom=False, left=False)
    ax2_parent.set_ylabel(r'$\|e_v(t)\|$ [m/s]', fontsize=26, labelpad=25)

    axins = ax2_bottom.inset_axes([0.31, 0.34, 0.55, 0.28])
    axins.plot(t[:10], CBF_ev_max[2][:10],  # reach
               color=colors["CBF_reach"], linestyle=(0, (6, 2, 2, 2)), linewidth=2.5)
    axins.plot(t[:10], CBF_ev_max[1][:10],  # safety
               color=colors["CBF_safety"], linestyle=(0, (6, 2, 2, 2)), linewidth=2.5)

    axins.tick_params(axis='both', which='major', labelsize=12)
    axins.tick_params(axis='both', which='minor', labelsize=10)
    axins.set_xlim([0, 0.035])
    axins.set_xticks([0, 0.01, 0.02, 0.03])

    def xfmt(x, pos):
        if abs(x) < 1e-9:
            return "0.00"
        return f"{x:.2f}"

    axins.xaxis.set_major_formatter(mtick.FuncFormatter(xfmt))
    axins.set_ylim([2, np.max(CBF_ev_max[1][:10]) + 0.3])
    axins.set_yticks([2.0, 4.0, 6.0, 8.0, 10.0])
    axins.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    for ax in [ax2_bottom, ax2_top]:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xticks(range(0, 21, 5))
        ax.tick_params(axis='both', which='major', labelsize=22)
        for spine in ax.spines.values():
            spine.set_linewidth(1)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    split_legend_NMPC(ax2_top, keyword="NMPC")

    fig.canvas.draw()
    draw_figlevel_zoom_indicator(
        fig,
        box_xy=(0.11, 0.832), 
        box_w=0.02,
        box_h=0.11,
        arrow_to=(0.3256, 0.53),
    )

    plt.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.12)
    fig.savefig(os.path.join(save_dir, 'Fig_ev.pdf'),format='pdf',bbox_inches='tight',pad_inches=0.1)
    fig.savefig(os.path.join(save_dir, 'Fig_ev.eps'),format='eps',bbox_inches='tight',pad_inches=0.1)



def draw_figlevel_zoom_indicator(fig,
                                 box_xy=(0.62, 0.43), 
                                 box_w=0.10,
                                 box_h=0.08,
                                 arrow_to=(0.78, 0.20)):



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

    # ===== 2) 使用一个隐藏的 Axes 来画箭头 =====
    ax_tmp = fig.add_axes([0, 0, 1, 1], frameon=False)
    ax_tmp.set_axis_off()      # 不显示任何坐标轴

    ax_tmp.annotate(
        "",
        xy=arrow_to,
        xytext=(box_xy[0] + box_w*0.85, box_xy[1]),
        xycoords=fig.transFigure,
        textcoords=fig.transFigure,
        arrowprops=dict(arrowstyle="->", ls="--", lw=1.1, alpha=0.8, color='black'),
        zorder=3000
    )

def split_legend_NMPC(ax, keyword="NMPC"):
    """
    左列：Lp/Lv 上界 + GC + CBFs
    右列：所有 NMPC(H=*) 系列
    Legend 属于 figure（非 axes），永不被遮挡。
    """

    fig = ax.figure                       # <<< 关键：从 ax 获取 figure

    handles, labels = ax.get_legend_handles_labels()

    left_h, left_l = [], []    # 上界 + GC + CBF
    right_h, right_l = [], []  # 所有 NMPC

    for h, lab in zip(handles, labels):
        if keyword in lab:
            right_h.append(h)
            right_l.append(lab)
        else:
            left_h.append(h)
            left_l.append(lab)

    # Figure-level legend（最上层）
    legend = fig.legend(
        left_h + right_h,
        left_l + right_l,
        ncol=2,
        columnspacing=1.5,
        labelspacing=0.3,
        fontsize=17,
        frameon=True,
        loc='upper right',               # 可以改为 'upper right'
        bbox_to_anchor=(1, 1),      #  ← 推荐位置，不被遮挡
        bbox_transform=ax.transAxes
    )

    #     ax.legend(
    #     left_h + right_h,
    #     left_l + right_l,
    #     ncol=2,
    #     columnspacing=1.5,
    #     labelspacing=0.3,
    #     fontsize=17,
    #     frameon=True,
    #     loc='upper right'
    # )


    legend.set_zorder(9999)               # <<< 永远置顶图层（关键）
    legend.get_frame().set_alpha(0.9)    # 可选：半透明背景

    return legend


def plot_ev(NMPC_data):
    
    NMPC_H10 = NMPC_data["H10"]
    # MPC_H10_ep = [NMPC_H10[f"traj_{i}_ep"]for i in range(1, 101)]
    MPC_H10_ep = np.stack([NMPC_H10[f"traj_{i}_ep"] for i in range(1, 101)], axis=0)
    MPC_H10_ep_norm = np.linalg.norm(MPC_H10_ep, axis=2)
    NMPC_H10_ep_max =NMPC_H10["max_ep_nmpc"].squeeze()

    NMPC_H15 = NMPC_data["H15"]
    # MPC_H15_ep = [NMPC_H15[f"traj_{i}_ep"]for i in range(1, 101)]
    MPC_H15_ep = np.stack([NMPC_H15[f"traj_{i}_ep"] for i in range(1, 101)], axis=0)
    MPC_H15_ep_norm = np.linalg.norm(MPC_H15_ep, axis=2)
    NMPC_H15_ep_max =NMPC_H15["max_ep_nmpc"].squeeze()

    NMPC_H20 = NMPC_data["H20"]
    # MPC_H20_ep = [NMPC_H20[f"traj_{i}_ep"]for i in range(1, 101)]
    MPC_H20_ep = np.stack([NMPC_H20[f"traj_{i}_ep"] for i in range(1, 101)], axis=0)
    MPC_H20_ep_norm = np.linalg.norm(MPC_H20_ep, axis=2)
    NMPC_H20_ep_max =NMPC_H20["max_ep_nmpc"].squeeze()

    NMPC_H25 = NMPC_data["H25"]
    # MPC_H25_ep = [NMPC_H25[f"traj_{i}_ep"]for i in range(1, 101)]
    MPC_H25_ep = np.stack([NMPC_H25[f"traj_{i}_ep"] for i in range(1, 101)], axis=0)
    MPC_H25_ep_norm = np.linalg.norm(MPC_H25_ep, axis=2)
    NMPC_H25_ep_max =NMPC_H25["max_ep_nmpc"].squeeze()


    t = NMPC_H10["traj_1_t"].squeeze()
    ep_norm = [MPC_H10_ep_norm, MPC_H15_ep_norm, MPC_H20_ep_norm, MPC_H25_ep_norm]
    ep_max = [NMPC_H10_ep_max, NMPC_H15_ep_max, NMPC_H20_ep_max, NMPC_H25_ep_max]
    fig, axs = plt.subplots(2, 2, figsize=(16, 12), tight_layout=True)
    for ax, ep, ep_m, H in zip(axs.flatten(), ep_norm, ep_max, [10, 15, 20, 25]):
        ax.plot(t, ep_m, color='r', linewidth=2.5, linestyle='-',label=rf'$\mathrm{{MPC}}~\max_{{i \in [1,100]}} \|e_p(t)\|$ (H={H})')
        for ep_norm_i in ep:
            ax.plot(t, ep_norm_i, color='gray', alpha=0.8, linewidth=1.5)

        ax.grid(True)
        ax.legend(loc='upper right', fontsize=17)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis='both', which='major', labelsize=22)
        for spine in ax.spines.values():
            spine.set_linewidth(1)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    

    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.12)
    plt.savefig('./MPC_ep_smoothnesscheck.pdf', format='pdf', bbox_inches='tight', pad_inches=0.1)
    plt.show()
    plt.show()





def plot_error(t, GC_data, NMPC_data, CBF_data, save_dir):

    NMPC_H10_ep_max = NMPC_data["H10"]["max_ep_nmpc"].squeeze()
    NMPC_H10_ev_max = NMPC_data["H10"]["max_ev_nmpc"].squeeze()

    NMPC_H15_ep_max = NMPC_data["H15"]["max_ep_nmpc"].squeeze()
    NMPC_H15_ev_max = NMPC_data["H15"]["max_ev_nmpc"].squeeze()

    NMPC_H20_ep_max = NMPC_data["H20"]["max_ep_nmpc"].squeeze()
    NMPC_H20_ev_max = NMPC_data["H20"]["max_ev_nmpc"].squeeze()

    NMPC_H25_ep_max = NMPC_data["H25"]["max_ep_nmpc"].squeeze()
    NMPC_H25_ev_max = NMPC_data["H25"]["max_ev_nmpc"].squeeze()


    CBF_epi_ep_max = CBF_data["epi"]["max_ep_cbf"].squeeze()
    CBF_epi_ev_max = CBF_data["epi"]["max_ev_cbf"].squeeze()
    CBF_safety_ep_max = CBF_data["safety"]["max_ep_cbf"].squeeze()
    CBF_safety_ev_max = CBF_data["safety"]["max_ev_cbf"].squeeze()
    CBF_reach_ep_max = CBF_data["reach"]["max_ep_cbf"].squeeze()
    CBF_reach_ev_max = CBF_data["reach"]["max_ev_cbf"].squeeze()


    GC_ep_list = GC_data['ep_list']
    GC_ev_list = GC_data['ev_list']
    lyap_V = GC_data['lyap_list_V']
    Lpt = GC_data['Lpt']
    Lvt = GC_data['Lvt']
    L = GC_data['L']

    save_subfigure_V(t, L, lyap_V, save_dir=save_dir)
    save_subfigure_ep(t, GC_ep_list, 
                    [NMPC_H10_ep_max, NMPC_H15_ep_max, NMPC_H20_ep_max, NMPC_H25_ep_max],
                    [CBF_epi_ep_max, CBF_safety_ep_max, CBF_reach_ep_max], Lpt, save_dir=save_dir)
    save_subfigure_ev(t, GC_ev_list,
                    [NMPC_H10_ev_max, NMPC_H15_ev_max, NMPC_H20_ev_max, NMPC_H25_ev_max],
                    [CBF_epi_ev_max, CBF_safety_ev_max, CBF_reach_ev_max], Lvt, save_dir=save_dir)

    plt.show()




def main():

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    save_dir = os.path.join(BASE_DIR, 'figures')
    os.makedirs(save_dir, exist_ok=True)

    GC_path = os.path.join(BASE_DIR, '..', 'GC', 'results_submit')
    GC_data = loadmat(os.path.join(GC_path, 'plot_error_satisfy.mat'))
    t = GC_data['t_span'].squeeze()

    NMPC_path = os.path.join(BASE_DIR, '..', 'NMPC', 'results_submit')
    NMPC_data = {
        "H10": loadmat(os.path.join(NMPC_path, 'NMPC_error_100_H10.mat')),
        "H15": loadmat(os.path.join(NMPC_path, 'NMPC_error_100_H15.mat')),
        "H20": loadmat(os.path.join(NMPC_path, 'NMPC_error_100_H20.mat')),
        "H25": loadmat(os.path.join(NMPC_path, 'NMPC_error_100_H25.mat')),
    }


    CBF_path = os.path.join(BASE_DIR, '..', 'CBFs', 'results_submit')
    CBF_data = {
        "epi": loadmat(os.path.join(CBF_path, "CBFs_error_100_initial.mat")),
        "safety": loadmat(os.path.join(CBF_path, "CBFs_error_100_safety.mat")),
        "reach": loadmat(os.path.join(CBF_path, "CBFs_error_100_reach.mat")),
    }

    plot_error(t, GC_data, NMPC_data, CBF_data, save_dir)



if __name__ == "__main__":
    main()
