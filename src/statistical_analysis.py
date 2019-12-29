from __future__ import print_function
import math
import matplotlib.pyplot as plt
from scipy.stats import *
import image_descriptors as idsc
from src.utils import enable_logger, cell_type_micropatterned, plot_colors, check_dir
from helpers import *

logger = logging.getLogger('DYPFISH_HELPERS')
np.set_printoptions(precision=2, suppress=True, linewidth=512)

# Implements equation 17 (supplemental)of padovan-merhar et al. 2015
def compare_volume_corrected_nm(file_handler, sec_file_handler, acquisition1, acquisition2):

    arhgdia = build_image_list(file_handler, 'mrna', acquisition1)
    arhgdia_cultured = build_image_list(file_handler, 'mrna', acquisition2)
    cell_area_arhgdia_cultured = [idsc.get_cell_area(sec_file_handler, image) for x in arhgdia_cultured]

    # build heightmap for cultured data before computing cell volume
    cell_volume_arhgdia=[idsc.compute_cell_volume(file_handler,image) for image in arhgdia]
    transcount_arhgdia = [len(idsc.get_spots(file_handler, image)) for image in arhgdia]
    transcount_arhgdia_cultured = [len(idsc.get_spots(file_handler, image)) for image in arhgdia_cultured]
    nms=[]

    for volumes, mrnas in [cell_volume_arhgdia,transcount_arhgdia],[cell_area_arhgdia_cultured,transcount_arhgdia_cultured]:
        coeffs = np.polyfit(volumes, mrnas, 1)
        a = coeffs[1]
        b = coeffs[0]
        expected_mrnas = [(a + (b * volume)) for volume in volumes]
        variance_expected_mrnas = np.std(expected_mrnas) ** 2
        variance_mrnas = np.std(mrnas) ** 2
        exp_mrnas = (np.mean(mrnas) ** 2)
        nm = (variance_mrnas - variance_expected_mrnas) / exp_mrnas
        print(variance_mrnas,variance_expected_mrnas,exp_mrnas,nm)
        nms.append(nm)
    return nms


# Implements equation 17 (supplemental)of padovan-merhar et al. 2015
def compute_volume_corrected_nm(file_handler, image_list):
    cell_volume=[idsc.compute_cell_volume(file_handler,image) for image in image_list]
    transcount = [len(idsc.get_spots(file_handler, image)) for image in image_list]
    coeffs = np.polyfit(cell_volume, transcount, 1)
    a = coeffs[1]
    b = coeffs[0]
    expected_mrnas = [(a + (b * volume)) for volume in cell_volume]
    variance_expected_mrnas = np.std(expected_mrnas) ** 2
    variance_mrnas = np.std(transcount) ** 2
    exp_mrnas = (np.mean(transcount) ** 2)
    nm = (variance_mrnas - variance_expected_mrnas) / exp_mrnas
    return nm

# Implements equation 17 (supplemental)of padovan-merhar et al. 2015
def compute_surface_corrected_nm(file_handler, image_list):
    cell_surface= [idsc.get_cell_area(file_handler, image) for image in image_list]
    transcount = [len(idsc.get_spots(file_handler, image)) for image in image_list]
    coeffs = np.polyfit(cell_surface, transcount, 1)
    a = coeffs[1]
    b = coeffs[0]
    expected_mrnas = [(a + (b * volume)) for volume in cell_surface]
    variance_expected_mrnas = np.std(expected_mrnas) ** 2
    variance_mrnas = np.std(transcount) ** 2
    exp_mrnas = (np.mean(transcount) ** 2)
    nm = (variance_mrnas - variance_expected_mrnas) / exp_mrnas
    return nm

def compute_mrna_peripheral_fraction(file_handler, image):
    spots_peripheral_flags = idsc.get_spots_peripheral_flags(file_handler, image)
    peripheral_fraction = spots_peripheral_flags.sum() / len(spots_peripheral_flags)
    return peripheral_fraction

def compute_protein_peripheral_fraction(file_handler, image):
    cell_mask=idsc.get_cell_mask(file_handler,image)
    nucleus_mask=idsc.get_nucleus_mask(file_handler,image)
    IF=idsc.get_IF(file_handler,image)
    peripheral_mask=idsc.get_peripheral_mask(file_handler,image)
    cytoplasm_mask = (cell_mask == 1) & (nucleus_mask == 0)
    peripheral_fraction = IF[peripheral_mask[:, :]].sum() / IF[cytoplasm_mask[:, :]].sum()
    return peripheral_fraction

def boxplot_profile(profiles, gene, slice_number, figname):
    logger.info('Build boxplot peripheral profile for gene ' + gene)
    plt.boxplot(np.matrix(profiles), 'gD')
    plt.ylabel('mRNA enrichment')
    plt.xlabel('Periph:cyt division, ' + str(int(np.around(100 / slice_number))) + '%')
    plt.title(gene + ' mrna enrichment profile')
    plt.axis([0, slice_number + 1, 0, np.around(np.amax(profiles)) + 0.1])
    plt.savefig(figname)
    plt.show()
    plt.close()

def compare_nucleus_area(file_handler, sec_file_handler, acquisition1, acquisition2):
    arhgdia = build_image_list(file_handler, 'mrna', acquisition1)
    arhgdia_cultured = build_image_list(file_handler, 'mrna', acquisition2)
    nucleus_area_arhgdia = [idsc.get_nucleus_area(sec_file_handler, x) for x in arhgdia]
    nucleus_area_arhgdia_cultured = [idsc.get_nucleus_area(sec_file_handler, x) for x in arhgdia_cultured]
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    ax.tick_params(right=False,top=False,direction='inout', length=8, width=3, colors='black')
    bp = ax.boxplot([nucleus_area_arhgdia, nucleus_area_arhgdia_cultured])
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3)
    plt.setp(bp['boxes'], color='black',linewidth=3)
    plt.setp(bp['whiskers'], color='black',linewidth=3)
    plt.setp(bp['medians'], color='black',linewidth=3)
    plt.setp(bp['caps'], color='black', linewidth=3)
    plt.setp(bp['fliers'], color='black', linewidth=3)
    plt.xticks([1, 2], ['Micropatterned cells', 'Standard cells'])
    plt.yticks(fontsize=50)
    plt.xticks(fontsize=26)
    plt.savefig(path.analysis_dir + "/analysis_spots_density/figures/nucleus_size_micropatterned_vs_standard.png",format='png')

def compare_cell_area(file_handler, sec_file_handler, acquisition1, acquisition2):
    arhgdia = build_image_list(file_handler, 'mrna', acquisition1)
    arhgdia_cultured = build_image_list(file_handler, 'mrna', acquisition2)
    cell_area_arhgdia = [idsc.get_cell_area(sec_file_handler, x) for x in arhgdia]
    cell_area_arhgdia_cultured = [idsc.get_cell_area(sec_file_handler, x) for x in arhgdia_cultured]
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    ax.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    bp = ax.boxplot([cell_area_arhgdia, cell_area_arhgdia_cultured])
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3)
    plt.setp(bp['boxes'], color='black',linewidth=3)
    plt.setp(bp['whiskers'], color='black',linewidth=3)
    plt.setp(bp['medians'], color='black',linewidth=3)
    plt.setp(bp['caps'], color='black', linewidth=3)
    plt.setp(bp['fliers'], color='black', linewidth=3)
    axes = plt.gca()
    axes.set_ylim([0, 2500])
    plt.xticks([1, 2], ['Micropatterned cells', 'Standard cells'])
    plt.yticks(fontsize=50)
    plt.xticks(fontsize=26)
    plt.savefig(path.analysis_dir + "/analysis_spots_density/figures/cell_size_micropatterned_vs_standard.png", format="png")
    plt.close()

def compare_cell_volume(file_handler, sec_file_handler, acquisition1, acquisition2):
    arhgdia = build_image_list(file_handler, 'mrna', acquisition1)
    arhgdia_cultured = build_image_list(file_handler, 'mrna', acquisition2)
    cell_volume_arhgdia = [idsc.compute_cell_volume(file_handler, image) for image in arhgdia]
    cell_volume_arhgdia_cultured = [idsc.compute_cell_volume(file_handler, image) for image in arhgdia_cultured]
    plt.boxplot([cell_volume_arhgdia, cell_volume_arhgdia_cultured])
    axes = plt.gca()
    axes.set_ylim([0, 2500])
    plt.xticks([1, 2], ['micropatterned', 'cultured'])
    plt.ylabel('Cell volume (um^3)')
    plt.title("Cell volume for micropatterned versus cultured")
    plt.savefig(path.analysis_dir + "/analysis_spots_density/figures/cell_volume_mic_vs_nuc.png",format='png')
    plt.close()



def compare_spots_density(file_handler, sec_file_handler, acquisition1, acquisition2):
    # Build all image for an acquisition
    arhgdia = build_image_list(file_handler, 'mrna', acquisition1)
    arhgdia_cultured = build_image_list(file_handler, 'mrna', acquisition2)
    num_spots_arhgdia = [len(idsc.get_spots(file_handler, x)) for x in arhgdia]
    num_spots_arhgdia_cultured = [len(idsc.get_spots(file_handler, x)) for x in arhgdia_cultured]
    cell_area_arhgdia = [idsc.get_cell_area(sec_file_handler, x) for x in arhgdia]
    cell_area_arhgdia_cultured=[]
    for x in arhgdia_cultured:
        nucleus_n=count_nucleus(file_handler, x)
        cell_area_arhgdia_cultured.append(idsc.get_cell_area(sec_file_handler, x)/nucleus_n)

    # display figures
    fig,(ax2,ax1)=plt.subplots(1,2,sharey=True,figsize=(15, 5))
    xs1 = cell_area_arhgdia
    ys1 = num_spots_arhgdia
    xs2 = cell_area_arhgdia_cultured
    ys2 = num_spots_arhgdia_cultured

    # Create linear regression object
    ax1.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.1)
    ax1.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    fit1 = np.polyfit(xs1, ys1, 1)
    fit_fn1 = np.poly1d(fit1)
    ax1.plot(xs1, ys1, 'yo', xs1, fit_fn1(xs1), '--k',c='#0080ff')
    ax2.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.1)
    ax2.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    fit2 = np.polyfit(xs2, ys2, 1)
    fit_fn2 = np.poly1d(fit2)
    ax2.plot(xs2, ys2, 'yo', xs2, fit_fn2(xs2), '--k',c='#646464')
    ax2.axis([200, 1500, -100, 600])
    ax1.axis([200, 800, -100, 600])

    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for axis in ['bottom', 'left']:
        ax1.spines[axis].set_linewidth(3)
        ax2.spines[axis].set_linewidth(3)

    plt.savefig(check_dir(path.analysis_dir + "/analysis_spots_density/figures/")+"mic_vs_cult.png", format="png")
    plt.close()

def compare_spots_volume_density(file_handler, sec_file_handler, acquisition1, acquisition2):
    # Build all image for an acquisition
    arhgdia = build_image_list(file_handler, 'mrna', acquisition1)
    arhgdia_cultured = build_image_list(file_handler, 'mrna', acquisition2)
    num_spots_arhgdia = [len(idsc.get_spots(file_handler, x)) for x in arhgdia]
    num_spots_arhgdia_cultured = [len(idsc.get_spots(file_handler, x)) for x in arhgdia_cultured]
    cell_volume_arhgdia = [idsc.compute_cell_volume(file_handler, image) for image in arhgdia]

    # cell_volume_arhgdia_cultured = [idsc.get_cell_area(file_handler, x) * const for x in arhgdia_cultured]
    # Use this cause we have sometime more than one nucleus, and consequently more than one cell
    # This will divide the area by the number of nucleus centroid
    cell_volume_arhgdia_cultured=[]
    for x in arhgdia_cultured:
        nucleus_n=count_nucleus(file_handler, x)
        cell_volume_arhgdia_cultured.append(idsc.compute_cell_volume(file_handler, x)/nucleus_n)

    # display figures
    fig,(ax2,ax1)=plt.subplots(1,2,sharey=True,figsize=(15, 5))
    xs1 = cell_volume_arhgdia
    ys1 = num_spots_arhgdia
    xs2 = cell_volume_arhgdia_cultured
    ys2 = num_spots_arhgdia_cultured

    # Create linear regression object
    ax1.scatter(xs1, ys1, color='black')
    fit1 = np.polyfit(xs1, ys1, 1)
    fit_fn1 = np.poly1d(fit1)
    ax1.plot(xs1, ys1, 'yo', xs1, fit_fn1(xs1), '--k')
    ax1.set_title('Micropatterned cells')
    ax2.scatter(xs2, ys2, color='black')
    fit2 = np.polyfit(xs2, ys2, 1)
    fit_fn2 = np.poly1d(fit2)
    ax2.set_title('Cultured cells')
    ax2.plot(xs2, ys2, 'yo', xs2, fit_fn2(xs2), '--k')
    ax1.set_xlabel('Cell volume (um^3)')
    ax2.axis([200, 1500, -100, 600])
    ax1.axis([200, 800, -100, 600])
    ax2.set_xlabel('Cell volume (um^3)')
    ax2.set_ylabel('Transcript number')
    plt.savefig(path.analysis_dir + "/analysis_spots_density/figures/mic_vs_cult.png", format="png")
    plt.show()

def compare_spots_density_by_gene_and_timepoint(file_handler, acquisition1, acquisition2, timepoints1, timepoints2):
    # build all image for an acquisition and a list of timepoint
    arhgdia_cultured = build_image_list_2(file_handler, 'mrna', acquisition2,timepoints2)
    arhgdia = build_image_list_2(file_handler, 'mrna', acquisition1,timepoints1)
    num_spots_arhgdia = [len(idsc.get_spots(file_handler, x)) for x in arhgdia]
    num_spots_arhgdia_cultured = [len(idsc.get_spots(file_handler, x)) for x in arhgdia_cultured]
    const = math.pow((1 / constants.SIZE_COEFFICIENT), 2) * 100
    cell_area_arhgdia = [idsc.get_cell_area(file_handler, x) * const for x in arhgdia]

    # cell_area_arhgdia_cultured = [idsc.get_cell_area(file_handler, x) * const for x in arhgdia_cultured]
    # Use this cause we have sometime more than one nucleus, and consequently more than one cell
    # This will divide the area by the number of nucleus centroid
    cell_area_arhgdia_cultured=[]
    for x in arhgdia_cultured:
        nucleus_n=count_nucleus(file_handler, x)
        cell_area_arhgdia_cultured.append(idsc.get_cell_area(file_handler, x)/nucleus_n)

    # display figures
    fig,(ax2,ax1)=plt.subplots(1,2,sharex=True,sharey=True,figsize=(15, 5))
    xs1 = cell_area_arhgdia
    ys1 = num_spots_arhgdia
    xs2 = cell_area_arhgdia_cultured
    ys2 = num_spots_arhgdia_cultured

    # Create linear regression object
    ax1.scatter(xs1, ys1, color='black')
    fit1 = np.polyfit(xs1, ys1, 1)
    fit_fn1 = np.poly1d(fit1)
    ax1.plot(xs1, ys1, 'yo', xs1, fit_fn1(xs1), '--k')
    ax1.set_title(acquisition1+ ' cells')
    ax2.scatter(xs2, ys2, color='black')
    fit2 = np.polyfit(xs2, ys2, 1)
    fit_fn2 = np.poly1d(fit2)
    ax2.set_title(acquisition2+' cells')
    ax2.plot(xs2, ys2, 'yo', xs2, fit_fn2(xs2), '--k')
    ax1.axis([200, 800, -100, 600])
    ax1.set_xlabel('Cell size (um^2)')
    ax2.axis([200, 800, -100, 600])
    ax2.set_xlabel('Cell size (um^2)')
    ax2.set_ylabel('Transcript number')
    plt.show()

def box_plot_fraction_profile(fractions1,fractions2,fraction1,fraction2,genes):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 10))
    bplot1 = axes[0].boxplot(fractions1,
                             vert=True,
                             patch_artist=True)
    bplot2 = axes[1].boxplot(fractions2,
                             vert=True,
                             patch_artist=True)
    colors =['#0A3950','#1E95BB','#A1BA6D','#F16C1B','#C02A18','#E9CB45']
    for bplot in (bplot1, bplot2):
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
    # adding horizontal grid lines
    axes[0].set_title(str(fraction1)+'%')
    axes[1].set_title(str(fraction2)+'%')
    for ax in axes:
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        ax.set_xticks([y + 1 for y in range(len(fractions1))], )
        ax.set_ylabel('Peripheral fraction')
    plt.setp(axes, xticks=[y + 1 for y in range(len(fractions1))],
             xticklabels=genes)
    plt.show()

def calculate_random_mpi(MTOC, nMTOC, bootstrap_mpi):
    scores=[]
    for value in MTOC:
        tmp_list=np.random.choice(np.array(nMTOC), bootstrap_mpi)
        med=np.median(tmp_list)
        scores.append(value-med)
    assert(len(MTOC)==len(scores))
    nneg = 0
    npos = 0
    for score in scores:
        if float(score) < 0:
            nneg += 1
        else:
            npos += +1
    mpi = (float(npos) / len(scores)) * 2 - 1
    p = binom.cdf(nneg, len(scores), 0.5)
    return (mpi, p)

def calculate_mpi(MTOC,nMTOC):
    scores = []
    tmp = MTOC - np.nanmedian(nMTOC)
    for elem in tmp.flatten():
        scores.append(elem)
    scores = np.array(scores)
    nneg = 0
    npos = 0
    for score in scores:
        if float(score) < 0:
            nneg += 1
        else:
            npos += 1
    mpi = (float(npos) / len(scores)) * 2 - 1
    p = binom.cdf(nneg, len(scores), 0.5)
    return (mpi, p)








