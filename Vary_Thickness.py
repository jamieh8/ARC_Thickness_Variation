import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binned_statistic

from solcore import material, si
from solcore.structure import Layer
from solcore.absorption_calculator import calculate_rat, OptiStack


def make_ARC_stack(ARC_materials, ARC_thicknesses_nm, base, substrate_mat):
    '''
    :param ARC_materials: Array of Solcore materials for ARC layers
    :param ARC_thicknesses_nm: Array of thicknesses, in nm, for each ARC layer
    :param base: List of layer(s) on which to add the ARC layers
    :param substrate_mat: Solcore material to use as semi-inf substrate behind stack
    :return: OptiStack object with ARC layers + base
    '''
    # make ARC layers
    ARC_layers = []
    for mat, l_thick_nm in zip(ARC_materials, ARC_thicknesses_nm):
        if type(mat) == list:
            ARC_layers += [[l_thick_nm]+mat]  # if material is a list [wls[nm], ns, ks], add thickness as first element to create layer
        else:
            ARC_layers += [Layer(l_thick_nm*1e-9, material = mat)]  # if


    # make and return stack
    return OptiStack(ARC_layers + base, substrate=substrate_mat)


def get_R_of_stack(stack, wls_nm):
    RAT = calculate_rat(stack, wls_nm, no_back_reflection=False)
    return RAT['R']

def fig_of_merit(R, wls, type_str):
    if type_str == 'integrated R':
        return np.trapz(y=R, x=wls_nm)

    elif type_str == 'averaged R':
        return 100*np.average(R)  # as [%]

    elif type_str == 'weighted averaged R':
        R_weights = np.where(wls<700, 4, 1)  # Rs for wls < 700 nm
        return 100*np.average(R, weights=R_weights)  # as [%]


def add_variation_to_layers(layer_ths, stand_dev):
    '''
    Returns a modified list of layers thicknesses, which is varied by some % from the original thicknesses provided in layer_ths.
    Variation, in % of original layer thickness, is sampled from normal distribution.

    :param layer_ths: Array of layer thicknesses to vary
    :param stand_dev: Standard deviation [%] of normal distribution from which to sample variation
    :return: Array of thicknesses with some random variation (added or subtracted)
    '''
    frac_deviation = np.random.normal(loc=0, scale=stand_dev, size=len(layer_ths))  #  get array of % deviations from nominal
    new_thicknesses = (0.01*frac_deviation+1)*layer_ths
    return new_thicknesses


def filename_from_design_dict(design_dict):
    label = design_dict['label']
    disc_stddev = design_dict['discrete stddevs']
    fig_of_merit = design_dict['figure of merit']

    filename = f'ThicknessVariation_LABEL-{label}_FIGOFMERIT-{fig_of_merit}_DISCRETE-{disc_stddev}.csv'

    return filename


DBR_list = []


# Set up materials, base
SiO2 = material('SiO2')()
TiO2 = material('TiO2')()
# print(SiO2.n_interpolated(500*1e-9))

window_mat = material('AlInP')(Al=0.52)
InGaP = material('GaInP')(In=0.49)
window_layer = [Layer(si('30nm'), material=window_mat)]

wls_nm = np.arange(350, 1600, 10)  # wavelengths in nm

# Designs to test
designs = []
# Note: to plot only existing data, set number of samples to 0

# 2 layer ARC
ARC_mats_double = [SiO2,TiO2]
ARC_ths_ref_double = [102,45]  # [nm], reference/nominal thicknesses
designs += [{'ARC mats':ARC_mats_double, 'ARC nom ths':ARC_ths_ref_double,
             'colour':'crimson', 'label':'2 layer',
             'figure of merit':'averaged R', 'discrete stddevs':False, 'number of samples':0}]  # Double layer

# 4 layer ARC
# l2_mat = [wls_nm, len(wls_nm)*[1.77], SiO2.k(wls_nm*1e-9)]
# l3_mat = [wls_nm, len(wls_nm)*[1.34], SiO2.k(wls_nm*1e-9)]
# l4_mat = [wls_nm, len(wls_nm)*[1.09], SiO2.k(wls_nm*1e-9)]
# ARC_mats_4l = [l4_mat, l3_mat, l2_mat, TiO2]
# ARC_ths_4l = [225,119,80,46]
# designs += [{'ARC mats':ARC_mats_4l, 'ARC nom ths':ARC_ths_4l,
#              'colour':'teal', 'label':'4 layer',
#             'figure of merit':'averaged R', 'discrete stddevs':False, 'number of samples':0}]  # 4 layer ARC


# 4 layer ARC, with variable n and k from composition
wls_m = wls_nm * 1e-9
# x calculation based on n at 500 nm wavelength:
# if n > n_SiO2 : n_new = x n_SiO2 + (1-x) n_TiO2
# if n < n_SiO2 : n_new = x n_SiO2 + (1-x) n_air --- n_air = 1, k_air = 0

x2 = (1.77-2.41)/(1.46-2.41)
l2_mat = [wls_nm, x2*SiO2.n(wls_m) + (1-x2)*TiO2.n(wls_m), x2*SiO2.k(wls_m) + (1-x2)*TiO2.k(wls_m)]

x3 = (1.34-1)/(1.46-1)
l3_mat = [wls_nm, x3*SiO2.n(wls_m) + (1-x3)*np.ones(len(wls_nm)), x3*SiO2.k(wls_m) + (1-x3)*np.zeros(len(wls_m))]

x4 = (1.09-1)/(1.46-1)
l4_mat = [wls_nm, x4*SiO2.n(wls_m) + (1-x4)*np.ones(len(wls_nm)), x4*SiO2.k(wls_m) + (1-x4)*np.zeros(len(wls_m))]

ARC_mats_4l = [l4_mat, l3_mat, l2_mat, TiO2]
ARC_ths_4l = [225,119,80,46]
designs += [{'ARC mats':ARC_mats_4l, 'ARC nom ths':ARC_ths_4l,
             'colour':'green', 'label':'4 layer - variable comp',
            'figure of merit':'averaged R', 'discrete stddevs':False, 'number of samples':0}]  # 4 layer ARC



# Set up plotting options
n_subplots = len(designs)
add_comparison_plot = True
plot_ref_R = True
if add_comparison_plot:
    compare_plt_index = n_subplots
    n_subplots += 1
if plot_ref_R:
    refR_plt_index = n_subplots
    n_subplots += 1

fig, axs = plt.subplots(1,n_subplots, sharey=False, layout='tight')


discrete_stddevs = np.arange(1, 16, 1)  # discrete standard deviations (1-15) to sample (if option is selected)


for i, design in enumerate(designs):
    # get and plot reference R (not varied)
    ARC_mats = design['ARC mats']
    ARC_ths = design['ARC nom ths']
    ref_stack = make_ARC_stack(ARC_materials=ARC_mats, ARC_thicknesses_nm=ARC_ths, base=window_layer, substrate_mat=InGaP)
    ref_R = get_R_of_stack(ref_stack, wls_nm)
    avgR_ref = fig_of_merit(ref_R, wls_nm, design['figure of merit'])
    axs[i].plot([-1,16], 2*[avgR_ref], '--', c=design['colour'])

    # vary layer thicknesses with diff standard deviations
    number_of_samples = design['number of samples']  # per stddev value, if discrete_stddevs is true. otherwise, total number of samples
    use_discrete_stddevs = design['discrete stddevs']
    new_rows = []

    # for each sample,
    for si in range(number_of_samples):
        if use_discrete_stddevs:
            d_stdvs = discrete_stddevs
        else:
             d_stdvs = [np.random.rand()*15]  # random stand dev between 0 and 15

        # 1 loop for random sampling. 1 loop PER stddev for discrete sampling.
        for dist_stddev in d_stdvs:
            new_thicknesses_nm = add_variation_to_layers(ARC_ths, stand_dev=dist_stddev)
            mod_stack = make_ARC_stack(ARC_materials=ARC_mats, ARC_thicknesses_nm=list(new_thicknesses_nm), base=window_layer, substrate_mat=InGaP)
            mod_R = get_R_of_stack(mod_stack, wls_nm)
            avgR_mod = fig_of_merit(mod_R, wls_nm, design['figure of merit'])
            row = [dist_stddev, avgR_mod]
            new_rows += [row]

    # retrieve and save data
    fnm = filename_from_design_dict(design)
    try:
        # check if there is existing data, and append new results to it
        prev_rows = np.loadtxt(fname=fnm, dtype=float, delimiter=',')
        if len(new_rows) == 0:
            # no new data, use only prev
            all_rows = prev_rows
        else:
            # new and prev data: append new to prev
            all_rows = np.append(prev_rows, np.array(new_rows), axis=0)
    except:
        # no prev data, use only new
        all_rows = np.array(new_rows)


    if len(all_rows)!= 0:
        np.savetxt(fname=fnm, X=all_rows, delimiter=',')

    stddevs = all_rows[:,0]
    fig_of_merits = all_rows[:,1]


    # scatter plot
    axs[i].plot(stddevs, fig_of_merits, '.', c=design['colour'], label=design['label'])

    # calculate binned mean, if random sampling was used
    if use_discrete_stddevs == False:
        bin_size = 1  # [%]
        bins = np.arange(0, 15 + bin_size, bin_size)  # bin edges
        dist_stddevs = np.arange(bin_size/2,15,bin_size)  # centers of bins
    else:
        # if discrete sampling, calculate mean for each discrete stddev
        bins = np.arange(-0.5, 16, 1)  # bins edges, so center is integers
        dist_stddevs = discrete_stddevs  # centers of bins = integer stddevs sampled

    binned_mean = binned_statistic(stddevs, fig_of_merits, statistic='mean', bins=bins)
    binned_stddev = binned_statistic(stddevs, fig_of_merits, statistic='std', bins=bins)
    jsc_means, jsc_stddevs = binned_mean.statistic, binned_stddev.statistic

    # mean and error bars in scatter plot
    axs[i].errorbar(dist_stddevs, jsc_means, yerr=jsc_stddevs, fmt = '-sk')
    axs[i].fill_between(dist_stddevs,
                        y1 = jsc_means + jsc_stddevs,
                        y2 = jsc_means - jsc_stddevs,
                        color = 'black', alpha=0.2)

    # mean and error "shadow" in "aside" comparison plot
    if add_comparison_plot:
        axs[compare_plt_index].plot(dist_stddevs, jsc_means, color=design['colour'])
        axs[compare_plt_index].fill_between(dist_stddevs,
                        y1 = jsc_means + jsc_stddevs,
                        y2 = jsc_means - jsc_stddevs,
                        color = design['colour'], alpha=0.2)

    if plot_ref_R:
        axs[refR_plt_index].plot(wls_nm, ref_R*100, '--', color=design['colour'])

    axs[i].set_xlabel('Standard Dev [%]')
    axs[i].set_xlim([0,15])
    yl_str = design['figure of merit']
    axs[i].set_ylabel(f'{yl_str} [%]')
    axs[i].legend()

if add_comparison_plot:
    axs[compare_plt_index].set_xlim([0,15])
    axs[compare_plt_index].set_xlabel('Standard Dev [%]')

if plot_ref_R:
    axs[refR_plt_index].set_xlabel('Wavelength [nm]')
    axs[refR_plt_index].set_ylabel('R [%]')
    axs[refR_plt_index].set_ylim([0,45])

fig.suptitle('Effect of Normally Distributed Layer Variation')


plt.show()

