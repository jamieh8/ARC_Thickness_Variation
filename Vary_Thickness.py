import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binned_statistic

from solcore import material, si
from solcore.material_system import create_new_material
from solcore.structure import Layer
from solcore.absorption_calculator import calculate_rat, OptiStack, download_db, search_db, create_nk_txt
from solcore.absorption_calculator import create_nk_txt, download_db, search_db
from solcore.config_tools import add_source
import os

def make_ARC_stack(ARC_materials: list, ARC_thicknesses_nm: list, base, substrate_mat):
    '''
    :param ARC_materials: Array of Solcore materials for ARC layers. A material can be a list of arrays: [wavelengths[in nm], ns, ks]
    :param ARC_thicknesses_nm: Array of thicknesses, in nm, for each ARC layer. Index 0 is the top layer (from air). Last element is above the "base" layer(s).
    :param base: List of layer(s) on which to add the ARC layers
    :param substrate_mat: Solcore material to use as semi-inf substrate behind stack
    :return: OptiStack object with ARC layers + base, and substrate set as material provided.
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
    '''
    Calculates RAT of the stack provided, returns R.
    :param stack: Solcore Structure or OptiStack
    :param wls_nm: Array of wavelengths, in nm
    :return: R of stack at each wavelength, as an array
    '''
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



wls_nm = np.arange(350, 1610, 10)  # wavelengths in nm. [350-1600] used in Yan et al 2012. 10 nm steps for optimization.


# Set up materials, base
SiO2 = material('SiO2')()
TiO2 = material('TiO2')()
# print(f'SiO2 n at 500nm: {SiO2.n_interpolated(500*1e-9)}')
# print(f'TiO2 n at 500nm: {TiO2.n_interpolated(500*1e-9)}')

window_mat = material('AlInP')(Al=0.5)
InGaP = material('GaInP')(In=0.5)

window_layer = [Layer(si('30nm'), material=window_mat)]
base_layers_reg = window_layer

# custom TiO2
# download_db()
# results = search_db('TiO2')
# create_nk_txt(pageid=479, file = 'TiO2_Devore-o')
# create_nk_txt(pageid=482, file = 'TiO2_Sakar')
# create_new_material(mat_name = 'TiO2_Devore-o', n_source='TiO2_Devore-o_n.txt', k_source='TiO2_Sakar_k.txt')
# create_new_material(mat_name = 'TiO2_Sakar', n_source='TiO2_Sakar_n.txt', k_source='TiO2_Sakar_k.txt')
# create_new_material(mat_name = 'TiO2_Zhukovsky', n_source='TiO2_Zhukovsky_n.txt', k_source='TiO2_Sakar_k.txt')
TiO2_Sakar = material('TiO2_Sakar')()
TiO2_Devore = material('TiO2_Devore-o')()  # doesn't extend to low enough wavelengths
TiO2_Zhukovsky = material('TiO2_Zhukovsky')()


# check refractive indices
plt.figure(3)
wls_nm_forn = np.arange(300,1000,10) * 1e-9
# plt.plot(wls_nm_forn*1e9, window_mat.n(wls_nm_forn), label='window, AlInP')
# plt.plot(wls_nm_forn*1e9, InGaP.n(wls_nm_forn), label='InGaP')

plt.plot(wls_nm_forn*1e9, TiO2.n(wls_nm_forn), label='TiO2 default', c='crimson')
plt.plot(wls_nm_forn*1e9, TiO2_Sakar.n(wls_nm_forn), label='TiO2 Sakar', c='purple')
# plt.plot(wls_nm_forn*1e9, TiO2_Zhukovsky.n(wls_nm_forn), label='TiO2 Zhukovsky', c='orangered')
plt.plot([500],[2.41], 'xk', label='Yan et al 2012')

plt.xlim([300,1000])
# plt.ylim([0,5])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Refractive index, n')
plt.legend()



# Set up designs (ARC materials & nominal thicknesses) to test

designs = []
# list of dictionaries. each dict is 1 "design".
# to plot only existing data, set number of samples to 0

# 2 layer ARC
ARC_mats_double = [SiO2,TiO2]
ARC_ths_ref_double = [102,45]  # [nm], reference/nominal thicknesses
designs += [{'ARC mats':ARC_mats_double, 'ARC nom ths':ARC_ths_ref_double,
             'base layers':base_layers_reg,
             'colour':'crimson', 'label':'2 layer',
             'figure of merit':'averaged R', 'discrete stddevs':False, 'number of samples':0}]  # Double layer

# 2 layer ARC - TiO2 sakar
ARC_mats_2lay_TiO2sakar = [SiO2,TiO2_Sakar]
ARC_ths_ref_double = [102,45]  # [nm], reference/nominal thicknesses
designs += [{'ARC mats':ARC_mats_2lay_TiO2sakar, 'ARC nom ths':ARC_ths_ref_double,
             'base layers':base_layers_reg,
             'colour':'purple', 'label':'2 layer - TiO2 Sakar',
             'figure of merit':'averaged R', 'discrete stddevs':False, 'number of samples':0}]  # Double layer



# 4 layer ARC - not varying n by wavelength
# l2_mat = [wls_nm, len(wls_nm)*[1.77], SiO2.k(wls_nm*1e-9)]
# l3_mat = [wls_nm, len(wls_nm)*[1.34], SiO2.k(wls_nm*1e-9)]
# l4_mat = [wls_nm, len(wls_nm)*[1.09], SiO2.k(wls_nm*1e-9)]
# ARC_mats_4l = [l4_mat, l3_mat, l2_mat, TiO2]
# ARC_ths_4l = [225,119,80,46]
# designs += [{'ARC mats':ARC_mats_4l, 'ARC nom ths':ARC_ths_4l,
#              'colour':'teal', 'label':'4 layer',
#             'figure of merit':'averaged R', 'discrete stddevs':False, 'number of samples':0}]  # 4 layer ARC


# 4 layer ARC, with variable n and k from composition
# wls_m = wls_nm * 1e-9
# # x calculation based on n at 500 nm wavelength:
# # if n > n_SiO2 : n_new = x n_SiO2 + (1-x) n_TiO2
# # if n < n_SiO2 : n_new = x n_SiO2 + (1-x) n_air --- n_air = 1, k_air = 0
#
# x2 = (1.77-2.41)/(1.46-2.41)  # n TiO2 at 500 nm was 2.41 in Yan et al, 3.03 from Solcore
# l2_mat = [wls_nm, x2*SiO2.n(wls_m) + (1-x2)*TiO2.n(wls_m), x2*SiO2.k(wls_m) + (1-x2)*TiO2.k(wls_m)]
#
# x3 = (1.34-1)/(1.46-1)
# l3_mat = [wls_nm, x3*SiO2.n(wls_m) + (1-x3)*np.ones(len(wls_nm)), x3*SiO2.k(wls_m) + (1-x3)*np.zeros(len(wls_m))]
#
# x4 = (1.09-1)/(1.46-1)
# l4_mat = [wls_nm, x4*SiO2.n(wls_m) + (1-x4)*np.ones(len(wls_nm)), x4*SiO2.k(wls_m) + (1-x4)*np.zeros(len(wls_m))]
#
# ARC_mats_4l = [l4_mat, l3_mat, l2_mat, TiO2]
# ARC_ths_4l = [225,119,80,46]
# designs += [{'ARC mats':ARC_mats_4l, 'ARC nom ths':ARC_ths_4l,
#              'base layers':base_layers_reg,
#              'colour':'green', 'label':'4 layer - variable comp',
#             'figure of merit':'averaged R', 'discrete stddevs':False, 'number of samples':5}]  # 4 layer ARC



# Set up plotting options
plot_ref_R = True  # creates a separate plot containing R vs wavelength for the reference/unvaried structures
add_comparison_plot = True  # adds a subplot to Fig 1, plotting mean and stddev of each design on the same y axis
add_comparison_delta_plot = True  # adds a subplot to Fig 1, plotting ^ but with Delta (mean-ref)

n_subplots = len(designs)
if add_comparison_plot:
    compare_plt_index = n_subplots
    n_subplots += 1
if add_comparison_delta_plot:
    compare_delta_plt_index = n_subplots
    n_subplots += 1


fig, axs = plt.subplots(1,n_subplots, sharey=False, layout='tight', num=1)


max_stddev = 15
discrete_stddevs = np.arange(1, max_stddev+1, 1)  # discrete standard deviations (1-15) to sample (if option is selected)


for i, design in enumerate(designs):
    # UNVARIED REFERENCE: Get and plot reference R, and reference figure of merit
    ARC_mats = design['ARC mats']
    ARC_ths = design['ARC nom ths']
    base_layers = design['base layers']
    ref_stack = make_ARC_stack(ARC_materials=ARC_mats, ARC_thicknesses_nm=ARC_ths, base=base_layers, substrate_mat=InGaP)
    ref_R = get_R_of_stack(ref_stack, wls_nm)
    avgR_ref = fig_of_merit(ref_R, wls_nm, design['figure of merit'])
    print(design['label'] + f' fig of merit ref: {avgR_ref}')
    if plot_ref_R:
        plt.figure(2)
        plt.plot(wls_nm, ref_R*100, '--', color=design['colour'], label=design['label'])


    # SAMPLE VARIATIONS: Vary layer thicknesses of reference with diff standard deviations
    number_of_samples = design['number of samples']  # per stddev value, if discrete_stddevs is true. otherwise, total number of samples
    use_discrete_stddevs = design['discrete stddevs']
    new_rows = []

    for si in range(number_of_samples):
        if use_discrete_stddevs:
            d_stdvs = discrete_stddevs
        else:
             d_stdvs = [np.random.rand()*max_stddev]  # random stand dev between 0 and 15

        # 1 loop for random sampling. 1 loop PER stddev for discrete sampling.
        for dist_stddev in d_stdvs:
            new_thicknesses_nm = add_variation_to_layers(ARC_ths, stand_dev=dist_stddev)
            mod_stack = make_ARC_stack(ARC_materials=ARC_mats, ARC_thicknesses_nm=list(new_thicknesses_nm), base=base_layers, substrate_mat=InGaP)
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


    # retrieve axis for scatter plot
    if n_subplots == 1:
        # only 1 design, no comparison plot.
        ax_scat = axs
    else:
        ax_scat = axs[i]

    # scatter plot
    ax_scat.plot(stddevs, fig_of_merits, '.', c=design['colour'], label=design['label'])

    # reference line for unvaried thickness fig of merit
    ax_scat.plot([-1, max_stddev+1], 2 * [avgR_ref], '--', c='black', lw=2)  # design['colour'])


    # BINNED STATISTICS: Calculate binned mean, if random sampling was used
    if use_discrete_stddevs == False:
        bin_size = 1  # [%]
        bins = np.arange(0, max_stddev + bin_size, bin_size)  # bin edges
        dist_stddevs = np.arange(bin_size/2,max_stddev,bin_size)  # centers of bins
    else:
        # if discrete sampling, calculate mean for each discrete stddev
        bins = np.arange(0.5, max_stddev+1, 1)  # bins edges, so center is integers
        dist_stddevs = discrete_stddevs  # centers of bins = integer stddevs sampled

    binned_mean = binned_statistic(stddevs, fig_of_merits, statistic='mean', bins=bins)
    binned_stddev = binned_statistic(stddevs, fig_of_merits, statistic='std', bins=bins)
    fom_means, fom_stddevs = binned_mean.statistic, binned_stddev.statistic

    # mean and error bars in scatter plot
    ax_scat.errorbar(dist_stddevs, fom_means, yerr=fom_stddevs, fmt = '-sk')
    ax_scat.fill_between(dist_stddevs,
                        y1 = fom_means + fom_stddevs,
                        y2 = fom_means - fom_stddevs,
                        color = 'black', alpha=0.1)
    ax_scat.set_xlabel('Standard Dev [%]')
    ax_scat.set_xlim([0, max_stddev])
    yl_str = design['figure of merit']
    ax_scat.set_ylabel(f'{yl_str} [%]')
    ax_scat.legend()


    # mean and error "shadow" in comparison plot
    if add_comparison_plot:
        axs[compare_plt_index].plot(dist_stddevs, fom_means, color=design['colour'])
        axs[compare_plt_index].fill_between(dist_stddevs,
                        y1 = fom_means + fom_stddevs,
                        y2 = fom_means - fom_stddevs,
                        color = design['colour'], alpha=0.2)
        axs[compare_plt_index].plot([-1, max_stddev+1], 2 * [avgR_ref], '--', c=design['colour'])
        axs[compare_plt_index].set_xlim([0,max_stddev])

    if add_comparison_delta_plot:
        delta_from_ref = fom_means - avgR_ref
        axs[compare_delta_plt_index].plot(dist_stddevs, delta_from_ref, c=design['colour'])
        axs[compare_delta_plt_index].fill_between(dist_stddevs,
                                                  y1 = delta_from_ref + fom_stddevs,
                                                  y2 = delta_from_ref- fom_stddevs,
                                                  color = design['colour'], alpha=0.2)


if add_comparison_plot:
    axs[compare_plt_index].set_xlim([0,max_stddev])
    axs[compare_plt_index].set_xlabel('Standard Dev [%]')
    axs[compare_plt_index].set_ylabel('averaged R [%]')  # label should depend on figure of merit used

if add_comparison_delta_plot:
    axs[compare_delta_plt_index].plot([-1, max_stddev + 1], [0, 0], '-k', lw=0.8)
    axs[compare_delta_plt_index].set_xlim([0,max_stddev])
    axs[compare_delta_plt_index].set_xlabel('Standard Dev [%]')
    axs[compare_delta_plt_index].set_ylabel('$\Delta$ averaged R [%]')


if plot_ref_R:
    plt.figure(2)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('R [%]')
    plt.ylim([0,45])
    plt.title('Reflectance with no variation')
    plt.legend()

fig.suptitle('Effect of Normally Distributed Layer Variation')


plt.show()

