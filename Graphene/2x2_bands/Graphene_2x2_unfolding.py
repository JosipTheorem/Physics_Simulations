import os
import banduppy
import pickle

# Define the simulation folder and file prefix
sim_folder = r'G:\Praksa\bands_2x2'
prefix = os.path.join(sim_folder, 'graphene')


results_dir = sim_folder

# Define Fermi energy and energy range (choose these based on your SCF output)
Efermi = -2.6975    # Fermi energy from the 2x2 (unfolded) calculation
Emin = -20           # Minimum energy to plot (relative to Efermi)
Emax = 10            # Maximum energy to plot (relative to Efermi)
save_file_name = 'unfolded_bandstructure.png'

super_cell_size = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]

# Initialize the Unfolding object
band_unfold = banduppy.Unfolding(supercell=super_cell_size, print_log='high')

# Assume your saved file is located in results_dir
with open(f'{results_dir}/KPOINTS_SpecialKpoints.pkl', 'rb') as handle:
    special_kpoints_pos_labels = pickle.load(handle)

PC_BZ_path = [
    [0.0, 0.0, 0.0],  # Γ
    [0.33333, 0.33333, 0.0],  # K
    [0, 0.5, 0.0],  # M
    [0.0, 0.0, 0.0]   # Γ
]
# Provide a label for each node in the detailed k-path.
# (You might use more labels if your path is divided into more segments.)
special_k_points = "GKMG"  # Adjust labels as needed

# Specify a detailed number of k-points for each segment to get a continuous path.
npoints_per_path_seg = (84, 65, 51)


kpointsPBZ_full, kpointsPBZ_unique, kpointsSBZ, SBZ_PBZ_kpts_mapping, special_kpoints_pos_labels = \
    band_unfold.generate_SC_Kpts_from_pc_k_path(
        pathPBZ=PC_BZ_path,
        nk=npoints_per_path_seg,
        labels=special_k_points,
        kpts_weights=1,
        save_all_kpts=True,
        save_sc_kpts=True,
        save_dir=sim_folder,
        file_name_suffix='',
        file_format='qe'
    )

bands = banduppy.BandStructure(code="espresso", spinor=False, prefix=prefix)


# Perform the unfolding
results_dir = sim_folder
unfolded_bandstructure, kpline = band_unfold.Unfold(
    bands,
    kline_discontinuity_threshold=0.1,
    save_unfolded_kpts={'save2file': True, 'fdir': results_dir, 'fname': 'kpoints_unfolded'},
    save_unfolded_bandstr={'save2file': True, 'fdir': results_dir, 'fname': 'bandstructure_unfolded'}
)

print("Unfolding completed.")



# You can call the built-in plot function from your Unfolding object
fig, ax, CountFig = band_unfold.plot_ebs(
    save_figure_dir=results_dir,
    save_file_name=save_file_name,
    CountFig=None, 
    Ef=Efermi,         # Here you pass the Fermi level of your 2x2 calculation
    Emin=Emin,
    Emax=Emax,
    pad_energy_scale=0.5, 
    mode="density",    # Use density mode
    special_kpoints=special_kpoints_pos_labels, 
    plotSC=True, 
    fatfactor=20, 
    nE=100,
    smear=0.2, 
    marker='o',
    threshold_weight=0.01, 
    show_legend=True,
    color='gray', 
    color_map='viridis'
)




Efermi = -2.6975    # Fermi energy from the 2x2 (unfolded) calculation
Emin = -20           # Minimum energy to plot (relative to Efermi)
Emax = 10          # Maximum energy to plot (relative to Efermi)
save_file_name = 'unfolded_bandstructure_spectral.png'

fig, ax, CountFig \
= band_unfold.plot_ebs(save_figure_dir=results_dir, save_file_name=save_file_name, CountFig=None, 
                      Ef=Efermi, Emin=Emin, Emax=Emax, pad_energy_scale=0.5, 
                      mode="density", special_kpoints=special_kpoints_pos_labels, 
                      plotSC=True, fatfactor=20, nE=100,smear=0.2, marker='o',
                      threshold_weight=0.01, show_legend=True,
                      color='gray', color_map='viridis')