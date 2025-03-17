import arc_2 as anp
import openmc
import numpy as np
import os
import time
import sys
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
from openmc_plasma_source import TokamakSource
from openmc_source_plotter import plot_source_direction, plot_source_position


openmc.config["cross_sections"] = '/home/hice1/kthompson309/endfb-viii.0-hdf5/cross_sections.xml'


# ==============================================================================
# Geometry
# ==============================================================================
def create_arc(Li6_enrichment):
    device = anp.generate_device("U", 0, Li6_enrichment = Li6_enrichment)
    
    # Plotting
    plot = openmc.Plot()
    plot.filename = 'geometry_plot'
    plot.basis = 'xz'
    plot.origin = (350, 0, 0)
    plot.width = (700, 800)
    plot.pixels = (plot.width[0]*10, plot.width[1]*10)
    
    color = ['beige', 'lightcoral', 'yellow', 'orange', 'lime', 'navy', 'lightcyan', 'black', 'salmon', 'magenta', 'green']
    color_dict = {cell.id: color[i] for i, cell in enumerate(device._cells)}
    plot.colors = color_dict
    for count, cell in enumerate(device._cells):
        print(f'Device name: {cell.name} with color: {color[count]}')
    
    #plot.highlight_domains(geometry=device.geometry, domains=device._cells)
    
    plots = openmc.Plots([plot])
    plots.export_to_xml()
    
    # ==============================================================================
    # Settings
    # ==============================================================================
    
    """ Source Definition """
    my_sources = TokamakSource(
    elongation=1.84, # jball code says 1.5
    ion_density_centre=1.8e20,# paper
    ion_density_peaking_factor=1,
    ion_density_pedestal=1.09e20,
    ion_density_separatrix=3e19,
    ion_temperature_centre=27,# paper
    ion_temperature_peaking_factor=8.06,
    ion_temperature_pedestal=6.09,
    ion_temperature_separatrix=0.1,
    major_radius=3.3,#4.0 jball in code
    minor_radius=1.15, #1.2 jball in code
    pedestal_radius=0.8 * 1.2,# jball
    mode="H",
    shafranov_factor=0.44789,# good enough
    triangularity=0.5,# jball
    ion_temperature_beta=6,
    sample_size=50,  # the number of individual sources to use to make a combined source
    # angles=( -3.141592 / 18 , 3.141592 / 18)  # angle in radians
    angles=(0.001, 3.14195 / 18)
    ).make_openmc_sources()  # returns a list of openmc sources
    # source = openmc.Source()
    # source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(400, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18), openmc.stats.Discrete(0, 1)) # original openmc.stats.Discrete(450, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18)
    # source.angles = openmc.stats.Isotropic()
    # source.energy = openmc.stats.Discrete([14.1E6], [1.0])


    # this is example code to plot the source
    # settings = openmc.Settings()
    # settings.particles = 1
    # settings.batches = 1
    # settings.source = my_sources
    # materials = openmc.Materials()
    # sph = openmc.Sphere(r=1000000, boundary_type="vacuum")
    # cell = openmc.Cell(region=-sph)
    # geometry = openmc.Geometry([cell])
    # model = openmc.Model(geometry, materials, settings)

    # plot = plot_source_position(this=model, n_samples=2000)
    # plot.show()
    # # plot.savefig('source_position.png')

    # plot = plot_source_direction(this=model, n_samples=500)
    # plot.show()
    # # plot.savefig('source_direction.png')`
    
    device.settings.source = my_sources
    # energy filter
    # energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")
    # ==============================================================================
    # Blanket Material
    # ==============================================================================
    # breeding_material = openmc.Material(material_id = 56)  # Pb84.2Li15.8
    # breeding_material.add_element('Pb', 84.2)
    # breeding_material.add_element('Li', 15.8)
    # breeding_material.set_density('g/cm3', 11.)
    # device.blanket.fill = breeding_material
    # ==============================================================================
    # Tallies
    # ==============================================================================
    # """ Cylindrical Mesh Tally """
    r_grid = np.linspace(0, 600, num=25) #[NEW] original: (25, 200, num=25), 0, 600
    z_grid = np.linspace(-700, 700, num=50) #[NEW] original: (-200, 200, num=50), -700, 700
    mesh = openmc.CylindricalMesh(r_grid=r_grid, z_grid=z_grid) #[NEW]
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)
    
    device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])
    
    # """ TBR Tally """
    # tbr_filter1 = openmc.MaterialFilter(anp.tungsten)
    # device.add_tally('Tbr Plasma-facing Component Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter1])
    
    # tbr_filter2 = openmc.MaterialFilter(device.vcrti_VV)
    # device.add_tally('Tbr Vacuum Vessel Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter2])
    
    tbr_filter3 = openmc.MaterialFilter(device.doped_flibe_channels)
    device.add_tally('Tbr Channel Tally ', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter3])
    
    # tbr_filter4 = openmc.MaterialFilter(device.vcrti_BI)
    # device.add_tally('Tbr Tank Inner Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter4])

    tbr_filter5 = openmc.CellFilter(device.get_cell(name = 'blanket'))
    # tbr_filter5 = openmc.MaterialFilter(device.doped_flibe_blanket) #device.doped_flibe, initially wanted 'doped_mat' which doesn't exist
    device.add_tally('Tbr Blanket Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter5])
    
    # tbr_filter6 = openmc.MaterialFilter(device.vcrti_BO)
    # device.add_tally('Tbr Tank Outer Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter6])
    
    # ==============================================================================
    # Run
    # ==============================================================================
    
    device.settings.photon_transport = True
    device.survival_biasing = True
    device.build()
    device.export_to_xml(remove_surfs=True)
    
    #openmc.plot_geometry()
    
    # set run parameters
    #device.settings.threads = 1 
    device.settings.particles = int(1e6)
    device.settings.batches = 10  
    device.settings.inactive = 1  
    
    # device.settings.mpi_args = ['mpiexec', '-np', '4']
    #remove old output files
    # for file in os.listdir('.'):
    #     if file.endswith('.h5'):
    #         os.remove(file)
    # # Run the simulation & set output file to a variable
    # out_file = device.run()
    return device
# ================================================================================
# additional runs
# ================================================================================

def make_materials_geometry_tallies(Li6_enrichment):
    """Makes a neutronics model of a blanket and simulates the TBR value.

    Arguments:
        enrichment (float): the enrichment percentage of Li6 in the breeder material
    
    Returns:
        resutsl (dict): simulation tally results for TBR along with the standard deviation and enrichment
    """
    start_time = time.time()
    # RUN OPENMC
    device = create_arc(Li6_enrichment)
    print(device.Li6_enrichment)
    
    #remove old output files
    #for file in os.listdir('.'):
    #    if file.endswith('.h5'):
    #        os.remove(file)
    sp_filename = device.run(output = False, threads=24)  # runs with reduced amount of output printing

    # OPEN OUPUT FILE
    sp = openmc.StatePoint(sp_filename)

    tbr_tally = sp.get_tally(name='Tbr Blanket Tally')

    df = tbr_tally.get_pandas_dataframe()
    print(df)
    df.to_csv(f'dataframe{device.Li6_enrichment}.csv')
    tbr_tally_result = df['mean'].sum()
    tbr_tally_std_dev = df['std. dev.'].sum()

    # command = ["openmc-plot-mesh-tally", sp_filename]
    # # Run the command
    # subprocess.run(command)
    print(f'time: {start_time - time.time()}')
    return {'enrichment': device.Li6_enrichment,
            'tbr_tally_result': tbr_tally_result,
            'tbr_tally_std_dev': tbr_tally_std_dev}

#results = []
#for enrichment in [0.01, 7.5, 15, 25, 50, 75, 99.99]:  # percentage enrichment from 0% Li6 to 100% Li6
#    results.append(make_materials_geometry_tallies(enrichment))
#print(results)
#results.append(make_materials_geometry_tallies(7.5))
# PLOTS RESULTS
#x = [entry['enrichment'] for entry in results]
#y = [entry['tbr_tally_result'] for entry in results]
#error_y = {'array': [entry['tbr_tally_std_dev'] for entry in results]}



result = make_materials_geometry_tallies(float(sys.argv[1]))
np.save(f'enrichment_{float(sys.argv[1])}.npy', result)
# read_dictionary = np.load('enrichment_7.5.npy',allow_pickle='TRUE').item()
# df = pd.DataFrame(data=result)
#df.to_csv(f'/home/hice1/kthompson309/TBR/output_df/enrichment_{sys.argv[1]}')
#plt.plot(x, y)
#plt.title="TBR as a function of Li6 enrichment",
#plt.xtitle="Li6 enrichment (%)",
#plt.ytitle="TBR"
#plt.show()
# ================================================================================
try:
    if sys.argv[1] is not None:
        os.mkdir(str(sys.argv[1]))
        device.move_files(str(sys.argv[1]))
        print("OpenMC files moved to new directory:", str(sys.argv[1]))

except:
    print("No directory specified, using this one")

# =============================================
# tally plot
# =============================================
# out_file = ""
# for file in os.listdir('.'):
#     if file.endswith('.h5'):
#         if file != "summary.h5":
#             out_file = file
# command = ["openmc-plot-mesh-tally", out_file]
# # Run the command
# subprocess.run(command)

# # open the results file
# sp = openmc.StatePoint(out_file)
# # access the tally using pandas dataframes
# tbr_tally = sp.get_tally(name='Tbr Blanket Tally')
# df = tbr_tally.get_pandas_dataframe()
# # prints the contents of the dataframe
# print(df)
