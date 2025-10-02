# Import Standard Modules
import numpy as np
import cupy as cp
import pandas as pd
import os

# Import self module relesed to PyPI(pip install massfunc bubblebarrier xxiop)
from xxiop.op import OpticalDepth
from .special import load_binary_data,TopHat_filter,xHII_field_update
from .ioninti_gpu import Ion

def reionization_calculator(fesc=0.2,A2byA1=1.0,kMpc_trans=1e6,alpha=0.0,beta=0.0,label = 'MH',DIM=256,box_length=800,save_on=True):

    ### Initialize the variables

    z_initial = np.arange(4.08, 30., 0.12)
    z_filtered = z_initial[(z_initial >= 5) & (z_initial <= 30)]
    z_value = np.concatenate((z_filtered[z_filtered <= 10], z_filtered[z_filtered > 10][::3]))
    z_value = z_value[::-1]  # reverse the order to have higher z first
    xHII_field = cp.zeros((DIM,DIM,DIM))
    nrec_field = cp.zeros((DIM,DIM,DIM))
    ionf = []
    
    ### main loop
    for i,z in enumerate(z_value):
        dz = z_value[i+1] - z_value[i] if i < len(z_value) - 1 else z_value[i] - z_value[i-1]
        ion = Ion(fesc=fesc,z=z,A2byA1=A2byA1,ktrans=kMpc_trans,alpha=alpha,beta=beta)
        read_path = f'df/updated_smoothed_deltax_z{z:06.2f}_{int(DIM)}_{int(box_length)}Mpc'
        ### load the density field
        delta_field_cpu = load_binary_data(read_path,DIM=DIM)
        delta_field = cp.asarray(delta_field_cpu)
        nrec_field += ion.dnrec_dz_path(delta_field,xHII_field) * dz
        del delta_field_cpu

        ### Partial ionization
        m_grid=ion.cosmo.rhom*(box_length/DIM)**3
        source = ion.nion_interp(m_grid,delta_field)
        igm = (1+nrec_field)*ion.n_HI(delta_field)
        minihalo = ion.nxi_interp(m_grid,delta_field)

        source_ratio = ion.nion_st(z) / cp.mean(source)
        minihalo_ratio = ion.nxi_st(z) / cp.mean(minihalo)

        partial_eff = source*source_ratio /(igm + minihalo*minihalo_ratio)
        del source,igm,minihalo

        ### main smoothing loop
        delta_field_ffted = cp.fft.rfftn(delta_field,norm="forward")
        nrec_field_ffted = cp.fft.rfftn(nrec_field,norm="forward")
        del delta_field

        rs = np.logspace(np.log10(box_length/DIM), np.log10(50), 50)
        rs=rs[::-1]
        ### Smooth the density field and recombination field at different scales
        for j,r in enumerate(rs):
            deltav_smoothed = TopHat_filter(delta_field_ffted,R=r,DIM=DIM,box_length=box_length)
            nrec_smoothed = TopHat_filter(nrec_field_ffted,R=r,DIM=DIM,box_length=box_length)
            mv=ion.cosmo.rhom*4*np.pi/3*r**3

            source_smoothed = ion.nion_interp(mv,deltav_smoothed)
            igm_smoothed = (1+nrec_smoothed)*ion.n_HI(deltav_smoothed)
            minihalo_smoothed = ion.nxi_interp(mv,deltav_smoothed)
            del deltav_smoothed,nrec_smoothed
            xHII_field[source_ratio*source_smoothed > ( igm_smoothed + minihalo_ratio*minihalo_smoothed)] = 1.0
            del source_smoothed,igm_smoothed,minihalo_smoothed


        os.makedirs(f'reionf/{label}', exist_ok=True)
        file_save_path = f'reionf/{label}/rf_{z:.2f}.npy'
        xHII_field = xHII_field_update(xHII_field,partial_eff)
        del partial_eff
        ionization_fraction_gpu = cp.mean(xHII_field)
        ionization_fraction = float(ionization_fraction_gpu)
        print(f"z = : {z:.2f} ({ionization_fraction*100:.2f}%)")
        if save_on:
            cp.save(file_save_path, xHII_field)
        ionf.append((z, ionization_fraction))

    final_z_values = np.array([item[0] for item in ionf])
    final_fractions = np.array([item[1] for item in ionf])

    # Calculate optical depth
    optical_depth = OpticalDepth(final_z_values[::-1], final_fractions[::-1])
    optical_depth_values = optical_depth.OpticalDepth(20)
    print(f"Optical depth at z=20: {optical_depth_values:.4f}")

    # 创建 DataFrame 并保存
    output_df = pd.DataFrame({
        'z': final_z_values,
        'ionf': final_fractions,
    })
    os.makedirs('csvfile', exist_ok=True)
    csv_save_path = f'csvfile/{label}.csv'
    output_df.to_csv(csv_save_path, index=False)
    return optical_depth_values
   
