# Code to apply fog
# Based on Koschmieder's principle
# https://www.researchgate.net/publication/320668349_Towards_Simulating_Foggy_and_Hazy_Images_and_Evaluating_Their_Authenticity


import numpy as np
import cv2
import pathlib
import time
import scipy.special


from colour_demosaicing import (
    demosaicing_CFA_Bayer_Malvar2004,
)

##########################################################################################################
## Initialise Classes for camera and weather
class Cam_Param:
    # Default parameters for Carissma Lucid Triton TRI028S-CC + Universe Optics BL060C (Lucid UC060-5M lens)
    def __init__(self, etime = float(0.03), fov = np.array([80.8, 61.6]), res = np.array([1936, 1464]), pix = np.array([4.5, 4.5]), flen = float(6.0), min_d = float(0.1), height = float(1), wavelength = np.array([460.0e-9, 530.0e-9, 600.0e-9])):
        self.etime = etime              # in s
        self.fov = np.deg2rad(fov)      # [hFov, vFov]
        self.res = res                  # [hRes, vRes]
        self.pix = pix                  # Pixel size in um [H V]
        self.flen = flen                # Focal length in mm
        self.min_d = min_d              # Min Focus Distance in m
        self.height = height            # Height of sensor above ground
        self.wavelength = wavelength    # Best efficiency of Wavelength of (B G R)    

'''
        # Camera Parameters - based on Lucid triton TRI054-S, Kowa LM3JC10M lens
        self.etime = float(0.003)                # in s
        self.fov = np.array([100.2, 83.7])      # [hFoV, vFoV]
        self.res = np.array([2880, 1860])       # [hRes, vRes]
        self.pix = np.array([3, 3])             # Pixel size in um [H V]
        self.flen = float(3.7)               # Focal length in mm
        self.mind = float(0.1)               # Min focus dllistance in m
'''
class Cam_Sek_Param:
    # Based on Sekonix SF3325 RCCB (NA6062 lens)
    # Wavelength are based on IMX429 - To be Updated!!!
    def __init__(self, etime = float(0.002), fov = np.array([60.0, 38.0]), res = np.array([1920, 1208]), pix = np.array([0, 0]), flen = float(5.49), min_d = float(0.1), height = float(1), wavelength = np.array([460.0e-9, 530.0e-9, 600.0e-9])):
        self.etime = etime              # in s
        self.fov = np.deg2rad(fov)      # [hFov, vFov]
        self.res = res                  # [hRes, vRes]
        self.pix = pix                  # Pixel size in um [H V]
        self.flen = flen                # Focal length in mm
        self.min_d = min_d              # Min Focus Distance in m
        self.height = height            # Height of sensor above ground
        self.wavelength = wavelength    # Best efficiency of Wavelength of (B G R) 


class Fog_Param:
    # Default parameters for fog
    def __init__(self, rate = int(10)):
        self.rate = rate                 # Fog Visibility in meters

    def drops_dist(self, drop_size):
        # DSD based on log normal law
        radius = drop_size*1e6            # Equation uses micrometers

        # Parameters from CE:
        if CE:
            if self.rate == 10:
                N = 6108.486
                sigma = 0.74651233
                ln_rt = 0.484235252
            elif self.rate == 20:
                N = 3559.641
                sigma = 0.777099786
                ln_rt = 0.362382125
            elif self.rate == 30:
                N = 2754.154
                sigma = 0.742155699
                ln_rt = 0.334042672
            elif self.rate == 50:
                N = 3189.711
                sigma = 0.656240841
                ln_rt = 0.118386081

        # Equation based on cm-3, need to convert to m-3
        return (N/(radius*sigma*np.power(2*np.pi,0.5)))*np.exp(-np.power(np.log(radius)-ln_rt,2)/(2*np.power(sigma,2)))*1e6

#########################################################################################################
# Supporting functions
def BGR2Bay(colour_img):
    # Extract colour intensity per pixel based on a CFA pattern
    (h, w) = colour_img.shape[:2]
    (B, G, R) = cv2.split(colour_img)

    bay_img = np.empty((h, w), np.uint8)
    # Pattern of RGGB
    bay_img[0::2, 0::2] = R[0::2, 0::2]
    bay_img[1::2, 0::2] = G[1::2, 0::2]
    bay_img[0::2, 1::2] = G[0::2, 1::2]
    bay_img[1::2, 1::2] = B[1::2, 1::2]

    cv2.imwrite(str(pathlib.Path.cwd() / 'debug_img/02_img_gray.png'), bay_img)

    return bay_img


def att_mask_lum(cam_param):
    # Compute the Luminosity of frame
    # https://www.pixelsham.com/2020/12/26/exposure-value-measurements/
    if CE:
        att_mat = 732       # Overcast/10k Lux
        att_mat = 255
    else:
        att_mat = 600      # Roughly for a sunny day

    return att_mat

def coc(mask):
    # Blur image for circle of confusion
    # Gaussian blur for now
    return cv2.GaussianBlur(mask, [3,3], cv2.BORDER_TRANSPARENT)


def water_ridx(blue = complex(1.3366, 9.86e-10), green = complex(1.3338, 1.448e-9),red = complex(1.3320, 1.09e-8)):
    # Complex number of water
    # Identify based on  wavelength from https://refractiveindex.info/?shelf=3d&book=liquids&page=water
    # lambda(460, 530, 600) nm - THI Lucid triton
    
    refract_idx = np.array([blue, green, red])

    return refract_idx

####################################################################################################################
# Main Functions
def debug_img(loc): 
    # Load Image
    return BGR2Bay(cv2.imread(str(pathlib.Path.cwd() / loc)))

def Q_ext(Diameter, step_size):
    # Based on Qext in A Methodology to Model the Rain and Fog Effect on the Performance of Automotive LiDAR Sensors
    # Converted from Jonathan's c_ext.m code in matlab
    # Calculate Dia +- half step size for use in sigma trapezium integration
    Dia = Diameter - (step_size/2)
    np.append(Dia, Diameter[-1] + ((Diameter[1]-Diameter[0])/2))
    size_param = (np.pi*Dia)/(cam.wavelength)
    ref_idx = water_ridx()
    N_Max = size_param + (4*np.power(size_param, (1/3))) + 10
    ext_c = np.zeros(3)

    # Rical Bessel Functions for Mie Scattering
    # https://www.researchgate.net/publication/41466575_Generating_Bessel_Functions_In_Mie_Scattering_Calculations_Using_Continued_Fractions
    for idx_k in range(3):
        # idx_k is the index for B, G, R respectively
        # Calculate intial N=0 psi and epsilon
        psi_x_prev = size_param[idx_k] * scipy.special.spherical_jn(0, size_param[idx_k])
        psi_mx_prev = ref_idx[idx_k] * size_param[idx_k] * scipy.special.spherical_jn(0, ref_idx[idx_k] * size_param[idx_k])
        epsilon_x_prev = complex(psi_x_prev.real, (size_param[idx_k] * scipy.special.spherical_yn(0, size_param[idx_k]).real))
                          
        for idx_j in range(int(N_Max[idx_k])):
            idx_i = idx_j + 1
            psi_x = size_param[idx_k] * scipy.special.spherical_jn(idx_i, size_param[idx_k])
            psi_mx = ref_idx[idx_k] * size_param[idx_k] * scipy.special.spherical_jn(idx_i, ref_idx[idx_k] * size_param[idx_k])

            epsilon_x = complex(psi_x.real, (size_param[idx_k] * scipy.special.spherical_yn(idx_i, size_param[idx_k]).real))

            rn_mx = psi_mx_prev / psi_mx

            mie_an = ((psi_x*((rn_mx/ref_idx[idx_k])+(idx_i*(1-(1/np.power(ref_idx[idx_k], 2)))/size_param[idx_k]))-psi_x_prev)/(epsilon_x*((rn_mx/ref_idx[idx_k])+(idx_i*(1-(1/np.power(ref_idx[idx_k], 2)))/size_param[idx_k]))-epsilon_x_prev))
            mie_bn = ((rn_mx*ref_idx[idx_k]*psi_x)-psi_x_prev)/((rn_mx*ref_idx[idx_k]*epsilon_x)-epsilon_x_prev)
            mie_ratio = mie_an + mie_bn

            ext_c[idx_k] = ext_c[idx_k] + ((2*idx_i) + 1)*mie_ratio.real

            # Update psi/epsilon previous value
            psi_x_prev = psi_x
            psi_mx_prev = psi_mx
            epsilon_x_prev = epsilon_x

    ext_q = (2/np.power(size_param, 2)) * ext_c
    return ext_q


def Q_ext_BGR(Dia):
    # Based on Qext in A Methodology to Model the Rain and Fog Effect on the Performance of Automotive LiDAR Sensors
    # Converted from Jonathan's c_ext.m code in matlab
    # Calculate Dia +- half step size for use in sigma trapezium integration
    size_param = (np.pi*Dia)/(cam.wavelength)
    ref_idx = water_ridx()

    N_Max = size_param + (4*np.power(size_param, (1/3))) + 10
    ext_c = np.zeros(3)

    # Rical Bessel Functions for Mie Scattering
    # https://www.researchgate.net/publication/41466575_Generating_Bessel_Functions_In_Mie_Scattering_Calculations_Using_Continued_Fractions
    for idx_k in range(3):
        # idx_k is the index for B, G, R respectively
        # Calculate intial N=0 psi and epsilon
        psi_x_prev = size_param[idx_k] * scipy.special.spherical_jn(0, size_param[idx_k])
        psi_mx_prev = ref_idx[idx_k] * size_param[idx_k] * scipy.special.spherical_jn(0, ref_idx[idx_k] * size_param[idx_k])
        epsilon_x_prev = complex(psi_x_prev.real, (size_param[idx_k] * scipy.special.spherical_yn(0, size_param[idx_k]).real))
                          
        for idx_j in range(int(N_Max[idx_k])):
            idx_i = idx_j + 1
            psi_x = size_param[idx_k] * scipy.special.spherical_jn(idx_i, size_param[idx_k])
            psi_mx = ref_idx[idx_k] * size_param[idx_k] * scipy.special.spherical_jn(idx_i, ref_idx[idx_k] * size_param[idx_k])

            epsilon_x = complex(psi_x.real, (size_param[idx_k] * scipy.special.spherical_yn(idx_i, size_param[idx_k]).real))

            rn_mx = psi_mx_prev / psi_mx

            mie_an = ((psi_x*((rn_mx/ref_idx[idx_k])+(idx_i*(1-(1/np.power(ref_idx[idx_k], 2)))/size_param[idx_k]))-psi_x_prev)/(epsilon_x*((rn_mx/ref_idx[idx_k])+(idx_i*(1-(1/np.power(ref_idx[idx_k], 2)))/size_param[idx_k]))-epsilon_x_prev))
            mie_bn = ((rn_mx*ref_idx[idx_k]*psi_x)-psi_x_prev)/((rn_mx*ref_idx[idx_k]*epsilon_x)-epsilon_x_prev)
            mie_ratio = mie_an + mie_bn

            ext_c[idx_k] = ext_c[idx_k] + ((2*idx_i) + 1)*mie_ratio.real

            # Update psi/epsilon previous value
            psi_x_prev = psi_x
            psi_mx_prev = psi_mx
            epsilon_x_prev = epsilon_x

    ext_q = (2/np.power(size_param, 2)) * ext_c
    return ext_q


def ext_co(ext_q, size):
    # integrate step size based on trapezium rule
    # Convert diameter into radius for extinction calculations
    size = size / 2
    step = size[1]-size[0]
    # N(D) extra step as N(D) will be used in another integration equation
    n_drop_size = np.zeros(len(size))               # Number of drops for each size
    sigma_d = np.zeros([len(size)-1])               # Extinction at different drop size
    ext_sigma = np.zeros(3)
    for idx in range(len(size)):
        n_drop_size[idx] = weather.drops_dist(size[idx])
        
    # Calculate the coefficient
    # step = step*1e3             # Convert Step size from m to mm
    for idx_i in range(3):
        # Calculate per colour channel
        for idx in range(len(size)-1):
            #sigma_d[idx] =  step * ((n_drop_size[idx]*ext_q[idx, idx_i]*np.power(size[idx],2))  + n_drop_size[idx+1]*ext_q[idx+1, idx_i]*np.power(size[idx+1],2) / 2)
            sigma_d[idx] =  step * (((n_drop_size[idx]*ext_q[idx, idx_i]*np.power(size[idx],2))  + n_drop_size[idx+1]*ext_q[idx+1, idx_i]*np.power(size[idx+1],2)) / 2)
        ext_sigma[idx_i] = np.sum(sigma_d) * np.pi
    return ext_sigma * 1e6


def sigma_BGR(sigma, res):
    sig_mat = np.empty(shape = res, dtype=np.double)
    sig_mat[0::2, 0::2] = sigma[2]
    sig_mat[1::2, 0::2] = sigma[1]
    sig_mat[0::2, 1::2] = sigma[1]
    sig_mat[1::2, 1::2] = sigma[0]

    return sig_mat


def fog_mask(img, depth, sig, cam):
    transmission = np.exp(-sigma_BGR(sig, np.shape(img))*depth)
    fog_img = transmission * img + (1- transmission)* att_mask_lum(cam)
    return fog_img


def Bay2BGR(bay):
    #BGR_img = cv2.cvtColor(bay, cv2.COLOR_BayerBG2BGR)
    BGR_img = cv2.cvtColor(np.clip(demosaicing_CFA_Bayer_Malvar2004(bay), 0, 255).astype(np.uint8), cv2.COLOR_RGB2BGR)

    return BGR_img

####################################################################################################################
if __name__ == '__main__':
    no_frames = 40
    CE = True
    t_start = time.time()
    # Initialise
    if CE:
        cam = Cam_Sek_Param(etime=0.004)
        weather = Fog_Param(rate=50)
    else:
        cam = Cam_Param(etime=0.002)                                   # Use Cam_Sek_Param for sekonix cam parameters
        weather = Fog_Param(rate=10)
    #weather = Snow_Param(rate=30, dis_range=[1.5, 50.0])                       # Rain/Snow_Param


    if CE:
        img_name = 'CE_03_car02_day04_clean_plate_15m_frame_51'
        #img_name = 'SIM_THI_Car_28m'
        raw_img = debug_img('input_img/'+ img_name + '.png')
    else:
        img_name = '03_car01_day01_clear_weather01_28m_45deg_frame_201'
        raw_img = debug_img('input_img/'+ img_name + '.png')
        
    depth = np.loadtxt('input_img/depth/' + img_name +'.txt')/3                # Estimated depth from txt file

    test_img_gray = raw_img.copy()
    
    # Pre Computing extinction parameter
    num_step = 199
    # num_step = 1981
    dia_step = np.linspace(0.2e-6, 20e-6, num = num_step)
    Qextinction = np.zeros([num_step, 3])
    for idx_x in range(len(Qextinction)):
        Qextinction[idx_x,:] = Q_ext_BGR(dia_step[idx_x]/2)
    sigma_ext = ext_co(Qextinction, dia_step)

    t_av_start = time.time()
    for idx_f in range(no_frames):
        print('Generating frame #' + str(idx_f+1) + ' of #' + str(no_frames) + ' ........................')
        

        noisy_img = fog_mask(test_img_gray, depth, sigma_ext, cam)
        t_end = time.time()
        print('Noise Applciation = ' + str(t_end-t_start))

        noisy_c_img = Bay2BGR(noisy_img)

        
        # if no_frames == 1:
        #     cv2.imwrite((str(pathlib.Path.cwd() / 'debug_img/03_noise_') + str(weather.rate) + '_mmh_'+ str(idx_f) + '.png'), noisy_img)
        #     cv2.imwrite((str(pathlib.Path.cwd() / 'debug_img/03_noise_colour_') + str(weather.rate) + '_mmh_' + str(idx_f) + '.png'), noisy_c_img)
        # else:
        #     cv2.imwrite((str(pathlib.Path.cwd() / 'debug_img/Multi/03_noise_') + str(weather.rate) + '_mmh_'+ str(idx_f) + '.png'), noisy_img)
        #     cv2.imwrite((str(pathlib.Path.cwd() / 'debug_img/Multi/03_noise_colour_') + str(weather.rate) + '_mmh_' + str(idx_f) + '.png'), noisy_c_img)
    t_av_end = time.time()
    print("Average frame time = " + str((t_av_end-t_av_start)/no_frames))