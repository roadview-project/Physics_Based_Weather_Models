import numpy as np
import cv2
import pathlib
import time
import imgaug.augmenters as iaa
import multiprocessing

from colour_demosaicing import (
    demosaicing_CFA_Bayer_Malvar2004,
)


##########################################################################################################
## Initialise Classes for camera and weather
class Cam_Param:
    # Default parameters for Carissma Lucid Triton TRI028S-CC + Universe Optics BL060C (Lucid UC060-5M lens)
    def __init__(self, etime = float(0.03), fov = np.array([80.8, 61.6]), res = np.array([1936, 1464]), pix = np.array([4.5, 4.5]), flen = float(6.0), f_dot = float(2.1), min_d = float(0.1), height = float(1)):
        self.etime = etime              # in s
        self.fov = np.deg2rad(fov)      # [hFov, vFov]
        self.res = res                  # [hRes, vRes]
        self.pix = pix                  # Pixel size in um [H V]
        self.flen = flen                # Focal length in mm
        self.f_dot = f_dot              # Aperture
        self.min_d = min_d              # Min Focus Distance in m
        self.height = height            # Height of sensor above ground

'''
        # Camera Parameters - based on Lucid triton TRI054-S, Kowa LM3JC10M lens
        self.etime = float(0.003)                # in s
        self.fov = np.array([100.2, 83.7])      # [hFoV, vFoV]
        self.res = np.array([2880, 1860])       # [hRes, vRes]
        self.pix = np.array([3, 3])             # Pixel size in um [H V]
        self.flen = float(3.7)               # Focal length in mm
        self.mind = float(0.1)               # Min focus distance in m
'''
class Cam_Sek_Param:
    # Based on Sekonix SF3325 RCCB (NA6062 lens)
    def __init__(self, etime = float(0.002), fov = np.array([60.0, 38.0]), res = np.array([1920, 1208]), pix = np.array([0, 0]), flen = float(5.49), min_d = float(0.1), height = float(1)):
        self.etime = etime              # in s
        self.fov = np.deg2rad(fov)      # [hFov, vFov]
        self.res = res                  # [hRes, vRes]
        self.pix = pix                  # Pixel size in um [H V]
        self.flen = flen                # Focal length in mm
        self.min_d = min_d              # Min Focus Distance in m
        self.height = height            # Height of sensor above ground

class Cam_ACDC_Param:
    # Parameters based on GOPro hero 5 used in ACDC dataset
    def __init__(self, etime = float(0.03), fov = np.array([118.2, 69.5]), res = np.array([1920, 1080]), pix = np.array([4.5, 4.5]), flen = float(14.0), f_dot = float(2.1), min_d = float(0.1), height = float(1)):
        self.etime = etime              # in s
        self.fov = np.deg2rad(fov)      # [hFov, vFov]
        self.res = res                  # [hRes, vRes]
        self.pix = pix                  # Pixel size in um [H V]
        self.flen = flen                # Focal length in mm
        self.f_dot = f_dot              # Aperture
        self.min_d = min_d              # Min Focus Distance in m
        self.height = height            # Height of sensor above ground

        sensor_size = np.array([6.17, 4.55])    # Sensor size for GoPro (not needed!)
        self.pix = self.res/sensor_size



class Rain_Param:
    # Default parameters for rain
    def __init__(self, rate = float(10), dia_range = np.array([0.5, 6.0]), dia_step_size = float(0.1), dis_range = np.array([1.00, 30.00]), dis_step_size = float(0.01)):
        self.rate = rate                    # Rainrate in mm/h
        self.dia_range = dia_range          # Range of raindrop size [min max], in mm
        self.dia_step_size = dia_step_size  # Step size for raindrop in lookup table
        self.dis_range = dis_range          # Range of distances to spawn raindrops [min max], in m
        self.dis_step_size = dis_step_size  # Step size for distance in lookup table


        self.dia_step_num = int(1 + ((dia_range[1] - dia_range[0]) / dia_step_size))    # Number of steps for linspace in lookup table
        self.dis_step_num = int(1 + ((dis_range[1] - dis_range[0]) / dis_step_size))    # Number of steps for linspace in lookup table

    def velocity(self, Diameter):
        # Based on empirical power law 
        # https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8936913 or https://journals.ametsoc.org/view/journals/apme/16/12/1520-0450_1977_016_1322_paairm_2_0_co_2.xml?tab_body=pdf 
        return 3.78 * (Diameter**0.67)
    
    def drops_dist(self, drop_size):
        return 8000 * np.exp(-4.1*np.power(self.rate, -0.21) * drop_size)

    def attenuation(self, img):
        # Create attenuation matrix based on image to apply noise on right colour channel
        # return 255
        return ((0.94 * att_mask_bayer(img)) + (0.06 * att_mask_lum(img)))


class Snow_Param:
    # Default parameters for snow
    # Use same variables as rain so the framework code can be used as it is
    # Assume snow are spherical initially - to be updated to oval
    def __init__(self, rate = float(10), dia_range = np.array([0.1, 10]), dia_step_size = float(0.1), dis_range = np.array([0.50, 30.09]), dis_step_size = float(0.01)):
        self.rate = rate                    # Rainrate in mm/h
        self.dia_range = dia_range          # Range of raindrop size [min max], in mm
        self.dia_step_size = dia_step_size  # Step size for raindrop in lookup table
        self.dis_range = dis_range          # Range of distances to spawn raindrops [min max], in m
        self.dis_step_size = dis_step_size  # Step size for distance in lookup table


        self.dia_step_num = int(1 + ((dia_range[1] - dia_range[0]) / dia_step_size))    # Number of steps for linspace in lookup table
        self.dis_step_num = int(1 + ((dis_range[1] - dis_range[0]) / dis_step_size))    # Number of steps for linspace in lookup table

    def velocity(self, Diameter):
        # Velocity based on "Fast Reactive Control for Illumination Through Rain and Snow"
        velocity = 0.65 - 10.*np.exp(-600 * Diameter)

        # Another accounting for snow properties can be found eq.9 https://journals.ametsoc.org/view/journals/apme/46/5/jam2489.1.xml

        return velocity
    
    def drops_dist(self, drop_size):
        # DSD based on "Fast Reactive Control for Illumination Through Rain and Snow"
        #n_drop = 3800*np.power(self.rate,-0.87)*np.exp(-25.5*np.power(self.rate,-0.48))

        # DSD based on A Statistical and Physical Description of Hydrometeor Distributions in Colorado Snowstorms Using a Video Disdrometer
        # Table 1 - Snowfall liquid equivanlent rate of 3.03 mm/h, using the Gamma model
        if self.rate == 3.03:
            N_O = 2.33e4
            mu = 1.23
            slope_lambda = 2.25
            gamma = True
        elif self.rate ==2.70:
            N_O = 1.39e3
            mu = -0.90
            slope_lambda = 1.40
            gamma = True
        elif self.rate ==2.59:
            N_O = 7.18e1
            mu = -0.78
            slope_lambda = 0.35
            gamma = True
        else:
            # Generalised equation for snow - Does not work well!!!
            # Uses liquid water content which is equivalent to the rainfall interpretation of snow
            rho_snow = 0.178 * np.power(drop_size, -0.992)
            D_O = 1e-3
            R_LWC = np.sqrt(np.power(self.rate/(487*rho_snow*D_O*self.velocity(drop_size)), 3))
            R_LWC = np.power(np.power(self.rate/(487*rho_snow*D_O*self.velocity(drop_size)), 3), 0.5)
            N_O = 7600 * np.power(R_LWC,-0.87)
            slope_lambda = 2.55 * np.power(R_LWC,-0.48)
            gamma = False
        if gamma:
            n_drop = N_O * np.power(drop_size, mu) * np.exp(-slope_lambda * drop_size)
        else:
            n_drop = N_O*np.exp(slope_lambda*drop_size)

        return n_drop
    
    def attenuation(self, img):
        # Create attenuation matrix based on image to apply noise on right colour channel
        # Snow is white
        # https://www.cambridge.org/core/journals/journal-of-glaciology/article/spectral-reflectances-of-snow-and-freshwater-ice-from-340-through-1-100-nm/EA46E8A87DADFBD874DC7986EE774236
        # Spectral Reflectivity roughly 0.85, based on New snow measured, will reflect global luminance
        # Assumes no localised sources of light
        # return 0.85 * att_mask_lum(img)
        intensity = 255
        # intensity = np.clip(int(np.random.default_rng().normal(225,15)), 200, 255)
        return intensity
    

#########################################################################################################
# Functions for the main functions
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


def att_mask_bayer(img):
    # For rain, average colour of the frame
    # Crop frame due to the limited rain area and enclosed
    if CE:
        # Reduce side wall colour (taking central half) at CE
        cr_left = int((img.shape[1])/4)
        cr_right = int((img.shape[1]) - ((img.shape[1])/4))
    else:
        cr_left = 0
        cr_right = img.shape[1]
    # Reduce the image by ignoring the bottom 2/5 of the image
    cr_bot = int(img.shape[0] - (2*(img.shape[0])/5))
    img = img[0:cr_bot, cr_left:cr_right]


    # Compute the average colour in of each channel for the frame
    av_r = int(img[0::2, 0::2].mean())
    av_g = int((img[0::2,1::2].mean() + img[1::2, 0::2].mean()) / 2)
    av_b = int(img[1::2, 1::2].mean())

    # Build the attenuation matrix
    att_mat = np.empty(shape = np.shape(raw_img), dtype=np.double)
    att_mat[0::2, 0::2] = av_r
    att_mat[1::2, 0::2] = av_g
    att_mat[0::2, 1::2] = av_g
    att_mat[1::2, 1::2] = av_b

    return att_mat


def att_mask_lum(img):
    # Comupte the Luminosity of frame
    # https://www.pixelsham.com/2020/12/26/exposure-value-measurements/
    att_mat = 6000      # Roughly for a sunny day
    # att_mat = 737       # Heavy Overcast

    return att_mat


def coc(mask):
    # Blur image for circle of confusion
    # Gaussian blur for now
    return cv2.GaussianBlur(mask, [3,3], cv2.BORDER_TRANSPARENT)


def contrast(img, level):
    # Level should be between 1 to 5?
    # Can also use pillow's PIL.ImageEnhance.Contrast
    return iaa.imgcorruptlike.apply_contrast(img, severity=level)


#####################################################################################################
# Main Functions
def lookup_tbl2():
    # Create the lookup tables
    rain_dia = np.linspace(weather.dia_range[0], weather.dia_range[1], num = weather.dia_step_num)                      # in mm array of sizes to test
    rain_dis = np.linspace(weather.dis_range[0], weather.dis_range[1], num = weather.dis_step_num)

    # Initialise look-up tables
    tbl_stk = np.zeros(shape = [rain_dis.shape[0], rain_dia.shape[0]])
    tbl_hor = tbl_stk.copy()
    tbl_time = tbl_stk.copy()

    
    for idx_i, dia in enumerate(rain_dia):
        # Set parameters per Diameter of rain
        # Speed of Raindrop
        vel_rain = weather.velocity(Diameter = dia)

        for idx_j, dis in enumerate(rain_dis):
            pix_height = ((2 * dis * np.tan(cam.fov[1]/2)) / cam.res[1])
            
            num_pix_vert = dia*1e-3 / pix_height
            vel_rain_pix = vel_rain / pix_height
            
            tbl_stk[idx_j, idx_i] = (vel_rain_pix * cam.etime) + num_pix_vert
            tbl_hor[idx_j, idx_i] = dia*1e-3 / ((2 * dis * np.tan(cam.fov[0]/2)) / cam.res[0])
            if tbl_stk[idx_j, idx_i] < 1:
                tbl_time[idx_j, idx_i] = cam.etime
            elif tbl_stk[idx_j, idx_i] < num_pix_vert + 1:
                tbl_time[idx_j, idx_i] = tbl_stk[idx_j, idx_i]/vel_rain_pix
            # if tbl_stk[idx_j,idx_i] < num_pix_vert * 1.5:
                # tbl_time[idx_j, idx_i] = cam.etime
            else:
                tbl_time[idx_j, idx_i] = (1 + num_pix_vert)/vel_rain_pix

    return tbl_stk, tbl_hor, tbl_time, rain_dia


def input_frame():
    # 

    return 


def debug_img(loc):
    # Load Image
    return BGR2Bay(cv2.imread(str(pathlib.Path.cwd() / loc)))



def integration_trap(size):
    # Integrate each step size for number of drops
    step = weather.dia_step_size
    n_drop_size = np.zeros(len(size))           # Drops array for all size
    n_drop = int(0)                             # Total drop number
    for idx in range(len(size)):
        n_drop_size[idx] = step * ((weather.drops_dist(drop_size = (size[idx]-(step/2))) + weather.drops_dist(drop_size = (size[idx]+(step/2)))) / 2)
        n_drop = n_drop + n_drop_size[idx]
    drop_ratio = np.asarray(n_drop_size / n_drop)
    return n_drop, drop_ratio


def fov_cuboid(range):
    # Compute Volume to create rain drop
    # Assumes ground flat and parallel

    vol = (2*range*np.tan(cam.fov[0]/2)) * ((range*np.tan(cam.fov[1]/2))+cam.height) * range 

    return vol


def poission(num_drop):
    poi_drops = np.random.poisson(num_drop, 1)

    return poi_drops


def rnd_loc_cuboid(num_rain_size, max_d):
    #Random location generation for cuboid volume
    rain_drops = np.sum(num_rain_size, dtype=int)
    rain_loc = np.zeros(shape= [rain_drops, 4], dtype=np.double)                           # Stores pointers to steps in lookup table
    half_width = max_d*np.tan(cam.fov[0]/2) * 1000                                       # At max range in mm half due to centerpoint
    half_height = max_d*np.tan(cam.fov[1]/2) * 1000                                       # At max range in mm half due to centerpoint
    idx_j = int(0)
    for idx_i in range(len(num_rain_size)):
        rain_loc[idx_j:(idx_j+num_rain_size[idx_i]), 0] = np.full((num_rain_size[idx_i]), idx_i)                                                            # Drop size pointer
        rain_loc[idx_j:(idx_j+num_rain_size[idx_i]), 1] = (np.random.default_rng().integers(-half_width, half_width, size=(num_rain_size[idx_i]))) / 1000   # Horizontal location Back to meters
        rain_loc[idx_j:(idx_j+num_rain_size[idx_i]), 2] = (np.random.default_rng().integers(-cam.height*1000, half_height, size=(num_rain_size[idx_i]))) / 1000  # Vertical starting locationback to meters
        rain_loc[idx_j:(idx_j+num_rain_size[idx_i]), 3] = np.random.default_rng().integers(0, weather.dis_step_num, size=(num_rain_size[idx_i]))                                 # uniform distribution for "Parallel" distance
        idx_j = idx_j + num_rain_size[idx_i]
    return rain_loc


def frame2pix(rain_loc):
    # Rain_Loc in [Drop size pointer, H_loc, V_loc, Distance pointer], please use fov_cuboid for this
    # Need to be updated if pixels are not even!
    # Takes rain coordinates frame and convert them into pixels. Assumes rain frame aligned to cam frame. Principle point in the centre of image
    # Convert rain pointer to real distance with hard coding
    distance = weather.dis_range[0] + (rain_loc[:, 3] * weather.dis_step_size)

    w_factor = (2 * np.tan(cam.fov[0]/2)) / cam.res[0]                  # Horizonal length per pixel based on distance
    v_factor = (2 * np.tan(cam.fov[1]/2)) / cam.res[1]                  # Vertical length per pixel based on distance
    w_pix = rain_loc[:,1] / (distance * w_factor)                      # Location horizontally based on principle point
    v_pix = rain_loc[:,2] / (distance * v_factor)                      # Location vertically based on principle point
    
    # Flip y axis to match camera pixel
    v_pix = np.negative(v_pix)

    # Find the principle point which is at the centre
    w_principle = cam.res[0] / 2
    v_principle = cam.res[1] / 2

    # Convert principle pixel to Image pixel, round down for pixel location, does not matter which part of pixel it is on
    rain_loc[:,1] = np.floor(w_pix + w_principle)
    rain_loc[:,2] = np.floor(v_pix + v_principle)
    
    # Remove drops outisde of H FoV and then V FoV
    rain_loc = np.delete(rain_loc, np.nonzero(np.bitwise_or(rain_loc[:,1]<0, rain_loc[:,1]>=cam.res[0])), 0)
    rain_loc = np.delete(rain_loc, np.nonzero(np.bitwise_or(rain_loc[:,2]<0, rain_loc[:,2]>=cam.res[1])), 0)

    rain_loc = rain_loc.astype(int)

    return rain_loc


def attenuation_mask(raw_img, rain_loc, tbl_stk, tbl_hor, tbl_time, depth_map):
    # Create a mask based on raindrop to apply noise
    # Assumes rain at the centre of the pixel
    # Initialise the mask
    (h, w) = raw_img.shape[:2]
    att_mask = np.zeros((h, w), np.double)
    
    # Create attenuation matrix based on image to apply noise on right colour channel
    att_mat = weather.attenuation(img=raw_img)

    filter_perc = 0.0000                                                      # Percent of drop in decimals covered to ignore
    for idx_i in range((np.shape(rain_loc))[0]):
        # Apply every drop into a mask based on total drops time on a pixel

        # Compare against estimated distance
        if rain_loc[idx_i, 3] > depth_map[rain_loc[idx_i, 2], rain_loc[idx_i, 1]]:
            # Skip if drop distance is more than depth. Depth map is based on "parallel" distance
            continue

        hor = (tbl_hor[rain_loc[idx_i, 3], rain_loc[idx_i, 0]])/ 2              # Horizontal pixels to apply noise left and right of rain location
        vert = tbl_stk[rain_loc[idx_i, 3], rain_loc[idx_i, 0]]                  # Vertical Streak Length
        rel_time = tbl_time[rain_loc[idx_i, 3], rain_loc[idx_i, 0]]/cam.etime   # Relative time on pixel

        if rel_time < filter_perc:                                               # Ignores drops which filled too little of a pixel
            continue

        # Find whole pixels and the remainder
        hor_int = int(np.floor(hor))
        vert_int = int(np.floor(vert))

        hor_rem = hor - hor_int
        vert_rem = vert - vert_int
        
        # Assume centered around left edge
        if vert <= 1:
            # Single pixel droplet
            # If vertical ratio covered is <1, hor ratio should also be <1
            rel_time_percent = rel_time * vert_rem * hor_rem
            if rel_time_percent <= filter_perc:
                continue
            else:
                att_mask[rain_loc[idx_i, 2], rain_loc[idx_i, 1]] = att_mask[rain_loc[idx_i,2], rain_loc[idx_i,1]] + rel_time_percent
        elif hor <= 1:
            # Single column droplet
            # approximate vertical always fully filled
            rel_time_percent = rel_time * hor_rem
            if rel_time_percent <= filter_perc:
                continue
            else:
                att_mask[rain_loc[idx_i,2]:np.clip(rain_loc[idx_i,2]+vert_int, 0, cam.res[1]-1), rain_loc[idx_i, 1]] = att_mask[rain_loc[idx_i,2]:np.clip(rain_loc[idx_i,2]+vert_int, 0, cam.res[1]-1), rain_loc[idx_i,1]] + rel_time_percent
        else:
            # Multi column droplet - apply left and right columns separate to middle
            # Left right has percent of time
            
            rel_time_percent = rel_time * hor_rem / 2       # Equally split for left and right remainder
            if rel_time_percent > filter_perc:
                # Left column
                att_mask[rain_loc[idx_i,2]:np.clip(rain_loc[idx_i,2]+vert_int, 0, cam.res[1]-1), np.clip(rain_loc[idx_i, 1] - hor_int - 1, 0, cam.res[0]-1)] = att_mask[rain_loc[idx_i,2]:np.clip(rain_loc[idx_i,2]+vert_int, 0, cam.res[1]-1), np.clip(rain_loc[idx_i, 1] - hor_int - 1, 0, cam.res[0]-1)] + rel_time_percent

                # Right column
                att_mask[rain_loc[idx_i,2]:np.clip(rain_loc[idx_i,2]+vert_int, 0, cam.res[1]-1), np.clip(rain_loc[idx_i, 1] + hor_int, 0, cam.res[0]-1)] = att_mask[rain_loc[idx_i,2]:np.clip(rain_loc[idx_i,2]+vert_int, 0, cam.res[1]-1), np.clip(rain_loc[idx_i, 1] + hor_int, 0, cam.res[0]-1)] + rel_time_percent
            '''
            # Apply to the central block
            att_mask[rain_loc[idx_i,2]:np.clip(rain_loc[idx_i,2]+vert_int, 0, cam.res[1]-1), np.clip(rain_loc[idx_i, 1] - hor_int, 0, cam.res[0]-1):np.clip(rain_loc[idx_i, 1] + hor_int, 0, cam.res[0]-1)] = att_mask[rain_loc[idx_i,2]:np.clip(rain_loc[idx_i,2]+vert_int, 0, cam.res[1]-1), np.clip(rain_loc[idx_i, 1] - hor_int, 0, cam.res[0]-1):np.clip(rain_loc[idx_i, 1] + hor_int, 0, cam.res[0]-1)] + rel_time
            '''
            # Make Circles for the top and bottomonly when it is more than 2 pixels wide?
            tmp_mask = np.zeros((h, w), np.uint8)
            if hor > 2:
                radius = hor_int - 1
                tmp_mask = cv2.ellipse(tmp_mask, (rain_loc[idx_i,1], rain_loc[idx_i,2]), (radius,radius), 0, 0, 360, 1, -1)
                tmp_mask = cv2.ellipse(tmp_mask, (rain_loc[idx_i,1], rain_loc[idx_i,2]+vert_int), (radius,radius), 0, 0, 360, 1, -1)
            # fill the middle
            tmp_mask[rain_loc[idx_i,2]:np.clip(rain_loc[idx_i,2]+vert_int, 0, cam.res[1]-1), np.clip(rain_loc[idx_i, 1] - hor_int, 0, cam.res[0]-1):np.clip(rain_loc[idx_i, 1] + hor_int, 0, cam.res[0]-1)] = 1

            # Apply circle into mask
            att_mask = att_mask + (tmp_mask*rel_time)

    # Apply the Atteunation
    # Blur
    att_mask = np.clip(att_mask,0,1)
    att_mask = coc(att_mask)
    #raw_img = contrast(raw_img, 1)

    #att_mask = np.clip(att_mask, 0, 1)
    raw_img = np.clip((att_mask * att_mat) + ((1- att_mask) * raw_img), 0, 255)

    # Cast array from float due to mask back to uint8
    raw_img = raw_img.astype(np.uint8)

    # Expand and save the mask image for Debugging
    att_mask = np.ceil(att_mask * 255)
    att_mask = att_mask.astype(np.uint8)
    cv2.imwrite(str(pathlib.Path.cwd() / 'debug_img/04_att_mask.png'), att_mask)

    return raw_img


def Bay2BGR(bay):
    #BGR_img = cv2.cvtColor(bay, cv2.COLOR_BayerBG2BGR)
    BGR_img = cv2.cvtColor(np.clip(demosaicing_CFA_Bayer_Malvar2004(bay), 0, 255).astype(np.uint8), cv2.COLOR_RGB2BGR)

    return BGR_img



if __name__ == '__main__':
    no_frames = 1
    CE = False                                           # Define global variable for if at Cerema facility or not
    t_start = time.time()
    # Initialise
    if CE:
        cam = Cam_Sek_Param(etime=0.004)
        weather = Rain_Param(rate=83, dis_range=[1.5, 30.0])
    else:
        cam = Cam_Param(etime=0.005)                                   # Use Cam_Sek_Param for sekonix cam parameters
        # weather = Snow_Param(rate=3.03, dis_range=[1.5, 30.0])
        weather = Rain_Param(rate=50, dis_range=[1.5, 80.0])
    weather = Snow_Param(rate=30, dis_range=[0.5, 10.0])                       # Rain/Snow_Param
        
    tbl_stk, tbl_hor, tbl_time, dsize = lookup_tbl2()
    t_end = time.time()
    print('Lookup-table = ' + str(t_end-t_start))


    if CE:
        img_name = 'CE_03_car02_day04_clean_plate_15m_frame_51'
        #img_name = 'THI_real'
        raw_img = debug_img('input_img/'+ img_name + '.png')
    else:
        # img_name = 'Car_56m_45deg_frame_051'
        img_name = 'SIM_THI_Car_28m'
        raw_img = debug_img('input_img/'+ img_name + '.png')
    depth = np.loadtxt('input_img/depth/' + img_name +'.txt')                 # Estimated depth from txt file

    test_img_gray = raw_img.copy()
    number_particle_unit, drop_ratio = integration_trap(dsize)            # Drops per unit volume and ratio of drops in certain sizes
    number_particle = number_particle_unit * fov_cuboid(weather.dis_range[1])


    t_end = time.time()
    print('Particle per unit volume calc (pre-processing complete) = ' + str(t_end-t_start))
    
    for idx_f in range(no_frames):
        print('Generating frame #' + str(idx_f+1) + ' of #' + str(no_frames) + ' ........................')
        poi_particle = poission(number_particle)
        print('Raindrops in Cuboid = ' + str(poi_particle))
        num_particle_size = poi_particle * drop_ratio
        num_particle_size = num_particle_size.astype(int)
        particle_loc = rnd_loc_cuboid(num_particle_size, weather.dis_range[1])
        particle_loc = frame2pix(particle_loc)
        print('Particles in frame = ' + str(len(particle_loc[:,0])))
        # rain_loc = rnd_distribution(number_rain)   # Create a random location function!

        t_end = time.time()
        print('Particle drop generation with location = ' + str(t_end-t_start))
        noisy_img = attenuation_mask(raw_img, particle_loc, tbl_stk, tbl_hor, tbl_time, depth)
        t_end = time.time()
        print('Noise Applciation = ' + str(t_end-t_start))

        noisy_c_img = Bay2BGR(noisy_img)


        # cv2.imshow("Grey noise",test_img_gray)
        if no_frames == 1:
            cv2.imwrite((str(pathlib.Path.cwd() / 'debug_img/03_noise_') + str(weather.rate) + '_mmh_'+ str(idx_f) + '.png'), noisy_img)
            cv2.imwrite((str(pathlib.Path.cwd() / 'debug_img/03_noise_colour_') + str(weather.rate) + '_mmh_' + str(idx_f) + '.png'), noisy_c_img)
        else:
            cv2.imwrite((str(pathlib.Path.cwd() / 'debug_img/Multi/03_noise_') + str(weather.rate) + '_mmh_'+ str(idx_f) + '.png'), noisy_img)
            cv2.imwrite((str(pathlib.Path.cwd() / 'debug_img/Multi/03_noise_colour_') + str(weather.rate) + '_mmh_' + str(idx_f) + '.png'), noisy_c_img)
