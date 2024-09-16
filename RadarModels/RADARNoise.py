def add_normal_noise(ampNoise, shape):
    noise = np.random.normal(0, ampNoise, shape) + 1j * np.random.normal(0, ampNoise, shape)
    return noise


def add_thermal_noise(signal_Rx_origin, f_start, f_end, P_Tx, T, Z_I):
    k_B = 1.380649E-23  # Boltzmann constant
    PNoise = 4 * Z_I * k_B * T * (f_end - f_start) / (10**(P_Tx / 10))
    ampNoise = np.sqrt(PNoise)
    return signal_Rx_origin + add_normal_noise(ampNoise, signal_Rx_origin.shape)


def rain_attenuation_function(rainrate):
    """
    Rain Attenuation for 90GHz
    J. Sander, "Rain attenuation of millimeter waves at λ = 5.77, 3.3, and 2 mm," in IEEE Transactions on Antennas and Propagation, vol. 23, no. 2, pp. 213-220, March 1975, doi: 10.1109/TAP.1975.1141059.
    keywords: {Rain;Attenuation;Millimeter wave technology;Laboratories;Electromagnetic propagation;Frequency;Mie scattering;Particle scattering;Electromagnetic scattering;Temperature},
    """
    a = 1.07
    b = 0.69
    alpha = 0
    beta = 0
    Y = a * (rainrate ** b) * (1 + alpha * (rainrate ** (-beta)))
    return Y


def fog_atenuation_function(foglevel):
    if foglevel == 10:
        N = 6108.486
        sigma = 0.74651233
        ln_tilde_r = 0.484235252
    elif foglevel == 20:
        N = 3559.641
        sigma = 0.777099786
        ln_tilde_r = 0.362382125
    elif foglevel == 30:
        N = 2754.154
        sigma = 0.742155699
        ln_tilde_r = 0.334042672
    elif foglevel == 50:
        N = 3189.711
        sigma = 0.656240841
        ln_tilde_r = 0.118386081
    else:
        return 0
    α = 0.182
    frequency = 77
    γ = α * calculate_lwc(N, sigma, ln_tilde_r) * frequency ** 2  # dB/km
    return γ


def calculate_lwc(N, sigma, ln_tilde_r):
    """
    Calculate the Liquid Water Content (LWC) based on a log-normal droplet size distribution.

    Parameters:
    N (float): Total number of droplets
    sigma (float): Standard deviation of the logarithmic radius distribution
    ln_tilde_r (float): Mean of the logarithmic radius distribution

    Returns:
    float: Calculated LWC in g/m³
    """
    # Convert ln_tilde_r to tilde_r
    tilde_r = np.exp(ln_tilde_r) * 1e-3

    # Calculate the LWC in m³/m³
    lwc_volume = (4 / 3) * np.pi * N * (tilde_r ** 3) * np.exp(9 * (sigma ** 2) / 2)

    # Convert volume concentration to mass concentration (assuming density of water = 1000 kg/m³ or 1 g/cm³)
    density_water = 1000  # kg/m³ or equivalently 1 g/cm³
    lwc_mass = lwc_volume * density_water  # g/m³

    return lwc_mass


def snow_atenuation_function(snowfall):  # mm/hr
    """
    F. Norouziari et al., "Low-THz Wave Snow Attenuation," 2018 International Conference on Radar (RADAR), Brisbane, QLD, Australia, 2018, pp. 1-4, doi: 10.1109/RADAR.2018.8557275. keywords: {specific attenuation;snowfall;Low-THz wave},
    attenuation through dry snow
    """
    return 0.8205 * snowfall  # dB/km
