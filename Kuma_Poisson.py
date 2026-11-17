import numpy as np
import matplotlib.pyplot as plt

# --- Settings ---
MODEL = "kuma"           # Choose "kuma" or "poisson"
P_input = 1.3e-3         # Pickup-cell pressure [mbar]
P_max = 3.0e-3           # Maximum plotting pressure [mbar]

# Pressures must already be calibrated and background-subtracted.

# --- Experimental parameters ---
N = 9678                # Initial number of He atoms
rho = 2.18e28           # He number density [atoms/m³]
T = 300.0               # Kr gas temperature [K]
v = 345.0               # Droplet velocity [m/s]
L = 0.018               # Pickup-cell length [m]

# --- Constants ---
kB = 1.380649e-23             # Boltzmann constant [J/K]
m = 84 * 1.66053906892e-27    # Kr atomic mass [kg]

# Energy conversion factor: 1 cm^-1 = 1.9864e-23 J.
J_PER_CM_INV = 6.62607015e-34 * 299792458 * 100

# --- Relative velocity: same RMS prescription as in the SM ---
v_dop = np.sqrt(3 * kB * T / m)
v_rel = np.sqrt(v**2 + v_dop**2)

# --- Evaporation energy balance ---
# Approximate kinetic energy deposited per captured Kr atom,
# converted from joules to cm^-1.
E_kin = (1.5 * kB * T + 0.5 * m * v**2) / J_PER_CM_INV

E_HeKr = 217.0          # Assumed He–Kr binding energy [cm^-1]
E_KrKr = 932.0          # Assumed bulk-like Kr binding energy [cm^-1]
E_evap = 5.0            # Energy per evaporated He atom [cm^-1]

# Constant number of He atoms evaporated per captured Kr atom.
c = (E_kin + E_HeKr + E_KrKr) / E_evap

# --- Initial droplet radius and pickup coefficient ---
R = (3 * N / (4 * np.pi * rho))**(1 / 3)

# X has units of Pa^-1.
X = np.pi * R**2 * L / (kB * T) * (v_rel / v)

# Complete evaporation occurs when 1 - c*X*P/(3*N) = 0.
# At this pressure, the mean pickup reaches N/c and all initial
# He atoms have evaporated. Do not extrapolate Kuma beyond it.
# Divide by 100 to convert pressure from Pa to mbar.
P_limit = 3 * N / (c * X) / 100

# --- Selected model ---
MODEL = MODEL.strip().lower()

def pickup(P_mbar):
    # Accept either one pressure or an array of pressures.
    P = np.asarray(P_mbar, dtype=float)

    if np.any(~np.isfinite(P)) or np.any(P < 0):
        raise ValueError("Pressure must be finite and non-negative.")

    # Mean Poisson pickup, using pressure in Pa.
    mean_poisson = X * P * 100

    if MODEL == "poisson":
        return mean_poisson

    if MODEL == "kuma":
        if np.any(P > P_limit):
            raise ValueError(
                f"Kuma is valid only up to {P_limit:.3e} mbar."
            )

        return N / c * (
            1 - (1 - c * mean_poisson / (3 * N))**3
        )

    raise ValueError('Choose MODEL = "kuma" or "poisson".')


# --- Result ---
print(f"Model: {MODEL}")
print(f"Kr RMS thermal velocity: {v_dop:.2f} m/s")
print(f"Relative velocity: {v_rel:.2f} m/s")
print(f"Kinetic energy: {E_kin:.1f} cm^-1")
print(f"He atoms evaporated per captured Kr atom: {c:.1f}")
print(f"Kuma complete-evaporation limit: {P_limit:.3e} mbar")
print(f"Mean pickup at {P_input:.3e} mbar: {pickup(P_input):.2f}")

# --- Plot ---
if not np.isfinite(P_max) or P_max <= 0:
    raise ValueError("P_max must be finite and positive.")

# Restrict the Kuma curve to its physical domain.
P_plot_max = min(P_max, P_limit) if MODEL == "kuma" else P_max

P = np.linspace(0, P_plot_max, 500)

plt.figure(figsize=(6, 4))
plt.plot(P, pickup(P), label=MODEL.capitalize())

if P_input <= P_plot_max:
    plt.plot(P_input, pickup(P_input), "ko")

plt.xlabel("Kr pickup pressure (mbar)")
plt.ylabel("Mean number of captured Kr atoms")
plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
plt.xlim(0, P_plot_max)
plt.ylim(bottom=0)
plt.legend(frameon=False)
plt.tight_layout()
plt.show()