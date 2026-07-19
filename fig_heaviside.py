import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 11, 'figure.dpi': 150})

c0, GM, Z0 = 1.0, 2.0, 377.0
grid_size, resolution = 20, 400
x = np.linspace(-grid_size, grid_size, resolution)
X, Y = np.meshgrid(x, x)

r_rest = np.clip(np.sqrt(X**2 + Y**2), 0.5, None)
Z_rest = Z0 * np.exp(2*GM/(r_rest*c0**2))

v = 0.85*c0; gamma = 1/np.sqrt(1 - v**2/c0**2)
r_mvt = np.clip(np.sqrt(gamma**2*X**2 + Y**2), 0.5, None)
Z_mvt = Z0 * np.exp(2*GM/(r_mvt*c0**2))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2))
levels = np.linspace(Z0, Z0*5, 15)
cp1 = ax1.contourf(X, Y, Z_rest, levels=levels, cmap='magma', extend='max')
ax1.contour(X, Y, Z_rest, levels=levels, colors='white', alpha=0.3, linewidths=0.6)
ax1.plot(0, 0, 'wo', ms=7)
ax1.set_title(r'Soliton at rest ($v=0$): isotropic profile')
ax1.set_xlabel('x'); ax1.set_ylabel('y'); ax1.set_aspect('equal')

cp2 = ax2.contourf(X, Y, Z_mvt, levels=levels, cmap='magma', extend='max')
ax2.contour(X, Y, Z_mvt, levels=levels, colors='white', alpha=0.3, linewidths=0.6)
ax2.plot(0, 0, 'wo', ms=7)
ax2.arrow(-10, 15, 20, 0, head_width=1.8, head_length=2.6, fc='cyan', ec='cyan', width=0.45)
ax2.text(0, 17.3, r'$v = 0.85\,c_0$', color='cyan', ha='center', fontsize=12, fontweight='bold')
ax2.set_title(r'Soliton in motion: Heaviside ellipsoid')
ax2.set_xlabel('x'); ax2.set_aspect('equal')

cbar = fig.colorbar(cp2, ax=[ax1, ax2], fraction=0.025, pad=0.03)
cbar.set_label(r'Local impedance $Z(x,y)\ \ [\Omega]$')
fig.suptitle(r'Boosted solution of the scalar sector: $Z = Z_0\exp\!\left[2GM/c_0^2\sqrt{\gamma^2x^2+y^2+z^2}\right]$', fontsize=12)
fig.savefig('fig_heaviside.pdf', bbox_inches='tight')
fig.savefig('fig_heaviside.png', bbox_inches='tight')
print("Heaviside figure saved.")
