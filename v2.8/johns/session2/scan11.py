import numpy as np
from defect import spectrum
print(f"{'kz':>5} {'om_lie':>9} {'v=om/kz':>8} {'poids':>7} {'2e mode':>9} {'poids2':>7}")
rows=[]
for kz in (0.4,0.8,1.2,1.6,2.0,2.4,2.8):
    om,wt,_ = spectrum(11,11,kz,"D3")
    sel = (om>0.02)&(om<kz/2-1e-6)
    order = np.argsort(-(wt*sel))
    i,j = order[0], order[1]
    rows.append((kz,om[i],wt[i]))
    print(f"{kz:5.2f} {om[i]:9.4f} {om[i]/kz:8.4f} {wt[i]:7.3f} {om[j]:9.4f} {wt[j]:7.3f}")
np.save("disp11.npy", np.array(rows))
