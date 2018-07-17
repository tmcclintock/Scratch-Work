import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", size=30, family='serif')

labels = [r"$\Omega_c$", r"$\Omega_\Lambda$", r"$\Omega_b$"]
labels = ["Dark matter", "Dark energy", "Baryonic matter"]
sizes = [25, 70, 5]
colors = ['r','cyan','yellow']

fig, ax = plt.subplots()#figsize=(3,3))
ax.pie(sizes, labels=labels, colors=colors, autopct='%1.0f\%%')
ax.axis('equal')
ax.set_title("Energy density fractions")

fig.savefig("cosmo_pie.png", bbox_inches='tight', dpi=300, transparent=True)
plt.show()
