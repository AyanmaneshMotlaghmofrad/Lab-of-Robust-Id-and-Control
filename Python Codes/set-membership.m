# Create the figure with the corrected label for y_tilde
fig, ax = plt.subplots(2, 1, figsize=(10, 8))

# Input Signal Plot
ax[0].plot(k, u_k, marker='o', color='cornflowerblue', label='u(k)', linestyle='-', linewidth=2, markersize=6)
ax[0].set_title('Input Signal u(k)', fontsize=14, weight='bold')
ax[0].set_xlabel('k (Time)', fontsize=12)
ax[0].set_ylabel('u(k)', fontsize=12)
ax[0].grid(True)
ax[0].legend(fontsize=12)

# Output Signal Plot
ax[1].plot(k, y_k, marker='o', color='tomato', label='y(k)', linestyle='-', linewidth=2, markersize=6)
ax[1].plot(k, y_k_upper, linestyle='--', color='mediumseagreen', label='y(k) + Δη', linewidth=1.5)
ax[1].plot(k, y_k_lower, linestyle='--', color='mediumseagreen', label='y(k) - Δη', linewidth=1.5)
ax[1].fill_between(k, y_k_lower, y_k_upper, color='mediumseagreen', alpha=0.3, label='Bounds')
ax[1].set_title('Output Signal y(k) with Bounds', fontsize=14, weight='bold')
ax[1].set_xlabel('y_{tilde}(k)', fontsize=12)  # Changed horizontal label
ax[1].set_ylabel('y(k)', fontsize=12)
ax[1].grid(True)
ax[1].legend(fontsize=12)

# Adjust layout and save figure
plt.tight_layout()
plt.savefig('/mnt/data/input_output_signals.png', dpi=300)
plt.show()
