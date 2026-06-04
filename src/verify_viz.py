"""
Verification visualizations for Messier Marathon planning.
Loads existing results.json and generates plots to verify:
  1. RA-Dec sky map (Mercator) colored by visibility
  2. Altitude time-heatmap (Gantt-style)
  3. Pareto frontier (surface brightness vs altitude)
  4. Suggested observing timeline
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.patches as mpatches

# CJK font
for font_name in ['PingFang SC', 'Heiti SC', 'STHeiti', 'Arial Unicode MS']:
    try:
        matplotlib.font_manager.findfont(font_name, fallback_to_default=False)
        plt.rcParams['font.family'] = font_name
        break
    except Exception:
        continue

# Load results
with open('output/results.json') as f:
    results = json.load(f)

# Load full catalog for ALL objects (visible or not)
# We need to regenerate full data to show invisible ones
from messier_catalog import MESSIER
from planner import ra_dec_to_deg, compute_surface_brightness, compute_visibility, BEIJING, OBJ_ALT_MIN, SUN_ALT_LIMIT

# Get visibility computation (reuse from planner)
all_results, times_utc = compute_visibility()

# Build lookup
visible_m_nums = {r['m_num'] for r in results}
all_by_m = {r['m_num']: r for r in all_results}

UTC_OFFSET = 8


def plot_sky_map(output_file='output/verify_sky_map.png'):
    """Mercator-style sky map: RA vs Dec, colored by max altitude.
    Includes ALL Messier objects. Shows Beijing's latitude line and visibility zone."""
    fig, ax = plt.subplots(figsize=(16, 9))

    # All Messier objects
    all_ra, all_dec, all_alt, all_labels, colors = [], [], [], [], []
    for r in all_results:
        m = r['m_num']
        all_ra.append(r['ra_deg'])
        all_dec.append(r['dec_deg'])
        all_alt.append(r['max_alt'])
        label = f'M{m}' if r['type'] not in ('BIN',) else ''
        all_labels.append(label)
        # Color: red (high altitude) to gray (low/invisible)
        if r['max_alt'] >= 50:
            colors.append('#e74c3c')
        elif r['max_alt'] >= 30:
            colors.append('#e67e22')
        elif r['max_alt'] >= 15:
            colors.append('#7f8c8d')
        else:
            colors.append('#bdc3c7')

    sc = ax.scatter(all_ra, all_dec, c=colors, s=30, alpha=0.7, edgecolors='white', linewidth=0.3)

    # Labels for top-30 by altitude
    sorted_by_alt = sorted(all_results, key=lambda x: x['max_alt'], reverse=True)
    for r in sorted_by_alt[:30]:
        if r['max_alt'] > 15:
            ax.annotate(f'M{r["m_num"]}', (r['ra_deg'], r['dec_deg']),
                       fontsize=6, ha='left', va='bottom', alpha=0.8,
                       xytext=(3, 3), textcoords='offset points',
                       fontweight='bold')

    # Beijing declination zone: visible area is roughly dec > lat - 90 or dec < 90 - lat... 
    # Objects with dec > lat-90 (for northern hemisphere) are circumpolar
    # Objects with dec < -35 never rise in Beijing (dec < lat - 90 => dec < 39.9 - 90 = -50.1)
    # More practically: objects need dec > -(90 - lat) = -50.1° to rise at all in Beijing
    ax.axhline(y=-50.1, color='gray', linestyle='--', alpha=0.5, label='Beijing horizon limit (Dec > -50.1°)')
    ax.axhline(y=39.9, color='blue', linestyle=':', alpha=0.3, label='Beijing zenith (Dec = 39.9°)')

    # Galactic plane (approximate, sin curve in Mercator)
    gp_ra = np.linspace(0, 360, 360)
    gp_dec = -29.0 + 59.0 * np.sin(np.radians(gp_ra - 283))  # rough approximation
    ax.plot(gp_ra, gp_dec, 'r-', alpha=0.15, linewidth=2, label='Galactic plane (approx.)')

    # Legend
    legend_elements = [
        mpatches.Patch(color='#e74c3c', label='Max alt ≥ 50°'),
        mpatches.Patch(color='#e67e22', label='Max alt 30°–50°'),
        mpatches.Patch(color='#7f8c8d', label='Max alt 15°–30°'),
        mpatches.Patch(color='#bdc3c7', label='Max alt < 15°'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=8)

    ax.set_xlabel('Right Ascension (degrees)')
    ax.set_ylabel('Declination (degrees)')
    ax.set_title('Messier Objects Sky Map — Colored by Maximum Altitude from Beijing (N39.9°)')
    ax.set_xlim(370, -5)  # RA reversed (east to west)
    ax.set_ylim(-90, 90)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  -> {output_file}')


def plot_altitude_heatmap(output_file='output/verify_altitude_heatmap.png'):
    """Heatmap showing altitude of all visible objects over time.
    Y-axis = objects (sorted by best time), X-axis = time, Color = altitude."""
    # Filter visible objects
    visible = [r for r in all_results if r['max_alt'] >= OBJ_ALT_MIN]
    visible.sort(key=lambda x: x['best_time_local'])

    if len(visible) == 0:
        print("  No visible objects for heatmap")
        return

    # Time grid
    n_objs = len(visible)
    times_local = [times_utc[i].datetime + timedelta(hours=UTC_OFFSET) for i in range(len(times_utc))]
    n_times = min(len(times_utc), 288)

    # Build altitude matrix
    alt_matrix = np.zeros((n_objs, n_times))
    for i, r in enumerate(visible):
        alt_curve = r['alt_curve'][:n_times]
        alt_matrix[i, :] = alt_curve

    fig, ax = plt.subplots(figsize=(16, max(8, n_objs * 0.2)))

    im = ax.imshow(alt_matrix, aspect='auto', origin='upper',
                   cmap='inferno', vmin=0, vmax=90,
                   extent=[0, 24, n_objs, 0])

    # Y-axis labels
    y_labels = [f'M{r["m_num"]} {r["name"].split(" ")[1] if " " in r["name"] else ""}'[:18]
                for r in visible]
    ax.set_yticks(np.arange(n_objs) + 0.5)
    ax.set_yticklabels(y_labels, fontsize=6)

    # X-axis: hours
    ax.set_xticks(np.arange(0, 25, 2))
    ax.set_xticklabels([f'{h:02d}:00' for h in range(0, 25, 2)], fontsize=7)
    ax.set_xlabel('Local Time (UTC+8)')

    # Sun altitude contour (approximate)
    from planner import get_sun, BEIJING, SUN_ALT_LIMIT

    # Mark dark window
    sunset_hour = 19.58  # ~19:35
    sunrise_hour = 4.92  # ~04:55
    ax.axvline(x=sunset_hour, color='cyan', linestyle='--', alpha=0.6, linewidth=1)
    ax.axvline(x=sunrise_hour, color='cyan', linestyle='--', alpha=0.6, linewidth=1)
    ax.annotate('Sunset', (sunset_hour, 0), fontsize=6, color='cyan', ha='center', va='bottom')
    ax.annotate('Sunrise', (sunrise_hour, 0), fontsize=6, color='cyan', ha='center', va='bottom')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Altitude (degrees)', fontsize=9)

    ax.set_title(f'Messier Object Altitude Heatmap — Beijing 2026-06-04 ({len(visible)} objects visible)')
    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  -> {output_file}')


def plot_pareto_frontier(output_file='output/verify_pareto.png'):
    """Surface brightness vs max altitude scatter, with Pareto frontier highlighted."""
    visible = [r for r in results if r['max_alt'] >= OBJ_ALT_MIN]

    fig, ax = plt.subplots(figsize=(12, 8))

    # Color by type
    type_colors = {
        'OC': '#e74c3c', 'GC': '#3498db', 'GX': '#2ecc71',
        'EGX': '#27ae60', 'NB': '#f39c12', 'PN': '#9b59b6',
        'SR': '#e67e22', 'SGX': '#1abc9c', 'IGX': '#16a085',
    }
    type_markers = {'OC': 'o', 'GC': 's', 'GX': '^', 'EGX': 'D', 'NB': '*', 'PN': 'P', 'SR': 'X', 'SGX': 'v', 'IGX': '<'}

    for t in type_colors:
        subset = [r for r in visible if r['type'] == t]
        if not subset:
            continue
        x = [r['max_alt'] for r in subset]
        y = [r['sb_arcsec2'] for r in subset]
        sizes = [max(20, min(200, r['composite_score'] * 150)) for r in subset]
        ax.scatter(x, y, c=type_colors[t], marker=type_markers.get(t, 'o'),
                  s=sizes, alpha=0.7, edgecolors='white', linewidth=0.5,
                  label=t, zorder=5)

    # Label top 10
    for r in visible[:10]:
        ax.annotate(f'M{r["m_num"]}', (r['max_alt'], r['sb_arcsec2']),
                   fontsize=7, fontweight='bold', alpha=0.8,
                   xytext=(5, 5), textcoords='offset points')

    # Pareto frontier (lowest SB for each altitude bin)
    alt_bins = np.linspace(15, 90, 30)
    pareto_x, pareto_y = [], []
    for i in range(len(alt_bins) - 1):
        bin_objs = [r for r in visible if alt_bins[i] <= r['max_alt'] < alt_bins[i+1]]
        if bin_objs:
            best = min(bin_objs, key=lambda r: r['sb_arcsec2'])
            pareto_x.append(best['max_alt'])
            pareto_y.append(best['sb_arcsec2'])
    if pareto_x:
        ax.plot(pareto_x, pareto_y, 'k--', alpha=0.4, linewidth=1.5, label='Pareto frontier (best SB at each alt)')

    ax.set_xlabel('Maximum Altitude (degrees)')
    ax.set_ylabel('Surface Brightness (mag/arcsec², lower = brighter)')
    ax.set_title('Messier Object Pareto Analysis — Higher Altitude + Lower SB = Better Target')
    ax.invert_yaxis()
    ax.set_xlim(12, 92)
    ax.legend(fontsize=8, ncol=3, loc='upper left')
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  -> {output_file}')


def plot_fov_comparison(output_file='output/verify_fov_comparison.png'):
    """Show how top targets fit in S50 vs S30 Pro FOV."""
    from planner import S50, S30_PRO

    top_targets = results[:8]  # top 8 for clarity

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    fov_s50_w = 1.28
    fov_s50_h = 0.72
    fov_s30_w = 4.25
    fov_s30_h = 2.39
    pixel_scale = 2.39  # "/px for S50

    for i, r in enumerate(top_targets):
        if i >= 8:
            break
        ax = axes[i]
        size_w = r['size_major'] / 60.0  # arcmin -> deg
        size_h = r['size_minor'] / 60.0

        # S50 FOV box
        s50_rect = plt.Rectangle((-fov_s50_w/2, -fov_s50_h/2), fov_s50_w, fov_s50_h,
                                  fill=False, edgecolor='#3498db', linewidth=2, label='S50 FOV', zorder=5)
        ax.add_patch(s50_rect)

        # S30 Pro FOV box
        s30_rect = plt.Rectangle((-fov_s30_w/2, -fov_s30_h/2), fov_s30_w, fov_s30_h,
                                  fill=False, edgecolor='#e74c3c', linewidth=1.5, linestyle='--', label='S30 Pro FOV', zorder=4)
        ax.add_patch(s30_rect)

        # Target ellipse (approximate)
        target = plt.Circle((0, 0), size_w/2, fill=True, color='#2ecc71', alpha=0.3, zorder=3)
        ax.add_patch(target)
        if size_w != size_h:
            # Draw as ellipse if non-circular
            from matplotlib.patches import Ellipse
            ellipse = Ellipse((0, 0), size_w, size_h, fill=True, color='#2ecc71', alpha=0.5, zorder=3)
            ax.add_patch(ellipse)

        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.5, 2.5)
        ax.set_aspect('equal')
        ax.set_title(f'M{r["m_num"]} {r["name"].split(" ")[1] if " " in r["name"] else r["name"]}\n{r["size_major"]:.0f}\'×{r["size_minor"]:.0f}\'', fontsize=8)
        ax.grid(True, alpha=0.2)
        ax.axhline(y=0, color='gray', alpha=0.2)
        ax.axvline(x=0, color='gray', alpha=0.2)
        if i == 0:
            ax.legend(fontsize=7, loc='upper right')

    plt.suptitle('FOV Comparison: S50 (blue, 1.28°×0.72°) vs S30 Pro (red, 4.25°×2.39°)', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  -> {output_file}')


def plot_timeline_schedule(output_file='output/verify_timeline.png'):
    """Suggested observing timeline for top 15 targets."""
    top15 = results[:15]
    # Sort by best_time_local
    top15 = sorted(top15, key=lambda x: x['best_time_local'])

    fig, ax = plt.subplots(figsize=(14, 8))

    # Convert times to hours
    base_date = datetime(2026, 6, 4)
    times_h = []
    for r in top15:
        bt = datetime.fromisoformat(r['best_time_local'])
        hour = bt.hour + bt.minute / 60.0
        if bt.day > 4:  # next day
            hour += 24
        times_h.append(hour)

    # Plot timeline bars
    y_positions = list(range(len(top15)))
    bar_height = 0.6
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(top15)))

    for i, r in enumerate(top15):
        bt = datetime.fromisoformat(r['best_time_local'])
        hour = bt.hour + bt.minute / 60.0
        if bt.day > 4:
            hour += 24

        # Bar showing exposure window (±1 hour around best time)
        exp_min_s50 = r.get('est_exposure_s50', 10)
        half_width = max(0.25, exp_min_s50 / 120.0)  # at least 30 min window

        start = hour - half_width
        end = hour + half_width

        ax.barh(i, end - start, bar_height, left=start, color=colors[i], alpha=0.7, edgecolor='white')
        ax.plot(hour, i, 'o', color='white', markersize=6, markeredgecolor='black', markeredgewidth=1)

        # Label
        exp_str = f'{exp_min_s50:.0f}min' if exp_min_s50 < 60 else f'{exp_min_s50/60:.1f}h'
        ax.annotate(f'M{r["m_num"]} ({exp_str}, {r["max_alt"]:.0f}°)',
                   (hour, i), fontsize=8, ha='left', va='center',
                   xytext=(5, 0), textcoords='offset points',
                   fontweight='bold')

    # Shade dark period
    sunset_h = 19.58
    sunrise_h = 4.92 + 24  # next day
    ax.axvspan(sunset_h, 24, alpha=0.05, color='blue')
    ax.axvspan(24, sunrise_h, alpha=0.05, color='blue')
    ax.axvspan(0, sunrise_h - 24, alpha=0.05, color='blue')

    ax.set_yticks(y_positions)
    ax.set_yticklabels([f'{i+1}' for i in range(len(top15))])
    ax.set_ylabel('Suggested Order (by best observation time)')

    # X-axis: local time
    ax.set_xlim(19, 32)  # 19:00 to 08:00 next day
    xticks = [20, 22, 0, 2, 4, 6, 8]
    xtick_labels = ['20:00', '22:00', '00:00', '02:00', '04:00', '06:00', '08:00+1']
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.set_xlabel('Local Time (UTC+8)')

    ax.set_title('Suggested Observing Timeline — Top 15 Targets (bar width = S50 exposure window)')
    ax.grid(True, alpha=0.2, axis='x')

    # Annotate moon phase
    ax.annotate('No Moon\n(before 23:00)', (19, -0.5), fontsize=7, color='gray', va='top')

    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  -> {output_file}')


if __name__ == '__main__':
    print("Generating verification visualizations...")
    plot_sky_map()
    plot_altitude_heatmap()
    plot_pareto_frontier()
    plot_fov_comparison()
    plot_timeline_schedule()
    print("Done.")
