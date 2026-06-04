"""
Messier Marathon planner for Beijing using ZWO SeeStar S30 Pro / S50.
Computes visibility, optimal observation windows, and recommended exposure times.
"""
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, get_sun, Angle
from astropy import units as u
from datetime import datetime, timedelta, UTC
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Use CJK-capable font on macOS
for font_name in ['PingFang SC', 'Heiti SC', 'STHeiti', 'Arial Unicode MS']:
    try:
        matplotlib.font_manager.findfont(font_name, fallback_to_default=False)
        plt.rcParams['font.family'] = font_name
        break
    except Exception:
        continue
from messier_catalog import MESSIER
import json

# ===================== Configuration =====================

BEIJING = EarthLocation(lat=39.9042 * u.deg, lon=116.4074 * u.deg, height=50 * u.m)
DATE = "2026-06-04"  # Beijing local date
UTC_OFFSET = 8  # CST = UTC+8
SUN_ALT_LIMIT = -12  # Nautical twilight, acceptable for smart telescope
OBJ_ALT_MIN = 15  # Minimum altitude for observation (degrees)

# Equipment specs
# S50: 50mm, f/5, FL=250mm, IMX462 (1920x1080, 2.9µm, 5.6x3.2mm)
# S30 Pro: 30mm? (actually ~30mm from name, 150mm FL, IMX585, 4.6° FOV)
S50 = {
    "name": "S50",
    "aperture_mm": 50,
    "f_ratio": 5.0,
    "fl_mm": 250,
    "sensor_w_mm": 5.57,   # IMX462: 1920*2.9µm
    "sensor_h_mm": 3.13,  # IMX462: 1080*2.9µm
    "pixel_um": 2.9,
    "resolution_px": (1920, 1080),
    "qe_pct": 90,
}

S30_PRO = {
    "name": "S30 Pro",
    "aperture_mm": 30,
    "f_ratio": 5.0,
    "fl_mm": 150,
    "sensor_w_mm": 11.14,  # IMX585: 3840*2.9µm
    "sensor_h_mm": 6.26,  # IMX585: 2160*2.9µm
    "pixel_um": 2.9,
    "resolution_px": (3840, 2160),
    "qe_pct": 80,
}


def ra_dec_to_deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, dec_sign):
    """Convert RA/Dec to degrees."""
    ra_deg = (ra_h + ra_m / 60.0 + ra_s / 3600.0) * 15.0
    sign = -1 if dec_d < 0 or dec_sign < 0 else 1
    dec_deg = sign * (abs(dec_d) + dec_m / 60.0 + dec_s / 3600.0)
    return ra_deg, dec_deg


def format_ra_dec(ra_deg, dec_deg):
    """Format decimal RA/Dec as compact sexagesimal strings."""
    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
    ra = coord.ra.to_string(unit=u.hour, sep='hms', precision=0, pad=True)
    dec = coord.dec.to_string(unit=u.deg, sep='°\'"', precision=0, alwayssign=True, pad=True)
    return ra, dec


def compute_surface_brightness(magnitude, size_major, size_minor):
    """Compute average surface brightness in mag/arcsec^2.
    Uses the formula: SB = m + 2.5 * log10(area_in_arcmin^2)
    Then convert: SB(mag/arcsec^2) = SB(mag/arcmin^2) + 8.89
    """
    area_arcmin2 = np.pi * (size_major / 2.0) * (size_minor / 2.0)
    if area_arcmin2 <= 0:
        area_arcmin2 = 0.001
    sb_arcmin2 = magnitude + 2.5 * np.log10(area_arcmin2)
    sb_arcsec2 = sb_arcmin2 + 8.89
    return sb_arcsec2


def fov_calc(fl_mm, sensor_w_mm, sensor_h_mm):
    """Calculate field of view in degrees."""
    fov_h = 2 * np.degrees(np.arctan(sensor_w_mm / (2 * fl_mm)))
    fov_v = 2 * np.degrees(np.arctan(sensor_h_mm / (2 * fl_mm)))
    fov_d = 2 * np.degrees(np.arctan(np.sqrt(sensor_w_mm**2 + sensor_h_mm**2) / (2 * fl_mm)))
    return fov_h, fov_v, fov_d


def pixel_scale(fl_mm, pixel_um):
    """Calculate pixel scale in arcseconds/pixel."""
    return (pixel_um / fl_mm) * 206.265


def estimate_exposure_minutes(sb_arcsec2, aperture_mm, target_type):
    """
    Estimate recommended total integration time in minutes.
    Based on surface brightness and aperture, with type-specific adjustments.
    Brighter SB (< 20 mag/arcsec2) = less time; dimmer = more time.
    """
    # Base: reference SB of 20 mag/arcsec2, 50mm aperture needs ~10 min
    effective_aperture_factor = (2500.0 / (aperture_mm**2))  # normalized to 50mm
    sb_diff = sb_arcsec2 - 20.0
    base_minutes = effective_aperture_factor * 10 * (10 ** (0.4 * sb_diff))

    # Type adjustments
    type_multiplier = {
        "OC": 0.3,   # Open clusters are easy
        "GC": 0.6,   # Globulars need more time to resolve core
        "NB": 0.8,   # Nebulae are usually OK
        "PN": 0.7,   # Planetary nebulae are compact
        "GX": 1.5,   # Galaxies need the most time
        "EGX": 1.5,
        "IGX": 1.5,
        "SGX": 1.5,
        "SR": 1.0,
        "AST": 0.2,
        "BIN": 0.1,
    }
    multiplier = type_multiplier.get(target_type, 1.0)
    minutes = base_minutes * multiplier
    # Clamp to reasonable range
    return max(1, min(minutes, 240))


def compute_visibility():
    """Compute visibility of all Messier objects from Beijing on the given night."""
    # Time range: from 2 hours before sunset to 2 hours after sunrise
    t0 = Time(f"{DATE}T12:00:00") - 12 * u.hour  # 2026-06-04 08:00 Beijing time
    # Create time grid every 5 minutes
    n_steps = 24 * 12  # 24 hours at 5-min intervals
    times_utc = t0 + np.arange(n_steps) * 5 * u.minute

    # Sun position
    sun_altaz = get_sun(times_utc).transform_to(AltAz(obstime=times_utc, location=BEIJING))

    results = []
    for obj in MESSIER:
        m_num, name, otype = obj[0], obj[1], obj[2]
        ra_h, ra_m, ra_s = obj[3], obj[4], obj[5]
        dec_d, dec_m_arg, dec_s, dec_sign = obj[6], obj[7], obj[8], obj[9]
        vmag, size_major, size_minor = obj[10], obj[11], obj[12]

        ra_deg, dec_deg = ra_dec_to_deg(ra_h, ra_m, ra_s, dec_d, dec_m_arg, dec_s, dec_sign)
        coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
        altaz = coords.transform_to(AltAz(obstime=times_utc, location=BEIJING))
        altitudes = altaz.alt.deg
        azimuths = altaz.az.deg

        # Filter: object alt > limit AND sun alt < limit
        visible_mask = (altitudes > OBJ_ALT_MIN) & (sun_altaz.alt.deg < SUN_ALT_LIMIT)
        visible_indices = np.where(visible_mask)[0]

        if len(visible_indices) == 0:
            continue

        max_alt = np.max(altitudes[visible_mask])
        best_idx = visible_indices[np.argmax(altitudes[visible_mask])]
        best_time_utc = times_utc[best_idx]
        best_time_local = best_time_utc + UTC_OFFSET * u.hour
        best_az = azimuths[best_idx]

        # Duration above 15° with dark sky
        dark_mask = (altitudes > OBJ_ALT_MIN) & (sun_altaz.alt.deg < SUN_ALT_LIMIT)
        dark_indices = np.where(dark_mask)[0]
        visible_hours = len(dark_indices) * 5.0 / 60.0

        # Compute surface brightness
        sb = compute_surface_brightness(vmag, size_major, size_minor)

        # Track altitude curve
        alt_curve = altitudes.copy()
        sun_curve = sun_altaz.alt.deg.copy()

        results.append({
            "m_num": m_num,
            "name": name,
            "type": otype,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "vmag": vmag,
            "size_major": size_major,
            "size_minor": size_minor,
            "sb_arcsec2": sb,
            "max_alt": max_alt,
            "best_time_local": best_time_local.datetime,
            "visible_hours": visible_hours,
            "best_az": best_az,
            "alt_curve": alt_curve,
            "sun_curve": sun_curve,
        })

    # Sort by max altitude descending
    results.sort(key=lambda x: x["max_alt"], reverse=True)
    return results, times_utc


def rank_objects(results, equipment):
    """
    Rank objects by a composite score:
    - altitude bonus (higher = better)
    - brightness bonus (lower sb = better)
    - visibility duration bonus
    - size-to-FOV fit bonus
    """
    fov_h, fov_v, fov_d = fov_calc(equipment["fl_mm"], equipment["sensor_w_mm"], equipment["sensor_h_mm"])

    for r in results:
        sb = r["sb_arcsec2"]
        max_alt = r["max_alt"]
        hours = r["visible_hours"]
        size_major = r["size_major"] / 60.0  # convert arcmin to degrees
        size_minor = r["size_minor"] / 60.0

        # FOV fit: how well the object fits in frame (1.0 = perfect fit)
        fov_fit_w = size_major / (fov_h * 0.8)
        fov_fit_h = size_minor / (fov_v * 0.8)
        fov_fit = 1.0 / (1.0 + max(0, max(fov_fit_w, fov_fit_h) - 1.0) * 3)

        # Altitude score: 0 at 15°, 1 at 90°
        alt_score = (max_alt - OBJ_ALT_MIN) / (90.0 - OBJ_ALT_MIN)
        alt_score = max(0, min(1, alt_score))

        # Brightness score: lower SB = better. 18 = excellent, 24 = poor
        # Normalize: 1.0 at SB=18, 0.0 at SB=24
        brightness_score = max(0, min(1, (24 - sb) / 6.0))

        # Duration score (hours visible)
        dur_score = min(1.0, hours / 8.0)

        # Composite
        composite = 0.30 * alt_score + 0.35 * brightness_score + 0.10 * dur_score + 0.25 * fov_fit
        r["composite_score"] = composite
        r["fov_fit"] = fov_fit
        r["alt_score"] = alt_score
        r["brightness_score"] = brightness_score
        r["dur_score"] = dur_score
        r["est_exposure_s50"] = estimate_exposure_minutes(sb, S50["aperture_mm"], r["type"])
        r["est_exposure_s30"] = estimate_exposure_minutes(sb, S30_PRO["aperture_mm"], r["type"])

    results.sort(key=lambda x: x["composite_score"], reverse=True)
    return results


def plot_altitude_curves(results, times_utc, top_n=20, output_file="altitude_curves.png"):
    """Plot altitude vs time for the top N objects."""
    fig, ax = plt.subplots(figsize=(14, 8))

    times_local = times_utc.datetime + timedelta(hours=UTC_OFFSET)
    times_num = mdates.date2num(times_local)

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, min(top_n, len(results))))
    for i, r in enumerate(results[:top_n]):
        alt = r["alt_curve"]
        mask = alt > 0  # Only plot above horizon
        if np.sum(mask) > 0:
            ax.plot(times_num[mask], alt[mask], color=colors[i],
                    label=f'{r["m_num"]} {r["name"]} ({r["max_alt"]:.0f}°)',
                    linewidth=1.5, alpha=0.85)

    # Sun altitude
    from astropy.coordinates import get_sun
    sun_altaz = get_sun(times_utc).transform_to(AltAz(obstime=times_utc, location=BEIJING))
    sun_alt = sun_altaz.alt.deg
    ax.plot(times_num, sun_alt, 'r--', linewidth=2, alpha=0.5, label='Sun')
    ax.axhline(y=SUN_ALT_LIMIT, color='gray', linestyle=':', alpha=0.5, label=f'Sun {SUN_ALT_LIMIT}°')

    ax.axhline(y=OBJ_ALT_MIN, color='gray', linestyle=':', alpha=0.3)

    ax.set_xlabel('Local Time (UTC+8)')
    ax.set_ylabel('Altitude (degrees)')
    ax.set_title(f'Messier Object Altitudes from Beijing - Night of {DATE}')
    ax.set_ylim(0, 90)

    # X-axis formatting
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    ax.legend(loc='upper left', fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    return output_file


def plot_top10_sun_altitude(results, times_utc, output_file="top10_sun_altitude.png"):
    """Plot Sun plus the 10 recommended targets on one altitude-vs-time chart."""
    fig, ax = plt.subplots(figsize=(14, 8))

    times_local = times_utc.datetime + timedelta(hours=UTC_OFFSET)
    times_num = mdates.date2num(times_local)

    sun_altaz = get_sun(times_utc).transform_to(AltAz(obstime=times_utc, location=BEIJING))
    sun_alt = sun_altaz.alt.deg
    ax.plot(times_num, sun_alt, color="#e74c3c", linewidth=2.5, linestyle="--", label="Sun")
    ax.axhline(y=SUN_ALT_LIMIT, color="#7f8c8d", linestyle=":", linewidth=1.5, label=f"Sun {SUN_ALT_LIMIT}°")
    ax.axhline(y=OBJ_ALT_MIN, color="#95a5a6", linestyle=":", linewidth=1.0, label=f"Target {OBJ_ALT_MIN}°")
    ax.fill_between(times_num, 0, 90, where=sun_alt < SUN_ALT_LIMIT, color="#2c3e50", alpha=0.08, label="usable dark window")

    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    for i, r in enumerate(results[:10]):
        alt = r["alt_curve"]
        mask = alt > 0
        label = f"{i+1}. M{r['m_num']} ({r['max_alt']:.0f}°, {r['best_time_local'].strftime('%H:%M')})"
        ax.plot(times_num[mask], alt[mask], color=colors[i], linewidth=1.8, alpha=0.9, label=label)

        best_time_num = mdates.date2num(r["best_time_local"])
        ax.scatter([best_time_num], [r["max_alt"]], color=colors[i], s=30, edgecolor="white", linewidth=0.8, zorder=5)

    ax.set_xlabel("Local Time (UTC+8)")
    ax.set_ylabel("Altitude (degrees)")
    ax.set_title(f"Sun and Top 10 Recommended Messier Targets — Beijing {DATE}")
    ax.set_ylim(-25, 90)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax.legend(loc="upper left", fontsize=7, ncol=2)
    ax.grid(True, alpha=0.25)

    plt.tight_layout()
    plt.savefig(output_file, dpi=180, bbox_inches="tight")
    plt.close()
    return output_file


def plot_scoring_bars(results, top_n=20, output_file="scoring_bars.png"):
    """Plot scoring breakdown for top N objects."""
    fig, ax = plt.subplots(figsize=(12, max(6, top_n * 0.35)))

    items = results[:top_n][::-1]  # Reverse for horizontal bar (bottom = best)
    labels = [f'M{r["m_num"]} {r["name"].split(" ")[1] if " " in r["name"] else r["name"]}'
              for r in items]

    alt_s = [r["alt_score"] for r in items]
    brt_s = [r["brightness_score"] for r in items]
    dur_s = [r["dur_score"] for r in items]
    fov_s = [r["fov_fit"] for r in items]

    x = np.arange(len(labels))
    width = 0.2
    ax.barh(x - 1.5*width, alt_s, width, label='Altitude', color='#3498db')
    ax.barh(x - 0.5*width, brt_s, width, label='Brightness', color='#e74c3c')
    ax.barh(x + 0.5*width, dur_s, width, label='Duration', color='#2ecc71')
    ax.barh(x + 1.5*width, fov_s, width, label='FOV Fit', color='#f39c12')

    ax.set_yticks(x)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('Score')
    ax.set_title(f'Messier Object Scoring Breakdown - Beijing {DATE}')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.2, axis='x')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    return output_file


def plot_sky_map(results, output_file="sky_map.png"):
    """Simple sky map of visible objects at their best observation time (alt/az)."""
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})

    for r in results:
        theta = np.radians(90 - r.get("best_az", 0))
        r_val = 90 - r["max_alt"]

        size = max(10, min(200, r["composite_score"] * 150))
        color = plt.cm.plasma(r["composite_score"])

        ax.scatter(theta, r_val, s=size, c=[color], alpha=0.8, edgecolors='white', linewidth=0.5)
        ax.annotate(f'M{r["m_num"]}', (theta, r_val), fontsize=5, ha='center', va='center',
                    color='white', fontweight='bold')

    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(90, 0)
    ax.set_yticks([0, 30, 60, 90])
    ax.set_yticklabels(['90°', '60°', '30°', '0°'])
    ax.set_title(f'Messier Visibility Map - Beijing {DATE}\n(radius = zenith distance, larger dot = higher score)')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    return output_file


def generate_report(results, times_utc):
    """Generate Markdown report."""
    fov_s50_h, fov_s50_v, fov_s50_d = fov_calc(S50["fl_mm"], S50["sensor_w_mm"], S50["sensor_h_mm"])
    fov_s30_h, fov_s30_v, fov_s30_d = fov_calc(S30_PRO["fl_mm"], S30_PRO["sensor_w_mm"], S30_PRO["sensor_h_mm"])
    ps_s50 = pixel_scale(S50["fl_mm"], S50["pixel_um"])
    ps_s30 = pixel_scale(S30_PRO["fl_mm"], S30_PRO["pixel_um"])

    # Sun times: compute then find sunset and sunrise
    from astropy.coordinates import get_sun
    t0 = Time(f"{DATE}T12:00:00") - 12 * u.hour
    n_steps = 24 * 12
    times_utc_full = t0 + np.arange(n_steps) * 5 * u.minute
    sun_alt = get_sun(times_utc_full).transform_to(AltAz(obstime=times_utc_full, location=BEIJING)).alt.deg
    total_steps = len(sun_alt)
    # Sunset: first index where sun alt < 0 after being positive
    sunset_idx = None
    for i in range(1, total_steps):
        if sun_alt[i-1] > 0 and sun_alt[i] <= 0:
            sunset_idx = i
            break
    if sunset_idx is None:
        sunset_idx = 0

    # Sunrise: first index where sun alt > 0 after being negative
    sunrise_idx = None
    for i in range(sunset_idx + 1, total_steps):
        if sun_alt[i-1] < 0 and sun_alt[i] >= 0:
            sunrise_idx = i
            break
    if sunrise_idx is None:
        sunrise_idx = total_steps - 1

    sunset_time = (times_utc_full[sunset_idx] + UTC_OFFSET * u.hour).datetime
    sunrise_time = (times_utc_full[sunrise_idx] + UTC_OFFSET * u.hour).datetime

    # Find dark window (sun alt < SUN_ALT_LIMIT)
    for i in range(total_steps):
        if sun_alt[i] < SUN_ALT_LIMIT:
            dark_start_idx = i
            break
    else:
        dark_start_idx = sunset_idx

    for i in range(dark_start_idx, total_steps):
        if sun_alt[i] >= SUN_ALT_LIMIT:
            dark_end_idx = i
            break
    else:
        dark_end_idx = total_steps - 1

    dark_start_time = (times_utc_full[dark_start_idx] + UTC_OFFSET * u.hour).datetime
    dark_end_time = (times_utc_full[dark_end_idx] + UTC_OFFSET * u.hour).datetime

    report = f"""# 北京今晚梅西耶马拉松观测规划

**日期**：2026年6月4日  
**地点**：北京（北纬 39.9°，东经 116.4°，海拔 50m）  
**设备**：ZWO SeeStar S30 Pro / S50 智能天文望远镜

---

## 设备参数

| 参数 | S50 | S30 Pro |
|------|-----|---------|
| 口径 | 50mm | 30mm |
| 焦比 | f/5 | f/5 |
| 焦距 | 250mm | 150mm |
| 传感器 | IMX462 | IMX585 |
| 分辨率 | 1920×1080 | 3840×2160 |
| 像元 | 2.9µm | 2.9µm |
| 视场 (H×V) | {fov_s50_h:.2f}°×{fov_s50_v:.2f}° | {fov_s30_h:.2f}°×{fov_s30_v:.2f}° |
| 视场 (对角线) | {fov_s50_d:.2f}° | {fov_s30_d:.2f}° |
| 像素分辨率 | {ps_s50:.2f}"/px | {ps_s30:.2f}"/px |
| 集光面积 | 1963mm² | 707mm² |

## 夜空条件

- 日落：{sunset_time.strftime('%H:%M')} CST
- 日出：{sunrise_time.strftime('%H:%M')} CST（次日）
- 可用暗夜窗口：{dark_start_time.strftime('%H:%M')} – {dark_end_time.strftime('%H:%M')} CST（太阳高度 < {SUN_ALT_LIMIT}°，航海晨昏）
- 月亮：6月4日为农历十九，下弦月附近，前半夜无月光干扰

## 六月梅西耶马拉松特点

传统梅西耶马拉松（3月下旬至4月初）可以在一夜之内拍到全部110个梅西耶天体。但六月的北京处于夏至前夕，黑夜时间短（约{((dark_end_time - dark_start_time).total_seconds() / 3600):.1f}小时），加上白昼夜晚的交替，很多冬季和春季的目标（如猎户座大星云、昴星团等）在太阳附近无法拍摄。

六月的优势集中在银河中心区域——人马座、天蝎座、蛇夫座一带，这是银河系最密集的天区，包含了大量明亮的星团、星云和球状星团。S30 Pro 的宽视场（{fov_s30_h:.1f}°×{fov_s30_v:.1f}°）在这个季节特别有优势，能一次性框住多个人马座星云。

## 观测规划方法

对每个梅西耶天体，我们计算了以下指标：

1. **可见窗口**：天体高度 > {OBJ_ALT_MIN}° 且太阳高度 < {SUN_ALT_LIMIT}° 的时间段
2. **暗夜窗口最高点时间**：在太阳高度 < {SUN_ALT_LIMIT}° 且目标高度 > {OBJ_ALT_MIN}° 的可拍摄窗口内，目标达到最大高度的时间
3. **面亮度**：根据总视星等和角面积计算的平均面亮度（mag/arcsec²）
4. **FOV匹配度**：天体角大小与设备视场的匹配程度
5. **综合评分**：加权组合以上因素（高度30% + 亮度35% + 可见时长10% + FOV匹配25%）
6. **推荐曝光**：根据面亮度、设备口径和天体类型估算的最佳累计曝光时间

## 主验证图：太阳与 Top 10 目标高度

这张图用于最快速地查证推荐是否合理。横轴是北京时间，纵轴是高度角。红色虚线是太阳，灰色虚线是太阳 -12° 暗夜线，淡蓝色阴影是可用拍摄窗口。每条彩色曲线是一个推荐目标，曲线上的圆点表示它在暗夜窗口内的最高位置。

读图方法很直接：一个目标越早进入阴影区、曲线越高、停留时间越长，就越适合今晚拍摄。Top 10 中 M57、M29、M92、M94、M13 都在暗夜窗口内达到 80° 以上，这也是它们进入推荐列表前列的主要原因。

![](top10_sun_altitude.png)

"""

    # Summary table
    report += "## 最佳观测目标排名（Top 25）\n\n"
    report += "| 排名 | 编号 | 名称 | 类型 | 星等 | 大小 | 面亮度 | 窗口最高 | 高度 | 可见 | S50曝光 | S30曝光 | 评分 |\n"
    report += "|------|------|------|------|------|------|--------|--------|------|------|----------|----------|------|\n"

    for i, r in enumerate(results[:25]):
        bt = r["best_time_local"]
        sb = r["sb_arcsec2"]
        size_str = f'{r["size_major"]:.0f}\'×{r["size_minor"]:.0f}\'' if r["size_major"] != r["size_minor"] else f'{r["size_major"]:.0f}\''
        report += f"| {i+1} | M{r['m_num']} | {r['name']} | {r['type']} | {r['vmag']:.1f} | {size_str} | {sb:.1f} | {bt.strftime('%H:%M')} | {r['max_alt']:.0f}° | {r['visible_hours']:.1f}h | {r['est_exposure_s50']:.0f}min | {r['est_exposure_s30']:.0f}min | {r['composite_score']:.2f} |\n"

    report += f"""
## 完整可见深空目标列表（共 {len(results)} 个）

{len(results)} 个梅西耶深空目标在今晚有观测窗口。M40 双星、M73 星群和 M24 恒星云已从推荐排序中过滤，因为它们不适合作为本设备的常规深空拍摄目标。按综合评分完整排序如下。

"""

    # All objects table
    report += "| 排名 | 编号 | 名称 | 类型 | 星等 | 面亮度 | 窗口最高 | 高度 | 可见 | 评分 |\n"
    report += "|------|------|------|------|------|--------|--------|------|------|------|\n"
    for i, r in enumerate(results):
        bt = r["best_time_local"]
        sb = r["sb_arcsec2"]
        report += f"| {i+1} | M{r['m_num']} | {r['name']} | {r['type']} | {r['vmag']:.1f} | {sb:.1f} | {bt.strftime('%H:%M')} | {r['max_alt']:.0f}° | {r['visible_hours']:.1f}h | {r['composite_score']:.2f} |\n"

    report += f"""
## 按类型分布

"""
    type_counts = {}
    for r in results:
        t = r["type"]
        type_counts[t] = type_counts.get(t, 0) + 1
    type_names = {"OC": "疏散星团", "GC": "球状星团", "GX": "旋涡星系", "EGX": "椭圆星系",
                  "NB": "弥漫星云", "PN": "行星状星云", "SR": "超新星遗迹", "IGX": "不规则星系",
                  "SGX": "透镜星系", "AST": "星云/恒星云", "BIN": "双星"}
    for t, c in sorted(type_counts.items(), key=lambda x: x[1], reverse=True):
        report += f"- **{type_names.get(t, t)}**（{t}）：{c} 个\n"

    # Top 5 recommendations with physical facts
    report += "\n## 今晚最佳 5 个推荐目标\n\n"

    top5 = results[:5]
    for i, r in enumerate(top5):
        report += f"### {i+1}. M{r['m_num']} {r['name']} — 综合评分 {r['composite_score']:.2f}\n\n"
        report += f"- **类型**：{type_names.get(r['type'], r['type'])}\n"
        report += f"- **最佳观测时间**：{r['best_time_local'].strftime('%H:%M')} CST，暗夜窗口高度 {r['max_alt']:.0f}°\n"
        report += f"- **视星等**：{r['vmag']:.1f}，面亮度：{r['sb_arcsec2']:.1f} mag/arcsec²\n"
        report += f"- **角大小**：{r['size_major']:.0f}'×{r['size_minor']:.0f}'\n"
        report += f"- **建议曝光**：S50 约 {r['est_exposure_s50']:.0f} 分钟 / S30 Pro 约 {r['est_exposure_s30']:.0f} 分钟\n"
        report += f"- **选择理由**：暗夜窗口高度 {r['max_alt']:.0f}°（高空有利减少大气消光），面亮度 {r['sb_arcsec2']:.1f} mag/arcsec²，可见窗口 {r['visible_hours']:.1f} 小时\n\n"

    report += """## 复用这个规划工具

这个报告背后的完整代码已经整理成公开仓库：[github.com/grapeot/messier-marathon-planner](https://github.com/grapeot/messier-marathon-planner)。仓库包含梅西耶星表、可见性计算、评分逻辑、测试和可视化脚本，`output/report.md` 是由代码直接生成的结果。

仓库里也包含一个可交给 AI agent 使用的 workflow skill：[`skills/messier_marathon_planner.md`](https://github.com/grapeot/messier-marathon-planner/blob/master/skills/messier_marathon_planner.md)。想换城市、日期、设备或筛选条件时，可以把仓库链接交给支持项目说明/skills 的 AI coding agent，并说明你的观测地点和设备。Agent 应先安装这个 skill，再按其中的流程修改配置、运行测试、重新生成报告和验证图。

一个可直接使用的提示词是：

> 请安装并使用这个 skill 帮我规划梅西耶马拉松：https://github.com/grapeot/messier-marathon-planner 。观测地点是[城市/经纬度]，日期是[当地日期]，设备是[望远镜/相机参数]。

"""

    # Cross-verification for top 5
    report += "---\n\n## 附录 A：5 个推荐目标的交叉验证\n\n"
    report += "以下是每个推荐目标的可公开验证链接。读者可以使用这些链接独立核对我们计算的高度、面亮度和可见性。\n\n"
    top5_verify = top_5_verification(top5)
    for item in top5_verify:
        report += item

    # Appendix B: Calculation details
    report += f"""---
## 附录 B：计算过程说明

### B.1 坐标与可见性

使用 Astropy 库计算天体在观察者地平坐标系中的高度角（Altitude）和方位角（Azimuth）。对每个梅西耶天体，以 5 分钟为步长计算当晚全程的高度变化曲线。

```python
# 核心计算（简化版）
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time

beijing = EarthLocation(lat=39.9, lon=116.4, height=50)
times = Time("2026-06-04T00:00:00") + np.arange(288) * 5*u.min  # UTC, 24h, 5min steps

for obj in messier_catalog:
    coord = SkyCoord(ra=ra_deg, dec=dec_deg, frame='icrs')
    altaz = coord.transform_to(AltAz(obstime=times, location=beijing))
    altitudes = altaz.alt.deg  # 天体的高度角（度）
```

可见条件：
- 天体高度 > 15°（低于此角度大气消光严重，跟踪也可能不稳）
- 太阳高度 < -12°（航海晨昏，对智能望远镜的 stacked imaging 足够）

### B.2 面亮度计算

平均面亮度（mag/arcmin²）= 视星等 + 2.5 × log₁₀(面积[arcmin²])

其中面积 = π × (长轴/2) × (短轴/2)，圆面天体取长短轴相等。

转换为 mag/arcsec²：加 8.89（因为 1 arcmin² = 3600 arcsec²，log₁₀(3600) ≈ 3.556，×2.5 = 8.89）。

这一面亮度是**平均**面亮度。实际观测中，星系中心比平均亮很多，星云的 HII 区域亮度也不均匀。因此这个数值适合做相对比较，但曝光时间要根据目标本身的特点调整。比如 M31 虽然平均面亮度不低，但核球附近非常亮，实际更容易拍到。

### B.3 曝光时间估算

基础公式参考了 Clark 的曝光模型，简化后的经验公式：

```
曝光时间（分钟）= 基准时间 × (2500 / 口径²) × 10^(0.4 × (面亮度 - 20))
```

基准时间为 50mm 口径、面亮度 20 mag/arcsec² 时约需 10 分钟。不同类型目标乘以经验系数：疏散星团 0.3×（亮星多），球状星团 0.6×，星系 1.5×（需分辨旋臂细节）。

对于 Seestar 这类智能望远镜，单张曝光 10 秒，通过 stacking 实现总曝光。推荐的分钟值是指总累计时间，由设备自动堆叠完成。

### B.4 综合评分公式

```
综合评分 = 0.30 × 高度得分 + 0.35 × 亮度得分 + 0.10 × 可见时长得分 + 0.25 × FOV匹配得分
```

各子得分均归一化到 [0, 1] 区间：
- 高度得分 = (高度 - 15) / 75（15° 为 0，90° 为 1）
- 亮度得分 = (24 - 面亮度) / 6（面亮度 18 对应满分，24 对应 0 分）
- 时长得分 = 可见小时 / 8
- FOV 匹配得分：天体角大小在设备视场 80% 以内为满分，超出按比例扣分

权重设计：亮度和 FOV 匹配占比较高（60%），因为对拍摄质量影响最大；高度占 30% 反映大气消光的重要性；可见时长权重较低是因为 Seestar 自动跟踪，只要窗口够拍即可。

### B.5 设备 FOV 计算

视场角 = 2 × arctan(传感器尺寸 / (2 × 焦距))

- S50：传感器 5.6×3.2mm，焦距 250mm → FOV 约 1.28°×0.73°
- S30 Pro：传感器约 11.1×6.3mm，焦距 150mm → FOV 约 4.25°×2.39°

像素分辨率 = 像元尺寸(µm) / 焦距(mm) × 206.265

- S50：2.9 / 250 × 206.265 ≈ 2.39"/px
- S30 Pro：2.9 / 150 × 206.265 ≈ 3.99"/px

---

*报告生成时间：{datetime.now(UTC).strftime('%Y-%m-%d %H:%M:%S')} UTC*  
*数据来源：[SEDS Messier Catalog](http://spider.seds.org/), ZWO 官方参数*  
*源代码：[github.com/grapeot/messier-marathon-planner](https://github.com/grapeot/messier-marathon-planner)*  
*本报告为临时版本，仅供内部参考*
"""
    return report


def top_5_verification(top5):
    """Generate cross-verification appendix: public URLs to independently check each target's visibility."""
    result = []
    for r in top5:
        m_num = r["m_num"]
        ra_str, dec_str = format_ra_dec(r["ra_deg"], r["dec_deg"])
        size_str = f'{r["size_major"]:.0f}\'×{r["size_minor"]:.0f}\'' if r["size_major"] != r["size_minor"] else f'{r["size_major"]:.0f}\''
        best_time = r["best_time_local"].strftime('%Y-%m-%d %H:%M CST')
        result.append(f"""
### M{m_num} {r['name']}

**坐标**：RA {ra_str}, Dec {dec_str}

**公开验证链接**：

| 工具 | 链接 | 验证要点 |
|------|------|----------|
| SEDS Messier Catalog | [https://messier.seds.org/m/m{m_num:03d}.html](https://messier.seds.org/m/m{m_num:03d}.html) | 查看视星等、坐标、角大小 |
| In-The-Sky (天体搜索) | [https://in-the-sky.org/search.php?term=M{m_num}](https://in-the-sky.org/search.php?term=M{m_num}) | 设置北京观测位置后可查看今晚的升起/中天/下落时间 |
| Stellarium Web | [https://stellarium-web.org](https://stellarium-web.org) | 设置 Beijing、时间 {best_time}，搜索 M{m_num} 验证高度 |

**我们的计算结论**：

{best_time}，M{m_num} 在暗夜窗口内高度约 {r['max_alt']:.0f}°。视星等 {r['vmag']:.1f}，角大小 {size_str}，平均面亮度约 {r['sb_arcsec2']:.1f} mag/arcsec²，可见窗口约 {r['visible_hours']:.1f} 小时。

**快速几何校验**：北京纬度 39.9°，该目标赤纬 {r['dec_deg']:+.1f}°，理论中天最高高度约 {90 - abs(39.9042 - r['dec_deg']):.0f}°。报告中的高度取自暗夜可拍摄窗口，因此可能低于全天中天高度。如果外部工具显示的高度与这一数量级不一致，应优先检查位置、日期和时区设置。

---

**补充说明**：

- 也可以使用 [Telescopius DSO Search](https://telescopius.com/deep-sky/search) 设置位置 Beijing、时间 2026-06-04 晚间，过滤高度 > 30° 的目标，查看我们的排名是否与 Telescopius 的推荐一致
- [Fourmilab Your Sky](https://www.fourmilab.ch/yoursky/) 可直接生成北京当晚的星图
- 所有高度计算基于 Astropy 的 `AltAz` 坐标系变换，太阳位置使用 `get_sun()` 函数。如果你得到的结果与我们有差异，请检查观测位置的经纬度设置是否一致（北京: 39.9042°N, 116.4074°E）
""")
    return result


def main():
    print("计算北京今晚梅西耶天体可见性...")
    results, times_utc = compute_visibility()
    print(f"  有观测窗口的目标：{len(results)} 个")

    # Rank for S50
    print("评估 S50 观测方案...")
    # Filter out non-deep-sky objects (binary stars, asterisms) as uninteresting
    results_dso = [r for r in results if r["type"] not in ("BIN", "AST")]
    ranked_s50 = rank_objects(results_dso, S50)

    # Rank for S30 Pro
    print("评估 S30 Pro 观测方案...")
    ranked_s30 = rank_objects(list(results), S30_PRO)

    # Generate plots (use S50 ranking as primary since S50 is the higher-spec device)
    print("生成可视化图表...")
    plot_altitude_curves(ranked_s50, times_utc, top_n=20, output_file="output/altitude_curves.png")
    plot_top10_sun_altitude(ranked_s50, times_utc, output_file="output/top10_sun_altitude.png")
    plot_scoring_bars(ranked_s50, top_n=20, output_file="output/scoring_bars.png")

    # Generate report
    print("生成报告...")
    report = generate_report(ranked_s50, times_utc)

    # Save report
    with open("output/report.md", "w") as f:
        f.write(report)

    # Save JSON results
    json_data = []
    for r in ranked_s50:
        entry = {k: v for k, v in r.items() if k not in ('alt_curve', 'sun_curve')}
        entry["best_time_local"] = entry["best_time_local"].strftime('%Y-%m-%dT%H:%M:%S')
        json_data.append(entry)
    with open("output/results.json", "w") as f:
        json.dump(json_data, f, indent=2, ensure_ascii=False)

    print(f"""
=== 完成 ===
- 报告: output/report.md
- 高度曲线图: output/altitude_curves.png
- 太阳+Top10高度图: output/top10_sun_altitude.png
- 评分图: output/scoring_bars.png
- JSON数据: output/results.json
""")


if __name__ == "__main__":
    main()
