"""Tests for the Messier catalog data and core calculation functions."""
import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from messier_catalog import MESSIER
from planner import (
    ra_dec_to_deg, compute_surface_brightness,
    fov_calc, pixel_scale, estimate_exposure_minutes,
    S50, S30_PRO, OBJ_ALT_MIN, SUN_ALT_LIMIT
)


def test_catalog_completeness():
    """All 110 Messier objects should be in the catalog."""
    assert len(MESSIER) == 110, f"Expected 110 objects, got {len(MESSIER)}"


def test_catalog_unique():
    """M numbers should be unique."""
    m_numbers = [obj[0] for obj in MESSIER]
    assert len(set(m_numbers)) == 110


def test_catalog_all_types():
    """Every object should have a known type."""
    valid_types = {"OC", "GC", "GX", "EGX", "IGX", "SGX", "NB", "PN", "SR", "AST", "BIN"}
    for obj in MESSIER:
        assert obj[2] in valid_types, f"M{obj[0]} has unknown type: {obj[2]}"


def test_spot_check_m31():
    """M31 Andromeda should be at approximately RA=0h42m, Dec=+41°."""
    for obj in MESSIER:
        if obj[0] == 31:
            ra_deg, dec_deg = ra_dec_to_deg(
                obj[3], obj[4], obj[5], obj[6], obj[7], obj[8], obj[9]
            )
            assert 8 < ra_deg < 13, f"M31 RA should be ~10.7°, got {ra_deg}"
            assert 40 < dec_deg < 42, f"M31 Dec should be ~41°, got {dec_deg}"
            assert obj[10] == 3.4, f"M31 Vmag should be 3.4"
            return
    pytest.fail("M31 not found in catalog")


def test_spot_check_m42():
    """M42 Orion Nebula basic data check."""
    for obj in MESSIER:
        if obj[0] == 42:
            assert obj[10] == 4.0, f"M42 Vmag should be 4.0"
            assert obj[11] == 85.0, f"M42 major axis should be 85'"
            return
    pytest.fail("M42 not found in catalog")


def test_spot_check_m57():
    """M57 Ring Nebula basic data check."""
    for obj in MESSIER:
        if obj[0] == 57:
            assert obj[10] == 8.8, f"M57 Vmag should be 8.8"
            assert 1.0 <= obj[11] <= 1.5, f"M57 size should be ~1.4'"
            assert obj[2] == "PN"
            return
    pytest.fail("M57 not found in catalog")


def test_surface_brightness_known():
    """Surface brightness formula gives expected ranges."""
    # M57: mag=8.8, size=1.4'x1.0' -> SB should be very bright (~17-18)
    sb = compute_surface_brightness(8.8, 1.4, 1.0)
    assert 17 < sb < 19, f"M57 SB should be ~17-18 mag/arcsec², got {sb:.1f}"

    # M33: mag=5.7, size=73'x45' -> low SB (spread out)
    sb = compute_surface_brightness(5.7, 73, 45)
    assert sb > 22, f"M33 SB should be >22 (diffuse galaxy), got {sb:.1f}"


def test_fov_calculation():
    """FOV calculation matches known values."""
    fov_h, fov_v, fov_d = fov_calc(S50["fl_mm"], S50["sensor_w_mm"], S50["sensor_h_mm"])
    assert 1.2 < fov_h < 1.4, f"S50 H FOV should be ~1.28°, got {fov_h:.2f}°"
    assert 0.6 < fov_v < 0.8, f"S50 V FOV should be ~0.72°, got {fov_v:.2f}°"

    fov_h, fov_v, fov_d = fov_calc(S30_PRO["fl_mm"], S30_PRO["sensor_w_mm"], S30_PRO["sensor_h_mm"])
    assert 4.0 < fov_h < 4.5, f"S30 Pro H FOV should be ~4.25°, got {fov_h:.2f}°"
    assert 2.0 < fov_v < 2.6, f"S30 Pro V FOV should be ~2.39°, got {fov_v:.2f}°"


def test_pixel_scale():
    """Pixel scale calculation."""
    ps_s50 = pixel_scale(S50["fl_mm"], S50["pixel_um"])
    assert 2.0 < ps_s50 < 3.0, f"S50 pixel scale should be ~2.39\"/px, got {ps_s50:.2f}"

    ps_s30 = pixel_scale(S30_PRO["fl_mm"], S30_PRO["pixel_um"])
    assert 3.5 < ps_s30 < 4.5, f"S30 Pro pixel scale should be ~3.99\"/px, got {ps_s30:.2f}"


def test_exposure_estimate_ranges():
    """Exposure estimates should be in reasonable ranges for different targets."""
    # Bright PN like M57: SB~17.8, S50 aperture 50mm -> should be short
    t = estimate_exposure_minutes(17.8, 50, "PN")
    assert 0.5 < t < 5, f"Bright PN exposure should be 1-3 min, got {t:.0f}"

    # Dim galaxy like M108: SB~22, S50 -> should be longer
    t = estimate_exposure_minutes(22.0, 50, "GX")
    assert t > 10, f"Dim galaxy exposure should be >10 min, got {t:.0f}"

    # Open cluster: should be very short
    t = estimate_exposure_minutes(19.0, 50, "OC")
    assert t < 10, f"OC exposure should be short, got {t:.0f}"


def test_exposure_aperture_scaling():
    """S30 Pro (30mm) should need ~2.8x more time than S50 (50mm) for same target."""
    t50 = estimate_exposure_minutes(20.0, 50, "GC")
    t30 = estimate_exposure_minutes(20.0, 30, "GC")
    ratio = t30 / t50
    expected_ratio = (50**2) / (30**2)
    assert 2.0 < ratio < 4.0, f"Exposure ratio S30/S50 should be ~{expected_ratio:.1f}, got {ratio:.1f}"


def test_dec_negative_handling():
    """M6 at Dec=-32°13' should give negative declination."""
    for obj in MESSIER:
        if obj[0] == 6:
            ra_deg, dec_deg = ra_dec_to_deg(
                obj[3], obj[4], obj[5], obj[6], obj[7], obj[8], obj[9]
            )
            assert dec_deg < -30, f"M6 Dec should be ~-32°, got {dec_deg}"
            return
    pytest.fail("M6 not found")
