"""Integration tests for visibility computation (depends on astropy)."""
import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from planner import compute_visibility, rank_objects, S50, OBJ_ALT_MIN


def test_visibility_returns_reasonable_count():
    """Beijing in June should have 50-100 Messier objects visible."""
    results, _ = compute_visibility()
    assert 50 <= len(results) <= 105, f"Expected 50-105 visible, got {len(results)}"


def test_all_results_have_required_fields():
    """Every result should have the standard fields."""
    results, _ = compute_visibility()
    required = {"m_num", "name", "type", "ra_deg", "dec_deg", "vmag",
                "size_major", "size_minor", "sb_arcsec2", "max_alt",
                "best_time_local", "visible_hours"}
    for r in results:
        for field in required:
            assert field in r, f"Missing field '{field}' in M{r.get('m_num')}"


def test_max_alt_above_min():
    """All visible objects should have max_alt >= OBJ_ALT_MIN."""
    results, _ = compute_visibility()
    for r in results:
        assert r["max_alt"] >= OBJ_ALT_MIN, \
            f"M{r['m_num']} max_alt={r['max_alt']} should be >= {OBJ_ALT_MIN}"


def test_surface_brightness_positive():
    """Surface brightness should be positive (typically 15-25 mag/arcsec²)."""
    results, _ = compute_visibility()
    for r in results:
        assert 10 < r["sb_arcsec2"] < 30, \
            f"M{r['m_num']} SB={r['sb_arcsec2']} out of expected range"


def test_rank_adds_scores():
    """Ranking should add composite score and exposure estimates."""
    results, _ = compute_visibility()
    scored = rank_objects(results, S50)
    assert len(scored) > 0
    for r in scored:
        assert "composite_score" in r
        assert "est_exposure_s50" in r
        assert "est_exposure_s30" in r
        assert 0 <= r["composite_score"] <= 1.0
