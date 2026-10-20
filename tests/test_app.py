"""Streamlit app renders its landing page without exceptions."""

from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parent.parent / "app" / "streamlit_app.py")


def test_app_landing_page_renders():
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert not at.exception
    assert any("Research use only" in w.value for w in at.warning)
