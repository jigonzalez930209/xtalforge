"""Pytest fixtures: a preview server and a Playwright browser."""

from __future__ import annotations

import os
import signal
import socket
import subprocess
import time
import urllib.request

import pytest
from playwright.sync_api import Error as PlaywrightError
from playwright.sync_api import sync_playwright

from support import ROOT

DIST = ROOT / "docs" / ".vitepress" / "dist"


def free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


@pytest.fixture(scope="session")
def base_url() -> str:
    if not (DIST / "playground.html").exists():
        pytest.fail("docs are not built; run scripts/e2e.sh first")
    port = free_port()
    process = subprocess.Popen(
        ["npm", "run", "docs:preview", "--", "--port", str(port)],
        cwd=ROOT / "docs",
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    url = f"http://127.0.0.1:{port}/xtalforge"
    deadline = time.time() + 30
    while time.time() < deadline:
        if process.poll() is not None:
            pytest.fail("the docs preview server exited early")
        try:
            with urllib.request.urlopen(f"{url}/playground", timeout=1) as response:
                if response.status == 200:
                    break
        except OSError:
            time.sleep(0.25)
    else:
        pytest.fail("the docs preview server did not start")
    yield url
    try:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGKILL)


@pytest.fixture(scope="session")
def browser():
    with sync_playwright() as playwright:
        args = ["--use-gl=swiftshader", "--disable-dev-shm-usage"]
        try:
            instance = playwright.chromium.launch(
                channel="chrome", headless=True, args=args
            )
        except PlaywrightError:
            instance = playwright.chromium.launch(headless=True, args=args)
        yield instance
        instance.close()


@pytest.fixture
def page(browser):
    context = browser.new_context(
        viewport={"width": 1440, "height": 1000}, device_scale_factor=1
    )
    current = context.new_page()
    errors: list[str] = []
    current.on("pageerror", lambda error: errors.append(str(error)))
    yield current
    context.close()
    assert not errors, f"page errors: {errors}"
