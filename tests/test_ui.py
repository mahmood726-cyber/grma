import contextlib
import http.server
import socket
import socketserver
import subprocess
import sys
import threading
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait


ROOT = Path(__file__).resolve().parents[1]
HOST = "127.0.0.1"
PORT = 8000


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, directory=None, **kwargs):
        super().__init__(*args, directory=str(ROOT), **kwargs)

    def log_message(self, format, *args):  # noqa: A003
        pass


class ReusableTCPServer(socketserver.TCPServer):
    allow_reuse_address = True


def _find_port():
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        sock.bind((HOST, 0))
        return sock.getsockname()[1]


def _base_url():
    port = PORT
    try:
        server = ReusableTCPServer((HOST, port), QuietHandler)
    except OSError:
        port = _find_port()
        server = ReusableTCPServer((HOST, port), QuietHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server, thread, f"http://{HOST}:{port}/index.html"


def _driver():
    options = Options()
    options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--disable-gpu")
    options.add_argument("--window-size=1440,1000")
    return webdriver.Chrome(options=options)


def test_landing_page_loads():
    server, thread, url = _base_url()
    browser = _driver()
    try:
        browser.get(url)
        WebDriverWait(browser, 20).until(lambda d: d.find_element(By.TAG_NAME, "h1"))
        assert "Grey Relational Meta-Analysis" in browser.title
        assert "Grey Relational Meta-Analysis" in browser.find_element(By.TAG_NAME, "h1").text
    finally:
        browser.quit()
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_landing_page_exposes_expected_links_and_figures():
    server, thread, url = _base_url()
    browser = _driver()
    try:
        browser.get(url)
        WebDriverWait(browser, 20).until(lambda d: len(d.find_elements(By.CSS_SELECTOR, ".links .card")) >= 3)
        hrefs = browser.execute_script(
            "return Array.from(document.querySelectorAll('.links .card')).map((el) => el.getAttribute('href'));"
        )
        assert "e156-submission/index.html" in hrefs
        assert "e156-submission/paper.md" in hrefs
        assert "e156-submission/protocol.md" in hrefs
        assert len(browser.find_elements(By.CSS_SELECTOR, ".figs img")) >= 5
    finally:
        browser.quit()
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_manuscript_verification_script_runs():
    result = subprocess.run(
        [sys.executable, str(ROOT / "verify_manuscript_numbers.py")],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stderr or result.stdout
    assert "OUTLIER BIAS REDUCTION CLAIMS" in result.stdout
    assert "Mean ratio across O1-O5: 1.022" in result.stdout
