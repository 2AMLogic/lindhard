"""`Run.run` releases the GIL while it transports."""

import threading
import time

import lindhard as lh
from conftest import EXAMPLES


def test_other_python_threads_progress_during_a_run():
    run = lh.Run.from_toml_file(EXAMPLES / "ar_1keV_cu.toml")
    run.threads = 1
    box = {}
    worker = threading.Thread(target=lambda: box.update(res=run.run(ions=3000)))
    worker.start()
    ticks = 0
    t0 = time.perf_counter()
    while worker.is_alive():
        time.sleep(0.001)
        ticks += 1
    worker.join()
    elapsed = time.perf_counter() - t0
    assert box["res"].histories == 3000
    # With the GIL held for the whole run this thread would tick about once.
    assert elapsed > 0.05, "run too short to tell; raise the ion count"
    assert ticks >= 10
