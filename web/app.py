# app.py - asteroid-hunter web app
#
# Run from the project root:   python web/app.py
# Then open http://127.0.0.1:5000

import os
import re
import sys
import uuid
import glob
import shutil
import zipfile
import subprocess
from flask import Flask, render_template, request, url_for

app = Flask(__name__)

BASE_DIR     = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(BASE_DIR)
PIPELINE     = os.path.join(PROJECT_ROOT, 'scripts', 'exp_set203_pipeline.py')
UPLOAD_ROOT  = os.path.join(BASE_DIR, 'uploads')
RESULTS_DIR  = os.path.join(BASE_DIR, 'static', 'results')
os.makedirs(UPLOAD_ROOT, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

MIN_FRAMES = 3
MAX_FRAMES = 10


@app.route('/')
def home():
    return render_template('index.html')


@app.route('/detect', methods=['GET', 'POST'])
def detect():
    if request.method == 'GET':
        return render_template('detect.html')

    uploads = request.files.getlist('frames')
    if not uploads or all(not u.filename for u in uploads):
        return render_template('detect.html', error="Please choose some files first.")

    session_id = uuid.uuid4().hex[:8]
    work_dir = os.path.join(UPLOAD_ROOT, session_id)
    os.makedirs(work_dir, exist_ok=True)

    for u in uploads:
        if not u.filename:
            continue
        fname = os.path.basename(u.filename)
        dest = os.path.join(work_dir, fname)
        u.save(dest)
        if fname.lower().endswith('.zip'):
            try:
                with zipfile.ZipFile(dest) as z:
                    z.extractall(work_dir)
            except Exception:
                pass
            os.remove(dest)

    fits_files = glob.glob(os.path.join(work_dir, '**', '*.fits'), recursive=True)
    n = len(fits_files)
    if n < MIN_FRAMES:
        shutil.rmtree(work_dir, ignore_errors=True)
        return render_template('detect.html',
            error=f"Found {n} FITS frame(s). You need at least {MIN_FRAMES} "
                  f"(4 or more gives the cleanest results).")

    for f in fits_files:
        if os.path.dirname(f) != work_dir:
            shutil.move(f, os.path.join(work_dir, os.path.basename(f)))

    env = dict(os.environ, MPLBACKEND='Agg')
    try:
        proc = subprocess.run(
            [sys.executable, PIPELINE,
             '--data', work_dir, '--output', RESULTS_DIR,
             '--no-skybot', '--save-frames'],
            cwd=PROJECT_ROOT, env=env,
            capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        shutil.rmtree(work_dir, ignore_errors=True)
        return render_template('detect.html',
            error="Detection timed out (took over 5 minutes). Try fewer frames.")

    stdout = proc.stdout or ""

    result_rel = f"results/{session_id}_tracks.png"
    result_abs = os.path.join(BASE_DIR, 'static', result_rel)
    result_image = url_for('static', filename=result_rel) if os.path.exists(result_abs) else None

    # individual aligned frames for the blink viewer (sorted frame1, frame2, ...)
    frame_files = sorted(glob.glob(os.path.join(RESULTS_DIR, f"{session_id}_frame*.png")),
                         key=lambda p: int(re.search(r'frame(\d+)', p).group(1)))
    frame_images = [url_for('static', filename=f"results/{os.path.basename(f)}") for f in frame_files]

    candidates = re.findall(r'CONFIRMED #\d+: .*', stdout)
    m = re.search(r'(\d+) confirmed candidate', stdout)
    count = m.group(1) if m else (str(len(candidates)) if candidates else "0")

    shutil.rmtree(work_dir, ignore_errors=True)

    if result_image is None:
        return render_template('detect.html',
            error="The pipeline ran but didn't produce a result image. "
                  "Check that your frames are valid FITS images of the same field.")

    return render_template('detect.html',
                           result_image=result_image,
                           frame_images=frame_images,
                           candidates=candidates,
                           count=count,
                           n_frames=n)


if __name__ == '__main__':
    app.run(debug=True, port=5000)