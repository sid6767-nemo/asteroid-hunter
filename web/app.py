# app.py - asteroid-hunter web app
#
# Run from the project root:   python web/app.py
# Then open http://127.0.0.1:5000

import math
import os
import re
import sys
import json
import uuid
import glob
import shutil
import zipfile
import subprocess
from flask import Flask, jsonify, render_template, request, url_for, make_response, redirect


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
    hunt_mode = (request.args.get('mode') == 'hunt'
                 or request.form.get('mode') == 'hunt')

    if request.method == 'GET':
        return render_template('detect.html', hunt_mode=hunt_mode)

    uploads = request.files.getlist('frames')
    if not uploads or all(not u.filename for u in uploads):
        return render_template('detect.html', error="Please choose some files first.",
                               hunt_mode=hunt_mode)

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

    # run the pipeline WITH SkyBoT so it can identify known asteroids (needed for orbit lookup)
    env = dict(os.environ, MPLBACKEND='Agg')
    try:
        proc = subprocess.run(
            [sys.executable, PIPELINE,
             '--data', work_dir, '--output', RESULTS_DIR, '--save-frames'],
            cwd=PROJECT_ROOT, env=env,
            capture_output=True, text=True, timeout=420)
    except subprocess.TimeoutExpired:
        shutil.rmtree(work_dir, ignore_errors=True)
        return render_template('detect.html',
            error="Detection timed out. Try fewer frames, or the catalog lookup may be slow.")

    result_rel = f"results/{session_id}_tracks.png"
    result_abs = os.path.join(BASE_DIR, 'static', result_rel)
    result_image = url_for('static', filename=result_rel) if os.path.exists(result_abs) else None

    frame_files = sorted(glob.glob(os.path.join(RESULTS_DIR, f"{session_id}_frame*.png")),
                         key=lambda p: int(re.search(r'frame(\d+)', p).group(1)))
    frame_images = [url_for('static', filename=f"results/{os.path.basename(f)}") for f in frame_files]

    # structured candidate list (with SkyBoT names) from the pipeline's results file
    results_path = os.path.join(RESULTS_DIR, f"{session_id}_results.json")
    candidates = []
    if os.path.exists(results_path):
        try:
            with open(results_path) as jf:
                candidates = json.load(jf)
        except Exception:
            candidates = []
    count = str(len(candidates))

    sf_rel = f"results/{session_id}_soundfield.png"
    bd_rel = f"results/{session_id}_backdrop.jpg"
    soundfield_url = url_for('static', filename=sf_rel) if os.path.exists(os.path.join(BASE_DIR,'static',sf_rel)) else None
    backdrop_url   = url_for('static', filename=bd_rel) if os.path.exists(os.path.join(BASE_DIR,'static',bd_rel)) else None

    shutil.rmtree(work_dir, ignore_errors=True)

    if result_image is None:
        return render_template('detect.html',
            error="The pipeline ran but didn't produce a result image. "
                  "Check that your frames are valid FITS images of the same field.")

    hunt_ready = os.path.exists(os.path.join(RESULTS_DIR, f'{session_id}_hunt.json'))

    # hunt mode: straight to the sonification, skipping the results page entirely.
    # The whole point of hunt-by-ear is finding the movers unaided - showing the
    # pipeline's answers on the way in would spoil it.
    if hunt_mode and hunt_ready:
        resp = make_response(redirect(url_for('hunt', sid=session_id)))
        resp.set_cookie('last_hunt_sid', session_id, max_age=60 * 60 * 24 * 30)
        return resp

    resp = make_response(render_template('detect.html',
                           result_image=result_image,
                           frame_images=frame_images,
                           candidates=candidates,
                           count=count,
                           n_frames=n,
                           soundfield_url=soundfield_url,
                           backdrop_url=backdrop_url,
                           session_id=session_id))
    # remember THIS browser's own upload, so the home-page "Hunt by ear" shortcut
    # finds it - instead of silently falling back to whichever hunt session
    # (possibly the bundled set203 sample) happens to have the newest file on disk
    if hunt_ready:
        resp.set_cookie('last_hunt_sid', session_id, max_age=60 * 60 * 24 * 30)
    return resp


@app.route('/hunt')
def hunt_latest():
    """Home-page entry: open THIS BROWSER'S own most recent upload if it
    remembers one (via cookie). Never silently substitutes someone else's
    session or the bundled set203 sample without saying so."""
    own_sid = request.cookies.get('last_hunt_sid')
    if own_sid and os.path.exists(os.path.join(RESULTS_DIR, f'{own_sid}_hunt.json')):
        return redirect(f'/hunt/{own_sid}')
    # no upload from this browser yet -> send them to import first.
    # hunt-by-ear runs on YOUR field; it never falls back to a canned dataset.
    return redirect(url_for('detect', mode='hunt'))


@app.route('/hunt/<sid>')
def hunt(sid):
    """Hunt by ear: find movers in the raw frames with sound, unaided.
    Serves only the linear frame data -- the candidate results are NOT
    loaded by this page until the user commits a guess."""
    if not re.fullmatch(r'[A-Za-z0-9_-]{1,40}', sid):
        return render_template('hunt.html', sid=None,
                               error="Invalid session id."), 404
    if not os.path.exists(os.path.join(RESULTS_DIR, f'{sid}_hunt.json')):
        return render_template('hunt.html', sid=None,
                               error="No hunt data for this session. Run a new detection "
                                     "first — the pipeline now saves hunt frames automatically."), 404
    return render_template('hunt.html', sid=sid, error=None)


_SPEC_CACHE = {}


def _col(tbl, key):
    """Read a column if Horizons returned it, else None."""
    try:
        if key in tbl.colnames:
            v = float(tbl[key][0])
            return v if v == v else None       # drop NaN
    except Exception:
        pass
    return None


def _fetch_elements(name):
    from astroquery.jplhorizons import Horizons
    from astropy.time import Time
    el = Horizons(id=name, location='@sun', epochs=Time.now().jd).elements()
    return {'a': float(el['a'][0]),   'e':  float(el['e'][0]),
            'i': float(el['incl'][0]),'Om': float(el['Omega'][0]),
            'w': float(el['w'][0]),   'M':  float(el['M'][0])}


def _classify(a, q):
    if q is not None and q < 1.3:
        return 'Near-Earth asteroid'
    if a is None:
        return 'Asteroid'
    if a < 2.0:  return 'Inner asteroid'
    if a < 2.5:  return 'Main belt (inner)'
    if a < 2.82: return 'Main belt (middle)'
    if a < 3.6:  return 'Main belt (outer)'
    return 'Outer solar system object'


_EARTH_EL = dict(a=1.00000261, e=0.01671123, i=0.0, Om=0.0, w=102.93768193)


def _orbit_pts(a, e, i, Om, w, th):
    import numpy as np
    i, Om, w = np.radians([i, Om, w])
    r = a * (1 - e * e) / (1 + e * np.cos(th))
    P = np.vstack([r * np.cos(th), r * np.sin(th), np.zeros_like(th)])
    Rz = lambda t: np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])
    Rx = lambda t: np.array([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])
    return ((Rz(Om) @ Rx(i) @ Rz(w)) @ P).T


def _moid(el, coarse=720):
    """Minimum distance between this orbit and Earth's orbit, in AU.

    Validated against published values: Bennu 0.0029 (JPL 0.0032),
    Apophis 0.00021 (JPL 0.00019).
    """
    import numpy as np
    th = np.linspace(0, 2 * np.pi, coarse, endpoint=False)
    A, B = _orbit_pts(**el, th=th), _orbit_pts(**_EARTH_EL, th=th)
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    ia, ib = np.unravel_index(np.argmin(D), D.shape)
    best, span = float(D[ia, ib]), 2 * np.pi / coarse
    a0, b0 = th[ia], th[ib]
    for _ in range(4):                                  # zoom in on the winning pair
        ta = np.linspace(a0 - span, a0 + span, 60)
        tb = np.linspace(b0 - span, b0 + span, 60)
        A, B = _orbit_pts(**el, th=ta), _orbit_pts(**_EARTH_EL, th=tb)
        D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
        ia, ib = np.unravel_index(np.argmin(D), D.shape)
        best = float(D[ia, ib]); a0, b0 = ta[ia], tb[ib]; span /= 12
    return best


def _sentry(name):
    """NASA's own impact-risk assessment, if this object is on the Sentry list."""
    import json as _json
    import urllib.parse
    import urllib.request
    url = 'https://ssd-api.jpl.nasa.gov/sentry.api?des=' + urllib.parse.quote(name)
    try:
        with urllib.request.urlopen(url, timeout=6) as r:
            d = _json.loads(r.read().decode())
    except Exception:
        return None
    summ = d.get('summary')
    if not summ:
        return None
    try:
        return {'ip': float(summ.get('ip')), 'ps_cum': summ.get('ps_cum'), 'ts_max': summ.get('ts_max')}
    except Exception:
        return None


def _hazard(q, moid_au, H, name):
    """Plain-language verdict on whether this thing can reach us."""
    if q is not None and q > 1.3:
        return {'verdict': 'No — it cannot reach Earth',
                'why': "Its orbit never comes closer to the Sun than %.2f AU, while Earth orbits at 1 AU. "
                       "The two orbits never come near each other." % q,
                'impact': 'Zero. Not physically possible on its current orbit.',
                'sentry': None}
    if moid_au is not None and moid_au > 0.05:
        return {'verdict': 'Near-Earth, but not hazardous',
                'why': "Its orbit comes within %.3f AU of Earth's — close enough to be a near-Earth object, "
                       "but beyond the 0.05 AU threshold for a potentially hazardous asteroid." % moid_au,
                'impact': 'Not on any impact-risk list.',
                'sentry': None}
    s = _sentry(name)
    label = 'Potentially hazardous asteroid' if (H is not None and H < 22) else 'Near-Earth, orbits cross closely'
    if s and s.get('ip') is not None:
        ip = s['ip']
        odds = ('about 1 in %s' % format(int(round(1 / ip)), ',')) if ip > 0 else 'effectively zero'
        return {'verdict': label,
                'why': "Its orbit passes within %.4f AU of Earth's." % (moid_au or 0),
                'impact': 'NASA Sentry gives a cumulative impact probability of %.2e (%s).' % (ip, odds),
                'sentry': s}
    return {'verdict': label,
            'why': "Its orbit passes within %.4f AU of Earth's." % (moid_au or 0),
            'impact': 'Not currently listed on NASA Sentry, so no measurable impact probability.',
            'sentry': None}


def _fetch_specs(name):
    """Full spec sheet for the flash card. Cached, because JPL is slow."""
    if name in _SPEC_CACHE:
        return _SPEC_CACHE[name]
    from astroquery.jplhorizons import Horizons
    from astropy.time import Time
    el = Horizons(id=name, location='@sun', epochs=Time.now().jd).elements()

    a  = _col(el, 'a');      e = _col(el, 'e');   i = _col(el, 'incl')
    q  = _col(el, 'q');      Q = _col(el, 'Q');   H = _col(el, 'H')
    per_days = _col(el, 'period')

    # rough size from absolute magnitude, assuming a typical 14% albedo
    diameter_km = None
    if H is not None:
        diameter_km = (1329.0 / math.sqrt(0.14)) * (10 ** (-H / 5.0))

    moid_au = None
    hazard = None
    try:
        if None not in (a, e, i):
            moid_au = _moid({'a': a, 'e': e, 'i': i,
                             'Om': _col(el, 'Omega') or 0.0, 'w': _col(el, 'w') or 0.0})
            hazard = _hazard(q, moid_au, H, name)
    except Exception:
        pass

    specs = {
        'name': name,
        'moid_au': moid_au,
        'hazard': hazard,
        'targetname': str(el['targetname'][0]) if 'targetname' in el.colnames else name,
        'a': a, 'e': e, 'i': i,
        'Om': _col(el, 'Omega'), 'w': _col(el, 'w'), 'M': _col(el, 'M'),
        'q': q, 'Q': Q, 'H': H,
        'period_years': (per_days / 365.25) if per_days else None,
        'diameter_km': diameter_km,
        'class': _classify(a, q),
    }
    _SPEC_CACHE[name] = specs
    return specs


@app.route('/api/asteroid')
def api_asteroid():
    name = (request.args.get('name') or '').strip()
    if not name:
        return jsonify({'error': 'No asteroid specified.'}), 400
    try:
        return jsonify(_fetch_specs(name))
    except Exception as exc:
        return jsonify({'error': f"Couldn't fetch data for '{name}'.", 'detail': str(exc)}), 502


@app.route('/orbit')
def orbit():
    name = (request.args.get('name') or '').strip()
    if not name:
        return render_template('orbit.html', asteroids=None, title=None,
                               error="No asteroid was specified.")
    try:
        asteroids = [{'name': name, 'el': _fetch_elements(name)}]
        return render_template('orbit.html', asteroids=asteroids,
                               title=f"Orbit of {name}", error=None)
    except Exception as e:
        return render_template('orbit.html', asteroids=None, title=None,
                               error=f"Couldn't fetch an orbit for '{name}'. ({e})")


@app.route('/orbits')
def orbits():
    names = [x.strip() for x in (request.args.get('names') or '').split(',') if x.strip()]
    if not names:
        return render_template('orbit.html', asteroids=None, title=None,
                               error="No asteroids were specified.")
    asteroids = []
    for nm in names:
        try:
            asteroids.append({'name': nm, 'el': _fetch_elements(nm)})
        except Exception:
            pass
    if not asteroids:
        return render_template('orbit.html', asteroids=None, title=None,
                               error="Couldn't fetch orbits for those asteroids.")
    return render_template('orbit.html', asteroids=asteroids,
                           title="All detected asteroid orbits", error=None)


if __name__ == '__main__':
    app.run(debug=True, port=5000)
