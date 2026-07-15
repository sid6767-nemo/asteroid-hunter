/* hunt.js - shift-and-stack asteroid hunting by ear.
 *
 * Every number that reaches the speakers is computed HERE, from the raw
 * pixel data + the user's guess. This file never loads the pipeline's
 * results until the user commits a guess (see commit(), which is the only
 * place the results JSON is fetched).
 *
 * The math is a line-for-line port of scripts/exp_hunt_score.py (Milestone 0),
 * which validated it on set203: worst asteroid > 5x best star/random on both
 * the roam and tune channels.
 *
 * ---- THE SCORE ENGINE (pure functions, console-testable) ----------------
 * After loading, window.huntDebug exposes:
 *   huntDebug.scoreAt(x, y, speedArcsecDay, angleDeg) -> full careful score
 *   huntDebug.roamAt(x, y)  -> {score, speed, angle} two-stage verified roam
 * Coordinates are BINNED pixels (the .bin frames are 2x2-binned).
 */

(function () {
  'use strict';

  /* ================= score engine ================= */

  const PATCH = 9;               // patch half-size: samples u,v in [-9, +9]
  const CORE_R = 2.5, ANN_IN = 6, ANN_OUT = 9, CONC_R = 7;
  const SPEED_MIN = 70, SPEED_MAX = 1250;   // arcsec/day; slower than 70 is
                                            // indistinguishable from a star in ~1h
  const S_REF = 0.8;             // score -> audio loudness scale (~1/3 of the
                                 // faintest validated asteroid's score)

  // patch geometry, precomputed once
  const coreOff = [], annOff = [], patchOff = [], patchIsCore = [], patchInConc = [];
  for (let v = -PATCH; v <= PATCH; v++) {
    for (let u = -PATCH; u <= PATCH; u++) {
      const r = Math.hypot(u, v);
      patchOff.push([u, v]);
      patchIsCore.push(r <= CORE_R);
      patchInConc.push(r <= CONC_R);
      if (r <= CORE_R) coreOff.push([u, v]);
      if (r > ANN_IN && r <= ANN_OUT) annOff.push([u, v]);
    }
  }
  const N_CORE = coreOff.length, N_PATCH = patchOff.length;

  // engine state, filled by init()
  let W = 0, H = 0, NF = 0, T = [], SIG = [], DIFF = [], PXSCALE = 0.256, BINF = 2;
  let BANK_VX = null, BANK_VY = null, BANK_SP = null, BANK_AN = null, NBANK = 0;

  function bilinear(f, x, y) {
    const x0 = Math.floor(x), y0 = Math.floor(y);
    if (x0 < 0 || y0 < 0 || x0 >= W - 1 || y0 >= H - 1) return NaN;
    const i = y0 * W + x0;
    const a = f[i], b = f[i + 1], c = f[i + W], d = f[i + W + 1];
    if (isNaN(a) || isNaN(b) || isNaN(c) || isNaN(d)) return NaN;
    const fx = x - x0, fy = y - y0;
    return a * (1 - fx) * (1 - fy) + b * fx * (1 - fy) + c * (1 - fx) * fy + d * fx * fy;
  }

  function median(arr) {
    const a = arr.filter(v => !isNaN(v)).sort((p, q) => p - q);
    if (!a.length) return NaN;
    const m = a.length >> 1;
    return a.length % 2 ? a[m] : 0.5 * (a[m - 1] + a[m]);
  }

  function velFrom(speed, angleDeg) {           // arcsec/day + deg -> binned px/day
    const v = speed / (PXSCALE * BINF), a = angleDeg * Math.PI / 180;
    return [v * Math.cos(a), v * Math.sin(a)];
  }

  /* Giant-star rejection zones: same idea as the main pipeline's bright-star
   * exclusion (see exp_set203_pipeline.py: SATURATION_LEVEL + giant_zones).
   * Bright stars leave big residuals after alignment (diffraction spikes,
   * imperfect subtraction), and roam's ~8640-velocity search can find a
   * velocity where those residuals happen to line up, scoring artifact
   * "detections" higher than real faint asteroids. The pipeline masks these
   * zones out; the hunt score engine did not - this brings it in line.
   * Threshold and radius are the pipeline's own formula (img_median + 200*sigma,
   * radius = max(180, 10*sqrt(size))), scaled from unbinned to binned pixels. */
  let GIANT_ZONES = [];
  function buildGiantZones(STATIC, W, H, binf) {
    const vals = [];
    for (let i = 0; i < STATIC.length; i++) if (!isNaN(STATIC[i])) vals.push(STATIC[i]);
    const med = median(vals);
    const mad = median(vals.map(v => Math.abs(v - med)));
    const sigma = Math.max(1e-6, 1.4826 * mad);
    const threshold = med + 200 * sigma;              // pipeline's SATURATION_LEVEL

    const visited = new Uint8Array(STATIC.length);
    const AREA_MIN = Math.max(1, Math.ceil(150 / (binf * binf)));   // pipeline: 150 (unbinned px)
    GIANT_ZONES = [];
    for (let start = 0; start < STATIC.length; start++) {
      if (visited[start] || !(STATIC[start] > threshold)) continue;
      // flood-fill this bright blob (4-connected) to find its full extent
      const stack = [start]; visited[start] = 1;
      let sx = 0, sy = 0, count = 0;
      while (stack.length) {
        const i = stack.pop();
        const px = i % W, py = (i / W) | 0;
        sx += px; sy += py; count++;
        const nbrs = [i - 1, i + 1, i - W, i + W];
        for (const ni of nbrs) {
          if (ni < 0 || ni >= STATIC.length || visited[ni]) continue;
          if (Math.abs((ni % W) - px) > 1) continue;    // guard row-wrap on i-1/i+1
          if (STATIC[ni] > threshold) { visited[ni] = 1; stack.push(ni); }
        }
      }
      if (count >= AREA_MIN) {
        const cx = sx / count, cy = sy / count;
        const r = Math.max(180 / binf, 10 * Math.sqrt(count));    // pipeline's radius formula
        GIANT_ZONES.push([cx, cy, r]);
      }
    }
  }
  function inGiantZone(x, y) {
    for (const [cx, cy, r] of GIANT_ZONES)
      if (Math.hypot(x - cx, y - cy) < r) return true;
    return false;
  }

  /* min-over-frames core signal along a track, in LOCAL noise units, with the
   * consistency factor (asteroid brightness is ~constant over an hour; junk
   * tracks stitched from unrelated residuals are wildly inconsistent). */
  function rawS(x, y, vxDay, vyDay) {
    if (inGiantZone(x, y)) return { s: 0, patches: null };   // bright-star zone: refuse
    const cs = [], patches = [];
    for (let k = 0; k < NF; k++) {
      const cx = x + vxDay * T[k], cy = y + vyDay * T[k];
      const p = new Float32Array(N_PATCH);
      let coreBad = false;
      for (let j = 0; j < N_PATCH; j++) {
        p[j] = bilinear(DIFF[k], cx + patchOff[j][0], cy + patchOff[j][1]);
        if (patchIsCore[j] && isNaN(p[j])) coreBad = true;
      }
      if (coreBad) continue;                     // off-image or saturated: drop frame
      const ann = annOff.map(o => bilinear(DIFF[k], cx + o[0], cy + o[1]));
      const b = median(ann);
      const mad = median(ann.map(v => Math.abs(v - b)));
      const noise = Math.max(1.4826 * mad, 0.6 * SIG[k]);   // local, floored
      let sum = 0;
      for (let j = 0; j < N_PATCH; j++) {
        if (patchIsCore[j]) sum += p[j] - b;
        p[j] -= b;
      }
      cs.push(sum / (noise * N_CORE));
      patches.push(p);
    }
    if (cs.length < 3) return { s: 0, patches: null };       // never invent data
    let s = Math.max(0, Math.min.apply(null, cs));
    if (s > 0) {
      const mean = cs.reduce((a, b) => a + b, 0) / cs.length;
      const sd = Math.sqrt(cs.reduce((a, b) => a + (b - mean) * (b - mean), 0) / cs.length);
      s *= Math.max(0, 1 - sd / Math.max(mean, 1e-9));
    }
    return { s: s, patches: patches };
  }

  /* full careful score: local noise + background + zero-velocity veto +
   * concentration + consistency. Returns {score, conc}. */
  function scoreFull(x, y, vxDay, vyDay, sStatic) {
    if (sStatic === undefined) sStatic = rawS(x, y, 0, 0).s;
    const r = rawS(x, y, vxDay, vyDay);
    const s = Math.max(0, r.s - sStatic);        // static junk explains itself
    if (s === 0 || !r.patches) return { score: 0, conc: 0 };
    // median-stack the patches, measure how point-like the pile-up is
    let coreSum = 0, outSum = 0;
    const tmp = new Float32Array(r.patches.length);
    for (let j = 0; j < N_PATCH; j++) {
      if (!patchInConc[j]) continue;
      let n = 0;
      for (let k = 0; k < r.patches.length; k++) {
        const v = r.patches[k][j];
        if (!isNaN(v)) tmp[n++] = v;
      }
      if (!n) continue;
      const sub = Array.prototype.slice.call(tmp, 0, n).sort((p, q) => p - q);
      const m = n >> 1;
      const med = Math.max(0, n % 2 ? sub[m] : 0.5 * (sub[m - 1] + sub[m]));
      outSum += med;
      if (patchIsCore[j]) coreSum += med;
    }
    const conc = outSum > 0 ? coreSum / outSum : 0;
    return { score: s * conc, conc: conc };
  }

  /* fast matched-filter bank: every velocity on a fine grid (25"/d x 2deg =
   * 8640 tracks), core pixels only, b=0 + global noise (the median
   * subtraction already removed the background; full local treatment happens
   * in the verify stage). */
  function buildBank() {
    const sp = [], an = [];
    for (let s = SPEED_MIN; s <= SPEED_MAX; s += 25) sp.push(s);
    for (let a = 0; a < 360; a += 2) an.push(a);
    NBANK = sp.length * an.length;
    BANK_VX = new Float32Array(NBANK); BANK_VY = new Float32Array(NBANK);
    BANK_SP = new Float32Array(NBANK); BANK_AN = new Float32Array(NBANK);
    let i = 0;
    for (const s of sp) for (const a of an) {
      const v = velFrom(s, a);
      BANK_VX[i] = v[0]; BANK_VY[i] = v[1]; BANK_SP[i] = s; BANK_AN[i] = a; i++;
    }
  }
function bank(x, y) {
    if (inGiantZone(x, y)) return new Float32Array(NBANK);   // bright-star zone: all-zero, skip the sweep
    const minC = new Float32Array(NBANK).fill(Infinity);
    const sum = new Float32Array(NBANK), sum2 = new Float32Array(NBANK);
    const nOK = new Uint8Array(NBANK);
    for (let k = 0; k < NF; k++) {
      const f = DIFF[k], tx = x + 0, ty = y + 0, t = T[k], sg = SIG[k];
      for (let i = 0; i < NBANK; i++) {
        const cx = tx + BANK_VX[i] * t, cy = ty + BANK_VY[i] * t;
        let s = 0, bad = false;
        for (let j = 0; j < N_CORE; j++) {
          const v = bilinear(f, cx + coreOff[j][0], cy + coreOff[j][1]);
          if (isNaN(v)) { bad = true; break; }
          s += v;
        }
        if (bad) continue;
        const c = s / (N_CORE * sg);
        nOK[i]++; sum[i] += c; sum2[i] += c * c;
        if (c < minC[i]) minC[i] = c;
      }
    }
    const out = new Float32Array(NBANK);
    for (let i = 0; i < NBANK; i++) {
      if (nOK[i] < 3) continue;
      let s = Math.max(0, minC[i]);
      if (s > 0) {
        const mean = sum[i] / nOK[i];
        const sd = Math.sqrt(Math.max(0, sum2[i] / nOK[i] - mean * mean));
        s *= Math.max(0, 1 - sd / Math.max(mean, 1e-9));
      }
      out[i] = s;
    }
    return out;
  }

  /* two-stage roam: the bank NOMINATES the best few tracks, the full careful
   * score VERIFIES them, a small refinement polishes the winner. Only the
   * verified value is played -- junk that fools the quick look cannot fool
   * the full treatment. */
  function roamVerified(x, y) {
    const s = bank(x, y);
    const idx = Array.from(s.keys()).sort((a, b) => s[b] - s[a]).slice(0, 12);
    const sStatic = rawS(x, y, 0, 0).s;
    let best = { score: 0, conc: 0 }, bestSp = 0, bestAn = 0;
    for (const i of idx) {
      if (s[i] <= 0) break;
      const r = scoreFull(x, y, BANK_VX[i], BANK_VY[i], sStatic);
      if (r.score > best.score) { best = r; bestSp = BANK_SP[i]; bestAn = BANK_AN[i]; }
    }
    if (best.score > 0) {
      for (const dsp of [-12.5, -6, 0, 6, 12.5]) {
        for (const dan of [-1, -0.5, 0, 0.5, 1]) {
          const v = velFrom(bestSp + dsp, bestAn + dan);
          const r = scoreFull(x, y, v[0], v[1], sStatic);
          if (r.score > best.score) { best = r; }
        }
      }
    }
    return { score: best.score, conc: best.conc, speed: bestSp, angle: bestAn };
  }

  /* ================= data loading ================= */

  const SID = window.HUNT_SID;
  const base = '/static/results/' + SID;
  let META = null;

  async function loadData(onProgress) {
    META = await (await fetch(base + '_hunt.json')).json();
    W = META.width; H = META.height; NF = META.nframes;
    T = META.t_days; PXSCALE = META.pixel_scale; BINF = META.bin;
    const frames = [];
    for (let k = 0; k < NF; k++) {
      onProgress('Loading frame ' + (k + 1) + ' of ' + NF + '…');
      const buf = await (await fetch(base + '_hunt_f' + k + '.bin')).arrayBuffer();
      const u16 = new Uint16Array(buf);
      const f = new Float32Array(u16.length);
      const g = META.gain[k], off = META.offset;
      for (let i = 0; i < u16.length; i++)
        f[i] = u16[i] === 0 ? NaN : (u16[i] - off) * g;    // 0 = no-data sentinel
      if (META.align_ok && META.align_ok[k] === false)
        f.fill(NaN);                                        // misregistered: unusable
      frames.push(f);
    }
    onProgress('Removing the static sky…');
    await new Promise(r => setTimeout(r, 30));              // let the message paint
    // static sky = per-pixel median across the unshifted frames: made of
    // things that sit still (stars). Subtracting it leaves only movers+noise.
    const n = W * H, vals = new Float32Array(NF);
    DIFF = frames.map(() => new Float32Array(n));
    const STATIC = new Float32Array(n).fill(NaN);   // per-pixel median: the frozen stars
    for (let i = 0; i < n; i++) {
      let m = 0;
      for (let k = 0; k < NF; k++) { const v = frames[k][i]; if (!isNaN(v)) vals[m++] = v; }
      let med = NaN;
      if (m >= 2) {
        const sub = Array.prototype.slice.call(vals, 0, m).sort((p, q) => p - q);
        const h2 = m >> 1;
        med = m % 2 ? sub[h2] : 0.5 * (sub[h2 - 1] + sub[h2]);
      }
      STATIC[i] = med;
      for (let k = 0; k < NF; k++) DIFF[k][i] = frames[k][i] - med;   // NaN propagates
    }
    onProgress('Finding bright stars to exclude…');
    await new Promise(r => setTimeout(r, 30));
    buildGiantZones(STATIC, W, H, BINF);
    // per-frame robust noise of the differenced data (1.4826 * MAD),
    // estimated from a pixel sample for speed
    SIG = [];
    for (let k = 0; k < NF; k++) {
      const sample = [];
      for (let i = 7; i < n; i += 131) { const v = DIFF[k][i]; if (!isNaN(v)) sample.push(v); }
      const med = median(sample);
      SIG.push(Math.max(1e-6, 1.4826 * median(sample.map(v => Math.abs(v - med)))));
    }
    buildBank();
  }

  /* ================= UI + audio ================= */

  const stage = document.getElementById('huntStage');
  const sky = document.getElementById('huntSky');
  const ov = document.getElementById('huntOv'), octx = ov.getContext('2d');
  const status = document.getElementById('huntStatus');
  const readout = document.getElementById('huntReadout');
  const live = document.getElementById('huntLive');
  const log = document.getElementById('huntLog');
  const startBtn = document.getElementById('huntStart');
  const speechBtn = document.getElementById('huntSpeech');
  const commitBtn = document.getElementById('huntCommit');

  /* speech (same pattern as the detect page) */
  const synth = window.speechSynthesis;
  let speechOn = true, voice = null;
  function pickVoice() {
    if (!synth) return;
    const vs = synth.getVoices();
    const prefer = ['Google UK English Male', 'Google US English', 'Daniel', 'Karen', 'Samantha', 'Alex'];
    for (const p of prefer) { const v = vs.find(x => x.name === p); if (v) { voice = v; return; } }
    voice = vs.find(x => x.lang && x.lang.indexOf('en') === 0) || vs[0] || null;
  }
  if (synth) { pickVoice(); synth.onvoiceschanged = pickVoice; }
  function say(text, interrupt) {
    live.textContent = text;
    if (!speechOn || !synth) return;
    if (interrupt) synth.cancel();
    const u = new SpeechSynthesisUtterance(text);
    if (voice) u.voice = voice;
    u.rate = 0.98;
    synth.speak(u);
  }
  speechBtn.onclick = () => {
    speechOn = !speechOn;
    speechBtn.textContent = 'Speech: ' + (speechOn ? 'on' : 'off');
    speechBtn.setAttribute('aria-pressed', speechOn);
    if (!speechOn && synth) synth.cancel(); else say('Speech on');
  };

  /* state */
  let started = false, mode = 'roam';            // 'roam' | 'tune'
  let cx = 0, cy = 0;                            // cursor, binned px
  let tuneSpeed = 400, tuneAngle = 90;           // manual velocity in tune mode
  let evalOut = { score: 0, conc: 0, speed: 0, angle: 0 };
  let needEval = false, lastEval = 0;
  let cur = 0;                                   // smoothed audio level
  let ac = null, osc = null, gain = null;
  let commits = [];

  function requestEval() { needEval = true; }

  function evaluate() {
    if (mode === 'roam') {
      evalOut = roamVerified(cx, cy);
    } else {
      const v = velFrom(tuneSpeed, tuneAngle);
      const r = scoreFull(cx, cy, v[0], v[1]);
      evalOut = { score: r.score, conc: r.conc, speed: tuneSpeed, angle: tuneAngle };
    }
    readout.textContent =
      mode + '  ·  x=' + cx.toFixed(0) + ' y=' + cy.toFixed(0) +
      '  ·  score=' + evalOut.score.toFixed(3) +
      '  ·  ' + (mode === 'roam' ? 'best track ' : '') +
      evalOut.speed.toFixed(0) + '"/day @ ' + evalOut.angle.toFixed(1) + '°';
  }

  /* canvas overlay */
  function sizeCanvas() { ov.width = sky.clientWidth; ov.height = sky.clientHeight; draw(); }
  sky.addEventListener('load', sizeCanvas);
  window.addEventListener('resize', sizeCanvas);

  function draw() {
    const Wc = ov.width, Hc = ov.height;
    if (!Wc) return;
    octx.clearRect(0, 0, Wc, Hc);
    const x = cx / W * Wc, y = cy / H * Hc;
    for (const c of commits) {
      const px = c.x / W * Wc, py = c.y / H * Hc;
      octx.strokeStyle = '#5ee0a0'; octx.lineWidth = 1.6;
      octx.beginPath(); octx.arc(px, py, 9, 0, 7); octx.stroke();
    }
    octx.strokeStyle = mode === 'roam' ? '#fff' : '#ffd24d';
    octx.globalAlpha = .8; octx.lineWidth = 2;
    octx.beginPath(); octx.arc(x, y, 12, 0, 7); octx.stroke();
    octx.beginPath();
    octx.moveTo(x - 16, y); octx.lineTo(x + 16, y);
    octx.moveTo(x, y - 16); octx.lineTo(x, y + 16);
    octx.stroke();
    if (mode === 'tune') {                        // show the guessed direction
      const a = evalOut.angle * Math.PI / 180;
      octx.beginPath(); octx.moveTo(x, y);
      octx.lineTo(x + 30 * Math.cos(a), y + 30 * Math.sin(a)); octx.stroke();
    }
    octx.globalAlpha = 1;
  }

  /* mouse + keyboard */
  sky.style.pointerEvents = 'auto';
  sky.addEventListener('pointermove', e => {
    if (!started || mode !== 'roam') return;
    const r = sky.getBoundingClientRect();
    cx = (e.clientX - r.left) / r.width * W;
    cy = (e.clientY - r.top) / r.height * H;
    requestEval(); draw();
  });
  sky.addEventListener('pointerdown', e => { stage.focus(); });

  stage.addEventListener('keydown', e => {
    if (!started) return;
    let handled = true;
    const big = e.shiftKey;
    switch (e.key) {
      case 'ArrowLeft':  cx = Math.max(0, cx - (big ? 16 : 3)); break;
      case 'ArrowRight': cx = Math.min(W - 1, cx + (big ? 16 : 3)); break;
      case 'ArrowUp':    cy = Math.max(0, cy - (big ? 16 : 3)); break;
      case 'ArrowDown':  cy = Math.min(H - 1, cy + (big ? 16 : 3)); break;
      case 'Tab':
        mode = mode === 'roam' ? 'tune' : 'roam';
        if (mode === 'tune' && evalOut.speed > 0) {
          tuneSpeed = evalOut.speed; tuneAngle = evalOut.angle;   // seed from pixels
        }
        say(mode === 'roam' ? 'Roam mode' : 'Tune mode. Speed ' +
            tuneSpeed.toFixed(0) + ' arcseconds per day, heading ' +
            tuneAngle.toFixed(0) + ' degrees.', true);
        break;
      case '[': tuneSpeed = Math.max(SPEED_MIN, tuneSpeed - (big ? 5 : 25)); break;
      case ']': tuneSpeed = Math.min(SPEED_MAX, tuneSpeed + (big ? 5 : 25)); break;
      case ';': tuneAngle = (tuneAngle - (big ? 0.5 : 2) + 360) % 360; break;
      case "'": tuneAngle = (tuneAngle + (big ? 0.5 : 2)) % 360; break;
      case 'v': case 'V':
        say('Speed ' + evalOut.speed.toFixed(0) + ' arcseconds per day, heading ' +
            evalOut.angle.toFixed(0) + ' degrees.', true);
        break;
      case 'p': case 'P':
        say('Position ' + cx.toFixed(0) + ', ' + cy.toFixed(0) + ' of ' + W + ' by ' + H + '.', true);
        break;
      case 'Enter': commit(); break;
      default: handled = false;
    }
    if (handled) { e.preventDefault(); requestEval(); draw(); }
  });

  /* commit: the ONLY place the pipeline's answers are consulted.
   * Three-band verdicts based on our measured empty-sky floor (~0.30):
   *   score < FLOOR_MIN     -> noise, refuse to log anything
   *   FLOOR_MIN..FLOOR_MAX  -> "maybe" band, amber commit, cautious wording
   *   score >= FLOOR_MAX    -> real candidate, green commit, confident wording */
  const FLOOR_MIN = 0.35;                    // must beat empty-sky noise floor
  const FLOOR_MAX = 0.45;                    // 2002 QK157 scored 0.46 - real recovery
  async function commit() {
    const g = { x: cx, y: cy, speed: evalOut.speed, angle: evalOut.angle, score: evalOut.score };
    if (g.score < FLOOR_MIN) {
      say('Below the noise floor. Score ' + g.score.toFixed(2) +
          '. Nothing detectable here - keep hunting.', true);
      return;                                // no circle, no log entry
    }
    const band = g.score < FLOOR_MAX ? 'maybe' : 'candidate';

    let cands = null;
    try { cands = await (await fetch(base + '_results.json')).json(); } catch (e) { /* offline is fine */ }
    let verdict;
    if (cands && cands.length) {
      let best = null, bestD = 1e9;
      for (const c of cands) {
        if (!c.fpos) continue;
        const d = Math.hypot(c.fpos[0][0] * W - g.x, c.fpos[0][1] * H - g.y);
        if (d < bestD) { bestD = d; best = c; }
      }
      const speedOK = best && Math.abs(best.rate - g.speed) / best.rate < 0.2;
      if (best && bestD < 6 && speedOK) {
        const who = best.name || ('candidate #' + best.id + ', not in the catalog - possibly undiscovered');
        verdict = 'Confirmed! You detected ' + who + ', moving ' + best.rate +
                  ' arcseconds per day. The pipeline agrees with your track.';
      } else if (best && bestD < 6) {
        verdict = 'You are on a real object, but your speed (' + g.speed.toFixed(0) +
                  ') differs from the pipeline\'s (' + best.rate + '). Tune it and commit again.';
      } else if (band === 'maybe') {
        verdict = 'Low-confidence candidate. Score ' + g.score.toFixed(2) +
                  ' is above the noise floor but below the reliable-detection threshold. ' +
                  'Could be a very faint real object - or an artifact. Blink to verify.';
      } else {
        verdict = 'Candidate at score ' + g.score.toFixed(2) +
                  '. Above threshold and unclaimed by the pipeline - a real recovery candidate. ' +
                  'Blink to verify.';
      }
    } else {
      verdict = 'Committed at score ' + g.score.toFixed(2) + '. No pipeline results to compare against.';
    }
    g.band = band;                           // stored so draw() can color the circle
    commits.push(g);
    const div = document.createElement('div');
    div.className = 'hit hit-' + band;       // CSS colors it by band
    div.textContent = 'x=' + g.x.toFixed(0) + ' y=' + g.y.toFixed(0) + ' · ' +
                      g.speed.toFixed(0) + '"/day @ ' + g.angle.toFixed(1) + '° · score ' +
                      g.score.toFixed(2) + ' — ' + verdict;
    log.prepend(div);
    say(verdict, true);
    draw();
  }
  commitBtn.onclick = commit;

  /* audio: tone = how strongly a real track sings; pulse rate = how sharply
   * the stack focuses (the "focusing a lens" channel) */
  function audioLoop() {
    requestAnimationFrame(audioLoop);
    const now = performance.now();
    if (needEval && now - lastEval > 90) {       // throttle the heavy math
      needEval = false; lastEval = now;
      evaluate();
    }
    if (!ac) return;
    const a = 1 - Math.exp(-evalOut.score / S_REF);
    cur += (a - cur) * 0.2;
    const t = ac.currentTime;
    const pulseRate = 2 + 12 * Math.min(1, evalOut.conc / 0.5);
    const pulse = 0.65 + 0.35 * Math.sin(t * pulseRate * 6.28);
    const duck = (synth && synth.speaking) ? 0.3 : 1.0;
    gain.gain.setTargetAtTime(cur * 0.5 * (cur > 0.04 ? pulse : 1) * duck, t, 0.03);
    osc.frequency.setTargetAtTime(220 + cur * 660, t, 0.03);
  }

  startBtn.onclick = async function () {
    if (started) return;
    startBtn.disabled = true;
    startBtn.textContent = 'Loading…';
    try {
      await loadData(m => { status.textContent = m; });
    } catch (e) {
      status.textContent = 'Could not load hunt data: ' + e;
      startBtn.textContent = 'Failed';
      return;
    }
    ac = new (window.AudioContext || window.webkitAudioContext)();
    osc = ac.createOscillator(); gain = ac.createGain();
    osc.type = 'sine'; osc.frequency.value = 220; gain.gain.value = 0;
    osc.connect(gain); gain.connect(ac.destination); osc.start();
    started = true;
    startBtn.textContent = 'Listening';
    commitBtn.disabled = false;
    cx = W / 2; cy = H / 2;
    sizeCanvas(); stage.focus();
    status.textContent = 'Hunting. Roam with the mouse or arrows; Tab to tune; Enter to commit.';
    say('Listening. Move across the field. When something real is moving under your cursor, ' +
        'you will hear it sing. Nothing is revealed until you commit.', true);
    requestEval();
    audioLoop();
  };

  /* console access for verification (Milestone 2: JS must match Python) */
  window.huntDebug = {
    scoreAt: (x, y, speed, angle) => {
      const v = velFrom(speed, angle);
      return scoreFull(x, y, v[0], v[1]).score;
    },
    roamAt: (x, y) => roamVerified(x, y),
    engine: () => ({ W, H, NF, T, SIG })
  };
})();
