/* asteroid_card.js — a flash card of specs for one asteroid.
 *
 * Used by both the blink viewer and the 3D orbit viewer.
 *   AsteroidCard.show('2002 GE56', { rate: 818, ra: 171.275, dec: 3.264 })
 *   AsteroidCard.toast('Pause first, then click an asteroid.')
 *
 * The rock drawn on the card is a *representative* asteroid, not a photograph:
 * its shape is generated from the asteroid's name, so each object gets its own
 * distinct lump, but nobody should mistake it for the real thing.
 */
(function () {
  'use strict';

  /* ---------- deterministic randomness, seeded from the name ---------- */
  function hashString(s) {
    let h = 2166136261;
    for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); }
    return h >>> 0;
  }
  function mulberry32(seed) {
    return function () {
      seed |= 0; seed = (seed + 0x6D2B79F5) | 0;
      let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /* ---------- build a lumpy rock: sphere vertices pushed outward only ---------- */
  function buildRock(name, radius) {
    const rnd = mulberry32(hashString(name));
    const NLAT = 22, NLON = 30;

    const bumps = [];
    for (let b = 0; b < 7; b++) {
      let x = rnd() * 2 - 1, y = rnd() * 2 - 1, z = rnd() * 2 - 1;
      const L = Math.hypot(x, y, z) || 1;
      bumps.push([x / L, y / L, z / L]);
    }

    const verts = [];
    for (let i = 0; i <= NLAT; i++) {
      const th = Math.PI * i / NLAT;
      for (let j = 0; j < NLON; j++) {
        const ph = 2 * Math.PI * j / NLON;
        const dx = Math.sin(th) * Math.cos(ph), dy = Math.cos(th), dz = Math.sin(th) * Math.sin(ph);
        let bump = 0;
        for (const b of bumps) {
          const d = dx * b[0] + dy * b[1] + dz * b[2];
          if (d > 0) bump += d * d;
        }
        const fine = (Math.sin(dx * 9) * Math.cos(dy * 11) + Math.sin(dz * 13)) * 0.05;
        const r = radius * (0.90 + 0.18 * bump + fine);   // never inward -> stays solid
        verts.push([dx * r, dy * r, dz * r]);
      }
    }

    const faces = [];
    const at = (i, j) => i * NLON + (j % NLON);
    for (let i = 0; i < NLAT; i++) {
      for (let j = 0; j < NLON; j++) {
        faces.push([at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)]);
      }
    }
    return { verts, faces };
  }

  /* ---------- draw it, rotating, into a canvas ---------- */
  function spinRock(canvas, name) {
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;
    const rock = buildRock(name, Math.min(W, H) * 0.36);
    const light = (() => { const v = [-0.4, 0.6, 1]; const L = Math.hypot(v[0], v[1], v[2]); return v.map(c => c / L); })();
    let angle = 0, raf = 0;

    function frame() {
      ctx.clearRect(0, 0, W, H);
      const ca = Math.cos(angle), sa = Math.sin(angle);
      const rot = rock.verts.map(v => [v[0] * ca + v[2] * sa, v[1], -v[0] * sa + v[2] * ca]);

      const drawable = [];
      for (const f of rock.faces) {
        const a = rot[f[0]], b = rot[f[1]], c = rot[f[2]];
        const u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        const w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        const L = Math.hypot(n[0], n[1], n[2]) || 1;
        n = [n[0] / L, n[1] / L, n[2] / L];
        if (n[2] <= 0) continue;                                   // back face, skip
        const lit = Math.max(0.14, n[0] * light[0] + n[1] * light[1] + n[2] * light[2]);
        const z = (a[2] + b[2] + c[2] + rot[f[3]][2]) / 4;
        drawable.push({ f, lit, z });
      }
      drawable.sort((p, q) => p.z - q.z);

      const cx = W / 2, cy = H / 2;
      for (const d of drawable) {
        const s = d.lit;
        const r = Math.round(148 * s), g = Math.round(138 * s), b = Math.round(124 * s);
        ctx.fillStyle = 'rgb(' + r + ',' + g + ',' + b + ')';
        ctx.strokeStyle = ctx.fillStyle;
        ctx.lineWidth = 0.6;
        ctx.beginPath();
        d.f.forEach((vi, k) => {
          const v = rot[vi], x = cx + v[0], y = cy - v[1];
          k ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
        });
        ctx.closePath(); ctx.fill(); ctx.stroke();
      }

      angle += 0.008;
      raf = requestAnimationFrame(frame);
    }
    frame();
    return () => cancelAnimationFrame(raf);
  }

  /* ---------- toast: a hint that fades after 5 seconds ---------- */
  let toastTimer = null;
  function toast(msg) {
    let t = document.getElementById('ac-toast');
    if (!t) {
      t = document.createElement('div');
      t.id = 'ac-toast';
      t.className = 'ac-toast';
      t.setAttribute('role', 'status');
      document.body.appendChild(t);
    }
    t.textContent = msg;
    t.classList.add('show');
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => t.classList.remove('show'), 5000);
  }

  /* ---------- the sliding panel ---------- */
  let stopSpin = null;

  /* The tiny handful of asteroids humanity has actually photographed up close. */
  const VISITED = {
    '433 Eros': 'NEAR Shoemaker', '101955 Bennu': 'OSIRIS-REx', '162173 Ryugu': 'Hayabusa2',
    '25143 Itokawa': 'Hayabusa', '4 Vesta': 'Dawn', '1 Ceres': 'Dawn', '243 Ida': 'Galileo',
    '951 Gaspra': 'Galileo', '253 Mathilde': 'NEAR Shoemaker', '21 Lutetia': 'Rosetta',
    '2867 Steins': 'Rosetta', '65803 Didymos': 'DART', '152830 Dinkinesh': 'Lucy',
    '99942 Apophis': 'OSIRIS-APEX, arriving 2029'
  };
  function visitedBy(specs) {
    const hay = ((specs.targetname || '') + ' ' + specs.name).toLowerCase();
    for (const k in VISITED) {
      const bits = k.split(' ');
      if (hay.indexOf(bits[bits.length - 1].toLowerCase()) !== -1) return { rock: k, mission: VISITED[k] };
    }
    return null;
  }

  function escHandler(e) { if (e.key === 'Escape') close(); }

  function ensureOverlay() {
    let ov = document.getElementById('ac-overlay');
    if (ov) return ov;
    ov = document.createElement('div');
    ov.id = 'ac-overlay';
    ov.className = 'ac-overlay';
    ov.addEventListener('click', e => { if (e.target === ov) close(); });
    document.body.appendChild(ov);
    requestAnimationFrame(() => ov.classList.add('ac-in'));
    document.addEventListener('keydown', escHandler);
    return ov;
  }

  function close() {
    const ov = document.getElementById('ac-overlay');
    if (!ov) return;
    document.removeEventListener('keydown', escHandler);
    if (stopSpin) { stopSpin(); stopSpin = null; }
    ov.classList.remove('ac-in');
    ov.removeAttribute('id');
    setTimeout(() => ov.remove(), 340);
  }

  function num(v, digits, unit) {
    if (v === null || v === undefined || Number.isNaN(v)) return '\u2014';
    return v.toFixed(digits) + (unit ? ' ' + unit : '');
  }

  function hazardClass(verdict) {
    if (!verdict) return 'ac-watch';
    if (verdict.indexOf('cannot reach') !== -1) return 'ac-safe';
    if (verdict.indexOf('not hazardous') !== -1) return 'ac-watch';
    return 'ac-danger';
  }

  function render(specs, extra) {
    const ov = ensureOverlay();

    const rows = [];
    const push = (k, v, note) => rows.push(
      '<div class="ac-row"><span class="ac-k">' + k + '</span>' +
      '<span class="ac-v">' + v + (note ? ' <em class="ac-note">' + note + '</em>' : '') + '</span></div>');

    push('Orbit size (a)', num(specs.a, 3, 'AU'));
    push('Eccentricity (e)', num(specs.e, 3));
    push('Inclination (i)', num(specs.i, 2, '&deg;'));
    push('Orbital period', num(specs.period_years, 2, 'years'));
    push('Closest to Sun', num(specs.q, 3, 'AU'));
    push('Farthest from Sun', num(specs.Q, 3, 'AU'));
    push('Absolute magnitude (H)', num(specs.H, 2));
    push('Diameter', num(specs.diameter_km, 1, 'km'), 'estimated, assumes 14% albedo');

    if (specs.moid_au !== null && specs.moid_au !== undefined) {
      const km = specs.moid_au * 149.598;
      push('Nearest approach to Earth\u2019s orbit',
           num(specs.moid_au, specs.moid_au < 0.01 ? 5 : 3, 'AU'),
           km.toFixed(km < 10 ? 2 : 0) + ' million km');
    }

    if (extra) {
      if (extra.rate !== undefined) push('Speed across the sky', Math.round(extra.rate) + ' &Prime;/day');
      if (extra.ra !== undefined)   push('Right ascension', num(extra.ra, 4, '&deg;'));
      if (extra.dec !== undefined)  push('Declination', num(extra.dec, 4, '&deg;'));
    }

    let hazardHtml = '';
    if (specs.hazard) {
      const h = specs.hazard;
      hazardHtml =
        '<div class="ac-hazard ' + hazardClass(h.verdict) + '">' +
          '<div class="ac-hz-title">Can it hit Earth?</div>' +
          '<div class="ac-verdict">' + h.verdict + '</div>' +
          '<p class="ac-why">' + h.why + '</p>' +
          '<p class="ac-impact">' + h.impact + '</p>' +
        '</div>';
    }

    const vis = visitedBy(specs);
    const imgNote = vis
      ? 'This asteroid <b>has</b> been photographed &mdash; by ' + vis.mission + '. The rock above is still a ' +
        'generated stand-in, not that photograph.'
      : 'No photograph of this asteroid exists. Only about twenty asteroids have ever been visited by ' +
        'spacecraft, and a few hundred more imaged by radar. Every other one \u2014 including this one \u2014 has ' +
        'never been seen as anything but a moving point of light. The rock above is a generated stand-in, ' +
        'shaped from this object\u2019s name.';

    ov.innerHTML =
      '<div class="ac-card" role="dialog" aria-modal="true" aria-label="Asteroid details">' +
        '<button class="ac-close" aria-label="Close">&times;</button>' +
        '<canvas class="ac-rock" width="150" height="150" aria-hidden="true"></canvas>' +
        '<h2 class="ac-name">' + specs.name + '</h2>' +
        '<div class="ac-class">' + specs.class + '</div>' +
        hazardHtml +
        '<div class="ac-rows">' + rows.join('') + '</div>' +
        '<p class="ac-src">' + imgNote + '</p>' +
        '<p class="ac-src">Orbit data from NASA JPL Horizons. Nearest-approach distance computed from the two ' +
          'orbits; impact risk, where it exists, is NASA\u2019s own Sentry assessment.</p>' +
      '</div>';

    stopSpin = spinRock(ov.querySelector('.ac-rock'), specs.name);
    ov.querySelector('.ac-close').onclick = close;
    ov.querySelector('.ac-close').focus();
  }

  function loading(name) {
    const ov = ensureOverlay();
    ov.innerHTML = '<div class="ac-card"><p class="ac-loading">Fetching ' + name + ' from JPL&hellip;</p></div>';
  }

  function show(name, extra) {
    loading(name);
    fetch('/api/asteroid?name=' + encodeURIComponent(name))
      .then(r => r.json())
      .then(d => {
        if (d.error) { close(); toast(d.error); return; }
        render(d, extra);
      })
      .catch(() => { close(); toast("Couldn't reach JPL just now. Try again."); });
  }

  window.AsteroidCard = { show, close, toast };
})();
