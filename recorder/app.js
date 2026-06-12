'use strict';

const STORAGE_KEY = 'roaaads_sessions';
const GPS_INTERVAL_MS = 20000;

let sessions = [];
let currentSession = null;
let pendingConfig = null; // { label, mode } — set before GPS phase
let gpsWatchId = null;
let lastGpsFix = null;
let lastTrackTimestamp = 0;
let timerInterval = null;
let selectedMode = 'dual';
let singleDirection = 'with';
let selectedIds = new Set();

// ── Storage ────────────────────────────────────────────────────────────────

function loadSessions() {
  try {
    sessions = JSON.parse(localStorage.getItem(STORAGE_KEY)) || [];
  } catch {
    sessions = [];
  }
}

function saveSessions() {
  localStorage.setItem(STORAGE_KEY, JSON.stringify(sessions));
}

// ── Screen navigation ──────────────────────────────────────────────────────

function showScreen(id) {
  document.querySelectorAll('.screen').forEach(s => s.classList.remove('active'));
  document.getElementById(id).classList.add('active');
}

// ── Home ───────────────────────────────────────────────────────────────────

function renderSessions() {
  const list  = document.getElementById('sessions-list');
  const empty = document.getElementById('no-sessions');
  list.innerHTML = '';

  // Remove any selectedIds that no longer exist
  const existingIds = new Set(sessions.map(s => s.id));
  for (const id of selectedIds) {
    if (!existingIds.has(id)) selectedIds.delete(id);
  }

  if (sessions.length === 0) {
    empty.style.display = 'block';
    updateBatchBar();
    return;
  }
  empty.style.display = 'none';

  [...sessions].reverse().forEach(session => {
    const withCount    = session.events.filter(e => e.direction === 'with').length;
    const againstCount = session.events.filter(e => e.direction === 'against').length;
    const date         = new Date(session.startTime).toLocaleString();
    const modeLabel    = session.mode === 'dual' ? 'Dual' : 'Single';
    const duration     = session.endTime
      ? formatDuration(new Date(session.endTime) - new Date(session.startTime))
      : 'in progress';
    const isChecked = selectedIds.has(session.id) ? 'checked' : '';

    const card = document.createElement('div');
    card.className = 'session-card';
    card.innerHTML = `
      <input type="checkbox" class="session-chk" data-id="${session.id}" ${isChecked}>
      <div class="session-info">
        <div class="session-title">${escapeHtml(session.label)}</div>
        <div class="session-meta">
          ${date} · ${modeLabel} · ${duration}<br>
          With: ${withCount} · Against: ${againstCount} · GPS pts: ${session.gpsTrack.length}
        </div>
      </div>
      <div class="session-actions">
        <button class="btn-sm danger" data-action="del" data-id="${session.id}">&#10005;</button>
      </div>`;
    list.appendChild(card);
  });

  updateBatchBar();
}

function updateBatchBar() {
  const bar       = document.getElementById('batch-bar');
  const exportBtn = document.getElementById('btn-export-selected');
  const selectAll = document.getElementById('chk-select-all');

  if (sessions.length === 0) {
    bar.style.display = 'none';
    return;
  }
  bar.style.display = 'flex';

  const count = selectedIds.size;
  exportBtn.disabled   = count === 0;
  exportBtn.textContent = count > 0 ? `Export CSV (${count})` : 'Export CSV';
  selectAll.checked       = count > 0 && count === sessions.length;
  selectAll.indeterminate = count > 0 && count < sessions.length;
}

document.getElementById('sessions-list').addEventListener('change', e => {
  const chk = e.target.closest('.session-chk');
  if (!chk) return;
  if (chk.checked) selectedIds.add(chk.dataset.id);
  else selectedIds.delete(chk.dataset.id);
  updateBatchBar();
});

document.getElementById('sessions-list').addEventListener('click', e => {
  const btn = e.target.closest('[data-action]');
  if (!btn) return;
  if (btn.dataset.action === 'del') deleteSession(btn.dataset.id);
});

document.getElementById('chk-select-all').addEventListener('change', e => {
  if (e.target.checked) sessions.forEach(s => selectedIds.add(s.id));
  else selectedIds.clear();
  renderSessions();
});

function deleteSession(id) {
  if (!confirm('Delete this session?')) return;
  sessions = sessions.filter(s => s.id !== id);
  selectedIds.delete(id);
  saveSessions();
  renderSessions();
}

function clearAllSessions() {
  if (!confirm(`Delete all ${sessions.length} session${sessions.length !== 1 ? 's' : ''}? This cannot be undone.`)) return;
  sessions = [];
  selectedIds.clear();
  saveSessions();
  renderSessions();
}

// ── Setup ──────────────────────────────────────────────────────────────────

function initSetupScreen() {
  document.getElementById('session-label').value = '';
  selectedMode = 'dual';
  document.querySelectorAll('.mode-btn').forEach(btn => {
    btn.classList.toggle('active', btn.dataset.mode === 'dual');
  });
}

document.querySelectorAll('.mode-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    selectedMode = btn.dataset.mode;
    document.querySelectorAll('.mode-btn').forEach(b =>
      b.classList.toggle('active', b === btn)
    );
  });
});

// ── Recording — GPS phase ──────────────────────────────────────────────────

function goToRecordingScreen() {
  pendingConfig = {
    label: document.getElementById('session-label').value.trim() || new Date().toLocaleString(),
    mode: selectedMode,
  };
  lastTrackTimestamp = 0;

  // If GPS is already running (came from New Session mid-recording), show current state
  if (gpsWatchId !== null) {
    document.getElementById('btn-start-gps').disabled = true;
    document.getElementById('btn-start-gps').textContent = 'GPS Active';
    if (lastGpsFix) {
      setGPSPhaseUI('ready', lastGpsFix.accuracy);
      document.getElementById('btn-start-recording').disabled = false;
    } else {
      setGPSPhaseUI('acquiring');
      document.getElementById('btn-start-recording').disabled = true;
    }
  } else {
    lastGpsFix = null;
    setGPSPhaseUI('idle');
    document.getElementById('btn-start-gps').disabled = false;
    document.getElementById('btn-start-gps').textContent = 'Start GPS';
    document.getElementById('btn-start-recording').disabled = true;
  }

  document.getElementById('rec-label').textContent = pendingConfig.label;
  document.getElementById('rec-timer').textContent = '00:00:00';
  document.getElementById('gps-panel').style.display      = 'flex';
  document.getElementById('counting-panel').style.display = 'none';
  document.getElementById('rec-timer').style.display      = 'none';
  document.getElementById('btn-back-recording').style.display = 'block';

  showScreen('screen-recording');
}

function handleStartGPS() {
  document.getElementById('btn-start-gps').disabled = true;
  document.getElementById('btn-start-gps').textContent = 'GPS Active';
  setGPSPhaseUI('acquiring');
  startGPS();
}

// ── Recording — counting phase ─────────────────────────────────────────────

function handleStartRecording() {
  currentSession = {
    id: crypto.randomUUID(),
    label: pendingConfig.label,
    startTime: new Date().toISOString(),
    endTime: null,
    mode: pendingConfig.mode,
    gpsTrack: [],
    events: [],
  };

  if (lastGpsFix) {
    currentSession.gpsTrack.push({ ...lastGpsFix });
    lastTrackTimestamp = Date.now();
  }

  singleDirection = 'with';
  updateSingleUI();

  const isDual = currentSession.mode === 'dual';
  document.getElementById('dual-ui').style.display   = isDual ? 'flex' : 'none';
  document.getElementById('single-ui').style.display = isDual ? 'none' : 'flex';

  ['count-with-dual', 'count-against-dual', 'count-with-single',
   'count-against-single', 'count-total-single'].forEach(id => {
    document.getElementById(id).textContent = '0';
  });

  document.getElementById('gps-panel').style.display      = 'none';
  document.getElementById('counting-panel').style.display = 'flex';
  document.getElementById('rec-timer').style.display      = 'inline';
  document.getElementById('btn-back-recording').style.display = 'none';

  startTimer();
}

function stopSession() {
  if (!confirm('Stop recording and save session?')) return;
  saveAndEndSession();
  stopGPS();
  renderSessions();
  showScreen('screen-home');
}

function newSessionMidRecording() {
  if (!confirm(`Save "${currentSession.label}" and start a new session?`)) return;
  saveAndEndSession();
  // GPS keeps running — carry lastGpsFix into the next session's GPS phase
  initSetupScreen();
  showScreen('screen-setup');
}

function saveAndEndSession() {
  currentSession.endTime = new Date().toISOString();
  stopTimer();
  sessions.push(currentSession);
  saveSessions();
  currentSession = null;
}

function cancelRecording() {
  stopGPS();
  pendingConfig = null;
  showScreen('screen-setup');
}

function recordTap(direction) {
  if (!currentSession) return;
  currentSession.events.push({
    timestamp: new Date().toISOString(),
    direction,
    directionMode: currentSession.mode === 'dual' ? 'both' : 'single',
    nearestGps: lastGpsFix && isRecentFix(lastGpsFix) ? { ...lastGpsFix } : null,
  });
  updateCounts();
}

function isRecentFix(fix) {
  return (Date.now() - new Date(fix.timestamp).getTime()) < 60000;
}

function updateCounts() {
  const w = currentSession.events.filter(e => e.direction === 'with').length;
  const a = currentSession.events.filter(e => e.direction === 'against').length;
  document.getElementById('count-with-dual').textContent      = w;
  document.getElementById('count-against-dual').textContent   = a;
  document.getElementById('count-with-single').textContent    = w;
  document.getElementById('count-against-single').textContent = a;
  document.getElementById('count-total-single').textContent   = w + a;
}

// ── GPS ────────────────────────────────────────────────────────────────────

function startGPS() {
  setGPSStatus('waiting');

  if (!navigator.geolocation) {
    setGPSStatus('unavailable');
    setGPSPhaseUI('unavailable');
    return;
  }

  gpsWatchId = navigator.geolocation.watchPosition(
    onGPSSuccess,
    onGPSError,
    { enableHighAccuracy: true, maximumAge: 20000, timeout: 30000 }
  );
}

function stopGPS() {
  if (gpsWatchId !== null) {
    navigator.geolocation.clearWatch(gpsWatchId);
    gpsWatchId = null;
  }
}

function onGPSSuccess(pos) {
  const fix = {
    lat: pos.coords.latitude,
    lng: pos.coords.longitude,
    accuracy: Math.round(pos.coords.accuracy),
    timestamp: new Date(pos.timestamp).toISOString(),
  };
  lastGpsFix = fix;
  setGPSStatus('ok', fix.accuracy);

  // GPS phase: update UI and enable Start Recording on first fix
  if (document.getElementById('gps-panel').style.display !== 'none') {
    setGPSPhaseUI('ready', fix.accuracy);
    document.getElementById('btn-start-recording').disabled = false;
  }

  // Counting phase: append to track at interval
  if (currentSession) {
    const now = Date.now();
    if (now - lastTrackTimestamp >= GPS_INTERVAL_MS) {
      currentSession.gpsTrack.push(fix);
      lastTrackTimestamp = now;
    }
  }
}

function onGPSError() {
  setGPSStatus('error');
  if (document.getElementById('gps-panel').style.display !== 'none') {
    setGPSPhaseUI('error');
  }
}

function setGPSStatus(state, accuracy) {
  const el = document.getElementById('gps-status');
  el.className = 'gps-status';
  if (state === 'ok') {
    el.textContent = `GPS ±${accuracy}m`;
    el.classList.add('ok');
  } else if (state === 'error') {
    el.textContent = 'GPS error';
    el.classList.add('err');
  } else if (state === 'unavailable') {
    el.textContent = 'No GPS';
    el.classList.add('err');
  } else {
    el.textContent = 'GPS...';
  }
}

function setGPSPhaseUI(state, accuracy) {
  const icon = document.getElementById('gps-big-icon');
  const text = document.getElementById('gps-big-text');
  const sub  = document.getElementById('gps-big-sub');
  if (state === 'idle') {
    icon.textContent = '📍';
    text.textContent = 'GPS not started';
    sub.textContent  = '';
    text.className   = 'gps-big-text';
  } else if (state === 'acquiring') {
    icon.textContent = '🔄';
    text.textContent = 'Acquiring position…';
    sub.textContent  = 'This can take up to 30 seconds';
    text.className   = 'gps-big-text';
  } else if (state === 'ready') {
    icon.textContent = '✅';
    text.textContent = 'Position locked';
    sub.textContent  = accuracy != null ? `±${accuracy}m accuracy` : '';
    text.className   = 'gps-big-text gps-ready';
  } else if (state === 'error') {
    icon.textContent = '❌';
    text.textContent = 'GPS unavailable';
    sub.textContent  = 'You can still record without location data';
    text.className   = 'gps-big-text gps-err';
    document.getElementById('btn-start-recording').disabled = false;
  } else if (state === 'unavailable') {
    icon.textContent = '❌';
    text.textContent = 'GPS not supported';
    sub.textContent  = 'Requires HTTPS. You can still record without location data.';
    text.className   = 'gps-big-text gps-err';
    document.getElementById('btn-start-recording').disabled = false;
  }
}

// ── Timer ──────────────────────────────────────────────────────────────────

function startTimer() {
  const t0 = Date.now();
  timerInterval = setInterval(() => {
    document.getElementById('rec-timer').textContent = formatDuration(Date.now() - t0);
  }, 1000);
}

function stopTimer() {
  clearInterval(timerInterval);
  timerInterval = null;
}

function formatDuration(ms) {
  const s   = Math.floor(ms / 1000);
  const h   = Math.floor(s / 3600);
  const m   = Math.floor((s % 3600) / 60);
  const sec = s % 60;
  return `${pad(h)}:${pad(m)}:${pad(sec)}`;
}

function pad(n) {
  return String(n).padStart(2, '0');
}

// ── Single mode toggle ─────────────────────────────────────────────────────

function toggleSingleDirection() {
  singleDirection = singleDirection === 'with' ? 'against' : 'with';
  updateSingleUI();
}

function updateSingleUI() {
  const isWith = singleDirection === 'with';
  document.getElementById('btn-toggle-direction').className    = `toggle-btn${isWith ? '' : ' against'}`;
  document.getElementById('toggle-current').textContent        = isWith ? 'WITH ↑' : 'AGAINST ↓';
  document.getElementById('single-tap-arrow').textContent      = isWith ? '↑' : '↓';
  document.getElementById('btn-tap-single').style.background   = isWith ? 'var(--with)' : 'var(--against)';
}

// ── Export ─────────────────────────────────────────────────────────────────

function csvEscape(val) {
  const s = String(val ?? '');
  return s.includes(',') || s.includes('"') || s.includes('\n')
    ? `"${s.replace(/"/g, '""')}"` : s;
}

function exportSelectedCSV() {
  const selected = sessions.filter(s => selectedIds.has(s.id));
  if (selected.length === 0) return;

  const header = [
    'session_id', 'road_name', 'session_mode', 'session_start', 'session_end',
    'timestamp', 'direction', 'directionMode', 'lat', 'lng', 'gps_accuracy_m',
  ];
  const rows = [header];

  for (const session of selected) {
    for (const e of session.events) {
      rows.push([
        session.id.slice(0, 8),
        session.label,
        session.mode,
        session.startTime,
        session.endTime ?? '',
        e.timestamp,
        e.direction,
        e.directionMode,
        e.nearestGps?.lat      ?? '',
        e.nearestGps?.lng      ?? '',
        e.nearestGps?.accuracy ?? '',
      ].map(csvEscape));
    }
  }

  const date = new Date().toISOString().slice(0, 10);
  downloadFile(
    `traffic_${date}.csv`,
    rows.map(r => r.join(',')).join('\n'),
    'text/csv'
  );
}

function downloadFile(filename, content, mimeType) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([content], { type: mimeType }));
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(a.href);
}

// ── Utilities ──────────────────────────────────────────────────────────────

function escapeHtml(str) {
  return String(str).replace(/[&<>"']/g, c =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])
  );
}

// ── Wire up events ─────────────────────────────────────────────────────────

document.getElementById('btn-new-session').addEventListener('click', () => {
  initSetupScreen();
  showScreen('screen-setup');
});

document.getElementById('btn-back-setup').addEventListener('click', () => {
  showScreen('screen-home');
});

document.getElementById('btn-export-selected').addEventListener('click', exportSelectedCSV);
document.getElementById('btn-clear-all').addEventListener('click', clearAllSessions);

document.getElementById('btn-go-to-recording').addEventListener('click', goToRecordingScreen);
document.getElementById('btn-back-recording').addEventListener('click', cancelRecording);

document.getElementById('btn-start-gps').addEventListener('click', handleStartGPS);
document.getElementById('btn-start-recording').addEventListener('click', handleStartRecording);
document.getElementById('btn-stop').addEventListener('click', stopSession);
document.getElementById('btn-new-session-rec').addEventListener('click', newSessionMidRecording);

document.getElementById('btn-with').addEventListener('click', () => recordTap('with'));
document.getElementById('btn-against').addEventListener('click', () => recordTap('against'));

document.getElementById('btn-tap-single').addEventListener('click', () => recordTap(singleDirection));
document.getElementById('btn-toggle-direction').addEventListener('click', toggleSingleDirection);

// ── Boot ───────────────────────────────────────────────────────────────────

loadSessions();
renderSessions();
showScreen('screen-home');
