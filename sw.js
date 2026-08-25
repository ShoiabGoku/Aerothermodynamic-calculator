/* Service worker for the Aerodynamics Calculator PWA.
   Strategy: network-first for index.html (so pushes propagate immediately),
   cache-first for static assets (icons, manifest). Bump VERSION on releases
   that change cached assets. */
const VERSION = 'aerocalc-v5';
const ASSETS = [
  './',
  './index.html',
  './manifest.json',
  './icon-192.png',
  './icon-512.png',
  './icon-512-maskable.png'
];

self.addEventListener('install', e => {
  e.waitUntil(caches.open(VERSION).then(c => c.addAll(ASSETS)).then(() => self.skipWaiting()));
});

self.addEventListener('activate', e => {
  e.waitUntil(
    caches.keys()
      .then(keys => Promise.all(keys.filter(k => k !== VERSION).map(k => caches.delete(k))))
      .then(() => self.clients.claim())
  );
});

self.addEventListener('fetch', e => {
  if (e.request.method !== 'GET') return;
  const url = new URL(e.request.url);
  if (url.origin !== location.origin) return;

  const isHTML = e.request.mode === 'navigate' || url.pathname.endsWith('index.html') || url.pathname.endsWith('/');
  if (isHTML) {
    // network-first: always try for the freshest app, fall back to cache offline
    e.respondWith(
      fetch(e.request)
        .then(r => { const copy = r.clone(); caches.open(VERSION).then(c => c.put(e.request, copy)); return r; })
        .catch(() => caches.match(e.request).then(r => r || caches.match('./index.html')))
    );
  } else {
    // cache-first for static assets
    e.respondWith(
      caches.match(e.request).then(r => r || fetch(e.request).then(resp => {
        const copy = resp.clone(); caches.open(VERSION).then(c => c.put(e.request, copy)); return resp;
      }))
    );
  }
});
