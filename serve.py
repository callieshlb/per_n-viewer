#!/usr/bin/env python3
"""Serve this folder on http://127.0.0.1:PORT and open the Per_n(0) atlas. usage: python3 serve.py [PORT]
Keep the window open while you use the site (Ctrl-C stops it). The pages need a server (they fetch their data files)."""
import functools, http.server, os, sys, threading, webbrowser
here = os.path.dirname(os.path.abspath(__file__))
handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=here)
ports = [int(sys.argv[1])] if len(sys.argv) > 1 else list(range(8000, 8030))
httpd = None
for port in ports:
    try:
        httpd = http.server.ThreadingHTTPServer(('127.0.0.1', port), handler); break
    except OSError as e:
        print(f'port {port} not available ({e}); trying the next one')
if httpd is None: sys.exit('no free port; try: python3 serve.py 8123')
url = f'http://127.0.0.1:{port}/index.html'
print(f'serving {here}\nopen {url}   (keep this window open; Ctrl-C to stop)')
threading.Timer(0.7, lambda: webbrowser.open(url)).start()
try: httpd.serve_forever()
except KeyboardInterrupt: print('stopped')
