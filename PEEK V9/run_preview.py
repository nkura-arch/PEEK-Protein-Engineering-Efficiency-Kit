"""Local preview launcher.

Serves the PEEK Flask app with waitress (a small production-grade WSGI server)
instead of the Werkzeug dev server. Waitress supports HTTP keep-alive, which
removes the per-navigation connection stall that made switching tabs feel slow.
"""
from app import app

if __name__ == "__main__":
    try:
        from waitress import serve
        serve(app, host="127.0.0.1", port=5000, threads=8, channel_timeout=120)
    except ImportError:
        # Fallback to the dev server if waitress isn't installed.
        app.run(host="127.0.0.1", port=5000, debug=False, use_reloader=False, threaded=True)
