import os
from app import app

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    try:
        from waitress import serve
        print(f"Serving PEEK on http://127.0.0.1:{port}  (waitress)")
        serve(app, host="0.0.0.0", port=port, threads=8, channel_timeout=120)
    except ImportError:
        app.run(host="0.0.0.0", port=port, debug=False, use_reloader=False, threaded=True)
