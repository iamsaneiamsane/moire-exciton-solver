import webview
import threading
import sys
import os
from backend import app

def start_server():
    app.run(port=5000)

def main():
    t = threading.Thread(target=start_server)
    t.daemon = True
    t.start()

    if getattr(sys, 'frozen', False):
        base_path = sys._MEIPASS
    else:
        base_path = os.path.dirname(os.path.abspath(__file__))
    
    file_path = os.path.join(base_path, "webpage.html")
    url = f"file://{file_path}"
    
    icon_path = os.path.join(base_path, "app_icon.ico")
    
    webview.create_window("Moiré Qubit Lab", url, width=1280, height=800)
    webview.start(icon=icon_path)

if __name__ == "__main__":
    main()
