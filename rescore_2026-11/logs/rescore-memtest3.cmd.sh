/home/freiburger/Documents/py_venv/bin/python -c 'b=bytearray(1500*1024*1024); b[::4096]=b"x"*len(b[::4096]); import time; time.sleep(3); print("alloc ok", len(b))'
