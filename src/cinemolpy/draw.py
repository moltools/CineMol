import os 

os.environ["PYTHONNET_RUNTIME"] = "coreclr"
runtime = os.environ.get("PYTHONNET_RUNTIME")
import clr 

abs_path = os.path.dirname(os.path.abspath(__file__))
name_dll = "CineMol.dll"
clr.AddReference(os.path.join(abs_path, name_dll))

import CineMol 

def test() -> int:

    
    
    return 0