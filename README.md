# Surface-and-Shor-Code-Demo
This contains files used for an interactive demo of Quantum Error Correction through the Shor and Surface Codes

The Shor_Code.py file contains the ShorCode class that implements encoding, error injection, error detection and correction and decoding using Shor code.

The MQ_Shor_Code.py file contains the MultipleLogicalQubitsShorCode(ShorCode) class that implements Shor code for multiple logical qubits, and implements the transversal gates X, Z and CNOT between logical qubits.

The Surface_Code.py file contains multiple classes to implement Surface code

The MQ_Surface_Code.py file implements the class LogicalOperations(Error_Correction), and implements the transversal gates X, Z and CNOT between logical qubits.


Flask for Shor Code and Flask for SurfaceCode can be downloaded for the visual implementations of the respective codes. In the main file, app.py links the HTML frontend to the Python backend, Shor_Code.py and Surface_Code.py contain the logic for the visualisations, and in the template folder, index is the HTML visual aspect.

Dependencies:

- qiskit
- qiskit-aer
- flask
- networkx
- matplotlib

Surface-and-Shor-Code-Demo/
├── Shor_Code.py
├── MQ_Shor_Code.py
├── Surface_Code.py
├── MQ_Surface_Code.py
├── Flask_ShorCode/
│   ├── app.py
│   └── templates/index.html
└── Flask_SurfaceCode/
    ├── app.py
    └── templates/index.html

To run the visualizers, in the terminal run python app.py for the respecite demos

Port for Shor code demo: http://127.0.0.1:5001

Port for Surface code demo: http://127.0.0.1:5000





