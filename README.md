# Surface-and-Shor-Code-Demo
This contains files used for an interactive demo of Quantum Error Correction through the Shor and Surface Codes

The Shor_Code.py file contains the ShorCode class that implements encoding, error injection, error detection and correction and decoding using Shor code.

The Multiple_Qubit_Shor_Code.py contains the MultipleLogicalQubitsShorCode(ShorCode) class that implements Shor code for multiple logical qubits, and implements the transversal gates X, Z and CNOT between logical qubits.

Flask for Shor Code and Flask for SurfaceCode can be downloaded for the visual implementations of the respective codes. In the main file, app.py links the HTML frontend to the Python backend, Shor_Code.py and Surface_Code.py contain the logic for the visualisations, and in the template folder, index is the HTML visual aspect.


