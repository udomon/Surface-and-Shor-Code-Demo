from unittest import case

from numpy.linalg import qr
import qiskit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import SamplerV2 as AerSampler
from qiskit import transpile
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import matplotlib.patches as mpatches
import networkx as nx 

class Surface_Code_Graph:
    """
    This maps out the graph of the surface code in a 2d-1 x 2d-1 lattice: 2d-1 as there are d data qubits and d-1 ancilla qubits
    per dimension
    
    layout: parameters - d: dimension of lattice
            variables - data_q: dictionary of all data qubits with their coordinates as the key, and the respective qubit in the circuit as a value
                        z_anc, x_anc: dictionary of ancilla qubits for respective syndrome measurements, same format as data_q
            returns - data_q, z_anc, x_anc: dictionaries of the data qubits, x ancilla qubits and z ancilla qubits
                    - qr: quantum register with all the qubits
    define_stabilizers: - This produces 2 dictionaries for both the X and Z stabilisers. Each dictionary maps each coordinate grid position of a stabilisers 
                        to the 4 adjacent data qubits that stabiliser acts upon, with one dictionary mapping for the X stabilisers and another dictionary mapping for the Z stabilisers.
                        i.e This shows the mapping of the vertex and plaquette stabilizers.
                      - assigns to self.x_stabilisers, and to self.z_stabilisers
    
    add_elements and visualise layout use matplotlib to layout the original grid of the lattice, colour coded for each qubit type
    with connections between stabilisers and data qubits similar to that of the tanner graph
    """
    def __init__(self, d):
        if d < 2:
            raise ValueError("d must be greater than or equal to 2") #otherwise there are no ancilla qubits
        self.d = d
        self.data_q, self.z_anc, self.x_anc, self.qr = self.layout(d)
        self.define_stabilisers()
        
        self.circuit = QuantumCircuit(self.qr)

    def layout(self,d):
        #store each qubit as its coordinates in the grid
        #key is the coordinate, value is the qubit
        data_q = {}
        z_anc = {}
        x_anc = {}

        #make the coordinate grid in the format (column, row)
        coords = []
        for i in range(2*d-1):
            for j in range(2*d-1):
                coords.append((i,j)) 
        
        #make a register for each qubit
        qr = QuantumRegister(len(coords),'q')
        qubit_coords = {}

        #goes through each coordinate pair and assigns a qubit that is in the qr register
        for idx, (i,j) in enumerate(coords):
            qubit_coords[(i,j)] = qr[idx]
      
    #assigning the different types of qubits their respective location on the lattice 
        for (i,j), qubit in qubit_coords.items():
            if (i+j) % 2 == 0:
                data_q[(i,j)] = qubit
            elif i % 2 == 0:
                z_anc[(i,j)] = qubit
            else:
                x_anc[(i,j)] = qubit

        return data_q, z_anc, x_anc, qr

    def define_stabilisers(self):
        # This defines the X and Z stabilisers for the lattice as dictionaries (is NOT ACTUALLY implementing them/ measuring)
        # Each dictionary maps a coordinate of the X/Z stabilise/X/Z measurement ancilla to the 4 adjacent data qubits
        
        self.x_stabilisers = {}
        self.z_stabilisers = {}

        #Each X stabiliser and Z stabiliser respectively correspond to a X or Z measurement qubit

        #Connecting the X measurement qubit with a given X stabiliser

        for ancilla in self.x_anc:
            i, j = ancilla
            neighbours = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]

        # Only include neighbors that are data qubits
            data_neighbours = []
            for n in neighbours:
                if n in self.data_q:
                    data_neighbours.append(n)
            
            if data_neighbours:  # Only add if there are valid neighbors
                self.x_stabilisers[ancilla] = data_neighbours

    #Connecting the Z measurement qubit with a given Z stabiliser

        for ancilla in self.z_anc:
            i, j = ancilla
            neighbours = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]

            # Only include neighbors that are data qubits
            data_neighbours = []
            for n in neighbours:
                if n in self.data_q:
                    data_neighbours.append(n)
            
            if data_neighbours:  # Only add if there are valid neighbors
                self.z_stabilisers[ancilla] = data_neighbours
  
    def _add_elements(self, shape, colour, label, qubit_dict, ax):
        match shape:
            case 'circle':
                for coord in qubit_dict:
                    i, j = coord #the coordinates in the data_q dictionary follow (y,x)
                    element = Circle((j,i), 0.3, color = colour, ec = 'white', linewidth = 2)
                    ax.add_patch(element)
                    ax.text(j, i, label, ha = 'center', va = 'center', fontweight = 'bold', fontsize = 10)
            case 'square':
                for coord in qubit_dict:
                    i, j = coord
                    element = Rectangle((j-0.3/2,i-0.3/2), 0.3, 0.3, color = colour, ec = 'black', linewidth = 2)
                    ax.add_patch(element)
                    ax.text(j, i, label, ha = 'center', va = 'center', fontweight = 'bold', fontsize = 10)
            case 'connection':
                for coord, qubits in qubit_dict:
                    i, j = coord
                    for qi, qj in qubits:
                        ax.plot([j, qj], [i, qi], color =colour, linestyle = label, alpha=0.7, linewidth = 2)   
  
    def visualise_layout(self):
        fig, ax = plt.subplots(1, 1, figsize=(10,8))

        size = 2 * self.d - 1 #dimensions of the grid
        
        #draw data qubits
        self._add_elements('connection', 'grey', '--', self.x_stabilisers.items(),ax)
        self._add_elements('connection', 'grey', '-', self.z_stabilisers.items(),ax)
        self._add_elements('circle', 'navy', 'D', self.data_q,ax)
        self._add_elements('square', 'yellowgreen', 'X', self.x_anc,ax)
        self._add_elements('square', 'lightgreen', 'Z', self.z_anc,ax)
        

        ax.set_xlim(-0.5, size-0.5)
        ax.set_ylim(-0.5, size-0.5)
        ax.set_aspect('equal')
        ax.invert_yaxis() 
        plt.tight_layout()
        plt.show()

class SyndromeExtraction(Surface_Code_Graph):
    """
    Implement syndrome extraction circuits.
    Create separate circuits for X and Z syndrome extraction.

    -apply_errors - x_pos: an array of coordinates where the x errors have been injected in the format (y,x)
                  - z_pos: an array of coordinates where the z errors have been injected in the format (y,x)
                  By taking in the coordinate positions of the errors to be carried out, this function enacts all the errors as either X, Z or Y (X & Z) errors

    - create_x_syndrome_circuit - x_pos, z_pos
                                - This function applies the Hadamard gates to transform to uniform superposition basis where CNOTS are performed between the ancilla qubits and the data qubits to act as 
                                  a parity check, which conveys whether the stabilizer has an eigenvalue of +1 or -1 through whether the parity check returns 0 or 1 respectively.
                                - We first reset the ancilla qubits. to |0>, then we apply Hamards to all data qubits. (This should move us to the eigenstates for the Z error, there is some confusion here!!!)
                                - Then we use the principle HZH = X, and we apply hadamards to the ancilla qubits, apply the cnot layers, and then reapply hadamards to the ancilla qubits
                                - (Q: Why does applying hamards to data qubits but then not back again work?)
    - create_z_syndrome_circuit - x_pos, z_pos
                                - This function work similarly to the function above but implements the Z stabilisers through using a paritty check with the CNOTS
    - organize_x_cnots - Uses a greedy algorithm to separate the CNOTS into different layers and apply the cnots in parallel to speed up implementation of cnots and ensure implementation is efficient

    - organize_z_cnots - Similarily uses a greedy algorithm to the implement CNOT layers.


    - draw_circuit - draws the syndrome extraction quantum circuit
    

    

    Returns:
    --------
    - X syndrome extraction circuit
    - Z syndrome extraction circuit
    - Circuit depth and total operations for both circuits

    """
    
    def __init__(self, distance):
        super().__init__(distance)
        self.x_measure_dict = {}
        self.z_measure_dict = {}
        self.classical_x = ClassicalRegister(len(self.x_anc), 'cx')
        self.classical_z = ClassicalRegister(len(self.z_anc), 'cz')
        self._sim = AerSimulator()


    def apply_errors(self, x_pos = None, z_pos = None):
        x_pos = set(x_pos or [])
        z_pos = set(z_pos or [])
        
        y_pos = x_pos & z_pos
        x_only = x_pos - y_pos
        z_only = z_pos - y_pos
        
        for pos in x_only:
            if pos in self.data_q:
                self.circuit.x(self.data_q[pos])
        for pos in z_only:
            if pos in self.data_q:
                self.circuit.z(self.data_q[pos])
        for pos in y_pos:
            if pos in self.data_q:
                self.circuit.y(self.data_q[pos])
            

    def create_x_syndrome_circuit(self, x_pos = None, z_pos = None):
        """
        Create X syndrome extraction circuit.
        The syndrome extraction circuit is optimized to
        implement conflict-free CNOT gates.
        """

        self.circuit = QuantumCircuit(self.qr, self.classical_x)
       
        for qubit_obj in self.x_anc.values():
            self.circuit.reset(qubit_obj)

        for qubit_obj in self.data_q.values():
            self.circuit.h(qubit_obj)

        self.apply_errors(x_pos, z_pos)
        
        for qubit_obj in self.x_anc.values():
            self.circuit.h(qubit_obj)
            
        self.organize_x_cnots()
        
        for qubit_obj in self.x_anc.values():
            self.circuit.h(qubit_obj)

        # Measurement of X measure qubits
        x_measure_dict = {}
        i = 0
        
        for measure_pos in self.x_anc:
            measure_qubit = self.x_anc[measure_pos]
            
            self.circuit.measure(measure_qubit,self.classical_x[i])
            x_measure_dict[measure_pos] = '0'
            i += 1

        #Returning the measurement of the ancilla qubits in the X syndrome extraction circuit as a BIT STRING
        transpiled = transpile(self.circuit)
        result = self._sim.run(transpiled, shots=1).result()
        counts = result.get_counts()
        bitstring = max(counts, key=counts.get)
        bitstring = bitstring[::-1]  # Reverse to read left-to-right matching index 0
        
        x_measure_dict = {}
        i = 0
        for measure_pos in self.x_anc:
            
            x_measure_dict[measure_pos] = bitstring[i]
            i += 1

        self.x_measure_dict = x_measure_dict
        return self.x_measure_dict
        

    def create_z_syndrome_circuit(self, x_pos, z_pos):
        """
        Create Z syndrome extraction circuit.
        The syndrome extraction circuit is optimized to
        implement conflict-free CNOT gates.
        
        """

        self.circuit = QuantumCircuit(self.qr, self.classical_z)
        for qubit_obj in self.z_anc.values():
            self.circuit.reset(qubit_obj)
            
        self.apply_errors(x_pos, z_pos)
        
        self.organize_z_cnots()

        # Measurement of Z measure qubits
        z_measure_dict = {}

        i = 0
        for measure_pos in self.z_anc:
            measure_qubit = self.z_anc[measure_pos]


            self.circuit.measure(measure_qubit, self.classical_z[i])
            z_measure_dict[measure_pos] = '0'
            i += 1
        
        #Returning the measurement of the ancilla qubits in the X syndrome extraction circuit as a BIT STRING
        transpiled = transpile(self.circuit)
        result = self._sim.run(transpiled, shots=1).result()

        counts = result.get_counts()
        bitstring = max(counts, key = counts.get)
        bitstring = bitstring[::-1]  # Reverse the bitstring to match the order of qubits in the circuit
        
        z_measure_dict = {}
        i = 0
        for measure_pos in self.z_anc:
            z_measure_dict[measure_pos] = bitstring[i]
            i += 1
        self.z_measure_dict = z_measure_dict

        return self.z_measure_dict
     
    def organize_x_cnots(self):
        """
        Organize CNOT gates of the x-stabilizer measurement into conflict-free layers.
        """
        layers = []
        remaining_cnots = []
        
        # Collect all CNOT operations for X stabilizers
        for measure_pos, data_positions in self.x_stabilisers.items():
            measure_qubit = self.x_anc[measure_pos]
            for data_pos in data_positions:
                data_qubit = self.data_q[data_pos]
                remaining_cnots.append([measure_qubit, data_qubit])
        #Storing a list of all the CNOTs that need to be implemented for the X stabilisers in a list called remaining_cnots
       
        # avoid conflicts
        while remaining_cnots:
            current_layer = []
            used_qubits = set()
            
            i = 0
            while i < len(remaining_cnots):
                cnot = remaining_cnots[i]
                control, target = cnot
                
                # Check if any qubit in this CNOT is already used in current layer
                if control not in used_qubits and target not in used_qubits:
                    current_layer.append(cnot)
                    used_qubits.add(control)
                    used_qubits.add(target)
                    remaining_cnots.pop(i)
                else:
                    i += 1
            
            if current_layer:
                layers.append(current_layer)
        
        #The above has created a list of lists called layers, where each inner list contains CNOTs that can be applied in parallel without conflict. This used a greedy algorithm to iteratively build layers of CNOTs, ensuring that no qubit is involved in multiple operations within the same layer.

        #Applying the CNOTs in each layer to the circuit

        for layer in layers:
            for pair in layer:
                control, target = pair
                self.circuit.cx(control, target)
        
    
    def organize_z_cnots(self):
        """
        Organize CNOT gates of the z-stabiliser measurement into conflict-free layers.
        """
        layers = []
        remaining_cnots = []
        
        # Collect all CNOT operations for Z stabilizers
        for measure_pos, data_positions in self.z_stabilisers.items():
            measure_qubit = self.z_anc[measure_pos]
            for data_pos in data_positions:
                data_qubit = self.data_q[data_pos]
                remaining_cnots.append([data_qubit, measure_qubit])
        
        # avoid conflicts
        while remaining_cnots:
            current_layer = []
            used_qubits = set()
            
            i = 0
            while i < len(remaining_cnots):
                cnot = remaining_cnots[i]
                control, target = cnot
                
                # Check if any qubit in this CNOT is already used in current layer
                if control not in used_qubits and target not in used_qubits:
                    current_layer.append(cnot)
                    used_qubits.add(control)
                    used_qubits.add(target)
                    remaining_cnots.pop(i)
                else:
                    i += 1
            
            if current_layer:
                layers.append(current_layer)
        
        #The above has created a list of lists called layers, where each inner list contains CNOTs that can be applied in parallel without conflict. This used a greedy algorithm to iteratively build layers of CNOTs, ensuring that no qubit is involved in multiple operations within the same layer.
        #Applying the CNOTs in each layer to the circuit

        for layer in layers:
            for pair in layer:
                control, target = pair
                self.circuit.cx(control, target)

    def draw_circuit(self):
        """
        Draw the circuit for the syndrome extraction.
        """
        fig = self.circuit.draw(output='mpl')
        
        plt.show()

class Error_Injection(Surface_Code_Graph):
    """
    This injects errors by taking input from the user and adding them onto the lattice - and the quantum circuit
    
    attributes ----
    x_errors, y_errors, z_errors: dictionaries to store respective errors with their coordinates and their qubits
    x_syndrome, z_syndrome: dictionaries to store the syndrome of coordinates in the format coordinate (x,y): syndrome
    
    inject_errors: parameters - x_pos: an array of coordinates where the x errors have been injected in the format (y,x)
                              - y_pos: an array of coordinates where the z errors have been injecte din the format (y,x)
                  returns - x_syndrome, z_syndrome: classical calculations of the syndromes in the patch
    _calculate_syndromes: parameters - none
                          variables - syndrome: initially 0, and flipped if odd number of neighbouring qubits have errors (anti-communte)
    visualise_error, plot_base_layout_errors, plot_layout_errors plots the errors on the original lattice showing the x and z syndromes using matplotlib
                              
    """
    def __init__(self, d):
        super().__init__(d)
        self.qc = QuantumCircuit(self.qr)
        
        self.x_errors = set()
        self.z_errors = set()
        self.y_errors = set()
        self.x_syndrome = {}
        self.z_syndrome = {}
        
    def inject_errors(self, x_pos = None, z_pos = None):
      #adding injected x errors into the x error dictionary  
        if x_pos:
            for pos in x_pos:
                if pos in self.data_q:
                    self.x_errors.add(pos)
      #adding injected z errors into the z error dictionary                
        if z_pos:
            for pos in z_pos:
                if pos in self.data_q:
                    self.z_errors.add(pos)
      #adding errors that appear in both x and z dictionaries as y errors
        if x_pos and z_pos:
            for pos in x_pos:
                if pos in z_pos and pos in self.data_q:
                    self.y_errors.add(pos)
                    self.x_errors.remove(pos)
                    self.z_errors.remove(pos)
      
        self._calculate_syndromes()
                    
    def _calculate_syndromes(self):
        #detect phase flips in neighbouring data qubits of ancilla classically
        for ancilla_pos, neighbour_pos in self.x_stabilisers.items():
            syndrome = 0
            for neighbour in neighbour_pos :
                if neighbour in self.z_errors or neighbour in self.y_errors:
                    syndrome ^= 1 # flip the syndrome bit if a z error is detected with XOR
            self.x_syndrome[ancilla_pos] = syndrome
        
        #detect bit flips in neighbouring data qubits of ancilla   
        for ancilla_pos, neighbour_pos in self.z_stabilisers.items():
            syndrome = 0
            for neighbour in neighbour_pos:
                if neighbour in self.x_errors or neighbour in self.y_errors:
                    syndrome ^= 1 # flip the syndrome bit if a x error is detected with XOR
            self.z_syndrome[ancilla_pos] = syndrome
            
        return self.x_syndrome, self.z_syndrome
    
    def visualise_errors(self, show_syndromes = True):
        if show_syndromes:
            fig, ax = plt.subplots(1, 3, figsize=(21, 7))
        else:
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            ax = [ax]
           
        size = 2 * self.d - 1
        
        self._plot_errors_layout(ax[0], size, show_syndromes = False)
        ax[0].set_title("Errors", fontsize = 14, fontweight = 'bold')
        
        if show_syndromes:
            self._plot_errors_layout(ax[1], size, syndrome_type = 'X')
            ax[1].set_title("X Syndromes", fontsize = 14, fontweight = 'bold')
            
            self._plot_errors_layout(ax[2], size, syndrome_type = 'Z')
            ax[2].set_title("Z Syndromes", fontsize = 14, fontweight = 'bold')
        
        for axis in ax:
            axis.set_xlim(-0.5, size-0.5)
            axis.set_ylim(-0.5, size-0.5)
            axis.set_aspect('equal')
            axis.grid(True, alpha=0.3)

        #for the key in the graph
        legend_elements = [
            Circle((0, 0), 0.1, color='navy', ec='white', label='Data Qubit'),
            Rectangle((0, 0), 0.1, 0.1, color='yellowgreen', ec='black', label='X Measure Qubit'),
            Rectangle((0, 0), 0.1, 0.1, color='lightgreen', ec='black', label='Z Measure Qubit'),
            Circle((0, 0), 0.1, color='red', ec='black', linewidth=3, label='X Error'),
            Circle((0, 0), 0.1, color='blue', ec='black', linewidth=3, label='Z Error'),
            Circle((0, 0), 0.1, color='purple', ec='black', linewidth=3, label='Y Error'),
        ]
        
        if show_syndromes:
            legend_elements.extend([
                mpatches.Patch(color='purple', alpha=0.7, label='Triggered Syndrome'),
                mpatches.Patch(color='white', ec='black', label='Quiet Syndrome'),
            ])
        
        
        fig.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 0.02), ncol=4)
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)
        plt.show()
        
    def _plot_errors_layout(self, ax, size, show_syndromes = True, syndrome_type = None):
        self._plot_base_layout_errors(ax, size, syndrome_type)
        
        error_colours = {'X': 'red', 'Z': 'blue', 'Y': 'purple'}
        
        #adding in the errors to the graph by highlighting in a coloured ring       
        for coord in self.x_errors:
            i,j = coord
            error_circle = Circle((j,i), 0.35, facecolor = 'none', edgecolor = error_colours['X'], linewidth = 5, alpha = 0.9)
            ax.add_patch(error_circle)
            ax.text(j+0.3, i+0.3, 'X', ha='center', va='center', 
                   fontweight='bold', fontsize=10, color= error_colours['X'])
            
        for coord in self.z_errors:
            i,j = coord
            error_circle = Circle((j,i), 0.35, facecolor = 'none', edgecolor = error_colours['Z'], linewidth = 5, alpha = 0.9)
            ax.add_patch(error_circle)
            ax.text(j+0.3, i+0.3, 'Z', ha='center', va='center', 
                   fontweight='bold', fontsize=10, color= error_colours['Z'])
            
        for coord in self.y_errors:
            i,j = coord
            error_circle = Circle((j,i), 0.35, facecolor = 'none', edgecolor = error_colours['Y'], linewidth = 5, alpha = 0.9)
            ax.add_patch(error_circle)
            ax.text(j+0.3, i+0.3, 'Y', ha='center', va='center', 
                   fontweight='bold', fontsize=10, color= error_colours['Y'])
            
    def _plot_base_layout_errors(self, ax, size, syndrome_type = None):
        self._add_elements('circle', 'navy', 'D', self.data_q, ax)
        
        for coord in self.x_anc:
            i, j = coord
            
            if syndrome_type == 'X' and coord in self.x_syndrome:
                colour = 'purple' if self.x_syndrome[coord] == 1 else 'yellowgreen' #change colour depending on if syndrome triggered
            elif syndrome_type is None:
                colour = 'yellowgreen'   #if it is in normal layout
            else:
                colour = 'white' #if show_syndromes = True but for the diagram the syndrome type is ignored
            
            element = Rectangle((j-0.3/2,i-0.3/2), 0.3, 0.3, color = colour, ec = 'black', linewidth = 2)
            ax.add_patch(element)
            ax.text(j, i, 'X', ha = 'center', va = 'center', fontweight = 'bold', fontsize = 10)

          #print 0 or 1 depending on syndrome
            if syndrome_type == 'X' and coord in self.x_syndrome:
                syndrome_val = self.x_syndrome[coord]
                ax.text(j+0.35, i+0.35, str(syndrome_val), ha='center', va='center', 
                       fontweight='bold', fontsize=12, color='purple' if syndrome_val else 'green',
                       bbox=dict(facecolor='black', alpha=0.2, edgecolor='white'))
        
            
        for coord in self.z_anc:
            i, j = coord
            
            if syndrome_type == 'Z' and coord in self.z_syndrome:
                colour = 'purple' if self.z_syndrome[coord] == 1 else 'lightgreen' #change colour depending on if syndrome triggered
            elif syndrome_type is None:
                colour = 'lightgreen'   ##if it is in normal layout
            else:
                colour = 'white' #if show_syndromes = True but for the diagram the syndrome type is ignored
            
            element = Rectangle((j-0.3/2,i-0.3/2), 0.3, 0.3, color = colour, ec = 'black', linewidth = 2)
            ax.add_patch(element)
            ax.text(j, i, 'Z', ha = 'center', va = 'center', fontweight = 'bold', fontsize = 10)

            #print 0 or 1 depending on syndrome
            if syndrome_type == 'Z' and coord in self.z_syndrome:
                syndrome_val = self.z_syndrome[coord]
                ax.text(j+0.35, i+0.35, str(syndrome_val), ha='center', va='center', 
                       fontweight='bold', fontsize=12, color='purple' if syndrome_val else 'green',
                       bbox=dict(facecolor='black', alpha=0.2, edgecolor='white'))   
        #making the grey lines between the qubits    
        if syndrome_type != 'X':
            self._add_elements('connection', 'grey', '--', self.x_stabilisers.items(), ax)
        if syndrome_type != 'Z':
            self._add_elements('connection', 'grey', '-', self.z_stabilisers.items(), ax)

class Error_Correction(Error_Injection):
    """
    This corrects the previously injected errors by connecting pairs of syndrome measurements that equal 1 (via minimum weight perfect matching)
    and then flips the bits along the path between the pair of qubits
    
    attribites:
    qc - quantum circuit
    syndrome_extractor - extracts x and z syndromes based off of the quantum circuit
    x_syndrome, z_syndrome - dictionary of classically extracted syndromes
    
    calculate_syndromes_from_circuits:  parameters - none
                                        variables - combined_x, combined_z: all x / z errors, splitting the y error dictionary into these groups
                                                  - x_bits, z_bits: dictionary of qubits and their syndrome measurements with the format coordinate: bit (as a string)
                                                  - x_syndrome, z_syndrome: dictionary of qubits and their syndrome measurements in the format coordinate: bit (as an integer)
    minimum_weight_matching: parameters - syndrome type: type of syndrome
                                        - syndrome_dict: respective dictionary for the syndromes
                            variables - errors: list of errors for the syndrome type
                                      - boundary_ids: list of boundary_id for each coordinate
                                      - G: graph of all erroneous syndrome qubits
                                      - boundary_distance: distance of coordinate to boundary
                                      - matches: list of pairs with the lowest overall cost
                                      - _touched: list of erroneous coordinates without any duplicates
    manhattan_distance: paramaters - coord1, coord2: two tuples that both hold an x,y value corresponding to a qubit each
                        returns - integer value of distance between points
    generate_coordinate_pairs: parameters - pairs : list of all pairs/matches
                                          - syndrome_type: type of syndrome being corrected
    get_boundary_coord: parameters - other_value: coordinate of qubit in the bulk of the lattice
                                   - syndrome_type
                        variables - i_2, j_2: i and j components of the 'other value' to be compared to
                                  - i_1, j_1: i and j components of the new boundary coordinate
                        returns - (i_1, j_1) the new boundary coordinate
    sort_pairs: parameters - value_1, value_2: coordinates of two qubits
                           - component: if ordering based off of horizontal or vertical component
                returns - value_1, value_2: the coordinate pairs in ascending order 
    flip_bits: parameters - value_1, value_2: the start and end coordinates respectively
                          - syndrome_type
              variables - path_cells: qubits that have been traversed
    apply_correction: parameters - coord: coordinate where correctoion will be applied
                                 - syndrome_type
                                 
                                        
                                                  
      
    """
    
    def __init__(self, d):
        super().__init__(d)
        self.data_q, self.z_anc, self.x_anc, self.qr = self.layout(d)
        self.qc = QuantumCircuit(self.qr)
        self.syndrome_extractor = SyndromeExtraction(d)
        self.x_syndrome, self.z_syndrome = self._calculate_syndromes()
        
    def calculate_syndromes_from_circuit(self):

        #gets the union of  errors from both sets
        combined_x = self.x_errors | self.y_errors
        combined_z = self.z_errors | self.y_errors

      #gets the syndromes for each coordinate as characters
        x_bits = self.syndrome_extractor.create_x_syndrome_circuit(
            x_pos=combined_x, z_pos=combined_z)
        z_bits = self.syndrome_extractor.create_z_syndrome_circuit(
            x_pos=combined_x, z_pos=combined_z)
        #
      #converts each syndrome from being a character to an integer
        self.x_syndrome = {pos: int(bit) for pos, bit in x_bits.items()}
        self.z_syndrome = {pos: int(bit) for pos, bit in z_bits.items()}
        return self.x_syndrome, self.z_syndrome
    
        
    def minimum_weight_matching(self, syndrome_type, syndrome_dict, resync = 'quantum'):
        # based on the provided syndrome type ('X' or 'Z').
        #resync means get syndromes either classically, with the quantum circuit or make no change
        errors = []
        
        #add the errors to the set of errors based on the syndrome dictionary
        for syndrome in syndrome_dict:
            if syndrome_dict[syndrome] == 1:
                errors.append(syndrome)

      #make a graph for the error coordinates
        G = nx.Graph()
        size = 2*self.d-1

      #store all unique identifiers of boundary coordinates for error nodes
        boundary_ids = []
        for error_coord in errors:
            G.add_node(error_coord)

          #make a unique boundary for each coordinate to match to if a pair has a higher cost for overall matching / odd number of 1 syndromes
            boundary_id = ('boundary', error_coord)
            boundary_ids.append(boundary_id)
            
            i, j = error_coord
          #calculate boundary distance for each coordinate
            if syndrome_type == 'X': #x goes left to right
                boundary_distance = min(j, size - 1 - j) + 0.5
            elif syndrome_type == 'Z': #z goes up and down
                boundary_distance = min(i, size - 1 - i) + 0.5
            G.add_edge(error_coord, boundary_id, weight=-boundary_distance)
                 
        #pair the coordinates to each other
        for i in range(len(errors)):
            for j in range(i + 1, len(errors)):
                G.add_edge(errors[i], errors[j], weight=-self._manhattan_distance(errors[i], errors[j]))
              
        #pair the boundaries together so every node in the graph will be paired
        for i in range(len(boundary_ids)):
            for j in range(i + 1, len(boundary_ids)):
                G.add_edge(boundary_ids[i], boundary_ids[j], weight=0) # so all nodes are matched

        self.matches = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)
        self._touched = set() #avoids correcting the same coordinate twice - making an error
        self._generate_coordinate_pairs(self.matches, syndrome_type) #produces actual pairs of coordinates
        
        for coord in self._touched:
            self._apply_correction(coord, syndrome_type)

      #recalculate syndromes if using the circuit
        if resync == 'quantum':
            self.x_syndrome, self.z_syndrome = self.calculate_syndromes_from_circuit()
        elif resync == 'classical': 
            self._calculate_syndromes() 
        else:
            pass               
            
    def _manhattan_distance(self, coord1, coord2):
        # Calculate the Manhattan distance between two coordinates
        return abs(coord1[0] - coord2[0]) + abs(coord1[1] - coord2[1]) / 2 #from one qubit to another should be 1, not 2

    def _generate_coordinate_pairs(self, pairs, syndrome_type):
      #get pairs of coordinates for the correction path, rather than a coordinate and a boundary id
        size = 2 * self.d - 1
        for pair in pairs:
            value_1, value_2 = pair
            i_1, j_1 = value_1
            i_2, j_2 = value_2
            
            #get two pairs of coordinates
            if i_1 == 'boundary' and i_2 == 'boundary':
                continue  # Skip if both are boundary nodes
            elif i_1 == 'boundary':
                value_1 = self._get_boundary_coord(value_2, syndrome_type)
            elif i_2 == 'boundary' :
                value_2 = self._get_boundary_coord(value_1, syndrome_type)

          #finding the path between the pair of coordinates
            self._flip_bits(value_1, value_2, syndrome_type)
    
    def _get_boundary_coord(self, other_value, syndrome_type):
      #generating a coordinate for the boundary id
      #other node is the error node it is paired to
        i_2, j_2 = other_value
        size = 2 * self.d - 1
        
        if syndrome_type == 'X':
            j_1 = j_2 #in same column
            if (size - 1 - i_2) < i_2: #checking if closer to the right or left boundary
                i_1 = size - 1
            else:
                i_1 = 0
        elif syndrome_type == 'Z':
            i_1 = i_2 #in the same row
            if (size - 1 - j_2) < j_2: #checking if closer to the top or bottom boundary
                j_1 = size - 1
            else:
                j_1 = 0
        
        return (i_1, j_1) #the new coordinate
            
    def _flip_bits(self, value_1, value_2, syndrome_type):
            path_cells = set() #done so paths do not cross and repeat
        
            #traverse x axis first
            left, right = self._sort_pairs(value_1, value_2, 'j') #sort them in ascending order
            
            row = left[0] #operate on the same row
            for j in range(left[1], right[1] + 1): #changing column
                coord = (row,j)
                if coord in self.data_q:
                    path_cells.add(coord)
                        
            #traverse y axis
            top, bottom = self._sort_pairs(left, right, 'i') #sort them in ascending order
            col = right[1] #operate on same column
            
            for i in range(top[0], bottom[0] + 1):
                coord = (i, col)
                if coord in self.data_q:
                    path_cells.add(coord)
                    
            self._touched.symmetric_difference_update(path_cells) #if in self touched and path cells, removed, otherwise it is added to touched
    
    def _apply_correction(self, coord, syndrome_type):
      #only flip data qubits
        if coord not in self.data_q:
            return
        qubit = self.data_q[coord]

      #correct the error and modify the respective dictionaries
        if syndrome_type == 'X':
            if coord in self.z_errors:
                self.z_errors.discard(coord) #error resolved so remove from dictionary
            elif coord in self.y_errors:
                self.y_errors.discard(coord) #error resolved so now only an x error
                self.x_errors.add(coord)
            else:
                self.z_errors.add(coord) #flip occured so is now a new 'error'
            self.qc.z(qubit)
            
        elif syndrome_type == 'Z':
            if coord in self.x_errors:
                self.x_errors.discard(coord)
            elif coord in self.y_errors:
                self.y_errors.discard(coord)
                self.z_errors.add(coord)
            else:
                self.x_errors.add(coord)
            self.qc.x(qubit)
                
                  
    def _sort_pairs(self, val_1, val_2, component):
      #if sorting them based on row or column
        if component == 'i':
            op_1, op_2 = val_1[0], val_2[0]
        else:
            op_1, op_2 = val_1[1], val_2[1]
        #if the first coordinate is large, swap the values around    
        if op_1 > op_2:
            return val_2, val_1
        else:
            return val_1, val_2
    
def visualise_correction(self, resync = 'quantum'):
        """ This is an unfinished function that should show pairwise correction of all errors across the patch"""
        #calculate syndrome each time
        #show path taken -> needs pairs
        #show flipped qubits
        self.visualise_errors()
    
def Execute(self, x_errors, z_errors, resync = 'quantum', visual = False): 
        #main function used to execute the program

        #visual can be used for debugging to see the steps taken and check logic
        if visual: self.visualise_layout()
        
        #validate the inputs
        for error_list in (x_errors or [], z_errors or []): #iterate through the lists together
            for error in error_list:
                if error not in self.data_q:
                    raise ValueError(f"{error} is not a data qubit, error cannot be inserted")
        
        self.inject_errors(x_errors, z_errors)
        
        if visual: self.visualise_errors()

  #calculating syndromes from the circuit, rather than classically
        if resync == 'quantum':
            self.x_syndrome, self.z_syndrome = self.calculate_syndromes_from_circuit()
        
        self.minimum_weight_matching('X', self.x_syndrome, resync = None)
        self.minimum_weight_matching('Z', self.z_syndrome, resync = None)

  #calculating syndromes after the correction algorithm
        if resync == 'quantum':
            self.x_syndrome, self.z_syndrome = self.calculate_syndromes_from_circuit()
            
        if visual: self.visualise_correction()        
                

if __name__ == "__main__":
    d = 3
   
    error_correction = Error_Correction(d)
    error_correction.Execute(x_errors= [(3,1), (1,3), (4,4)], z_errors=[(4,4)])
                
    
