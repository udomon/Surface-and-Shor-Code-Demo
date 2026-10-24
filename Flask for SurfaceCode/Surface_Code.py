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
##make the graph of the surface code, with the data qubits and the ancilla qubits in their respective positions
    def __init__(self, d):
        if d < 2:
            raise ValueError("d must be greater than or equal to 2")
        self.d = d
        self.data_q, self.z_anc, self.x_anc, self.qr = self.layout(d)
        self.define_stabilisers()
        
        self.circuit = QuantumCircuit(self.qr)

    def get_dictionaries(self):
        return self.data_q, self.z_anc, self.x_anc
    
    def get_stabilisers(self):
        return self.x_stabilisers, self.z_stabilisers
    
    def get_registers(self):
        return self.qr

    def layout(self,d):
        #store each qubit as its coordinates in the grid
        #key is the coordinate, value is the qubit
        data_q = {}
        z_anc = {}
        x_anc = {}

        coords = []
        for i in range(2*d-1):
            for j in range(2*d-1):
                coords.append((i,j)) 
        
        #make a register for each qubit
        qr = QuantumRegister(len(coords),'q')
        qubit_coords = {}

        #goes thru each coordinate pair and assigns a qubit that is in the qr register
        for idx, (i,j) in enumerate(coords):
            qubit_coords[(i,j)] = qr[idx]
      
    #assigning the different types of qubits their respective spots 
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
        # Each dictionary maps a coordinate of the X stabilise/X measuremnt ancilla to the 4 adjacent data qubits
        
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
                    i, j = coord
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

    Returns:
    --------
    - X syndrome extraction circuit
    - Z syndrome extraction circuit
    - Circuit depth and total operations for both circuits
    - No. of parallel CNOT layers
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
    """To inject errors into the surface code,
    an object must be instantiated with distance d, and then the inject_errors method can be called
    with the coordinates of the qubits to be affected by X or Z errors.
    
    To visualise, call visualise_errors with show_syndromes = True to see the triggered syndromes."""
    
    
    def __init__(self, d):
        super().__init__(d)
        self.qc = QuantumCircuit(self.qr)
        
        self.x_errors = set()
        self.z_errors = set()
        self.y_errors = set()
        self.x_syndrome = {}
        self.z_syndrome = {}
        
    def inject_errors(self, x_pos = None, z_pos = None):
        if x_pos:
            for pos in x_pos:
                if pos in self.data_q:
                    self.x_errors.add(pos)
                    
        if z_pos:
            for pos in z_pos:
                if pos in self.data_q:
                    self.z_errors.add(pos)
        
        if x_pos and z_pos:
            for pos in x_pos:
                if pos in z_pos and pos in self.data_q:
                    self.y_errors.add(pos)
                    self.x_errors.remove(pos)
                    self.z_errors.remove(pos)
      
        self._calculate_syndromes()
                    
    def _calculate_syndromes(self):
        #detect phase flips
        for ancilla_pos, neighbour_pos in self.x_stabilisers.items():
            syndrome = 0
            for neighbour in neighbour_pos :
                if neighbour in self.z_errors or neighbour in self.y_errors:
                    syndrome ^= 1 # flip the syndrome bit if a z error is detected
            self.x_syndrome[ancilla_pos] = syndrome
        
        #detect bit flips    
        for ancilla_pos, neighbour_pos in self.z_stabilisers.items():
            syndrome = 0
            for neighbour in neighbour_pos:
                if neighbour in self.x_errors or neighbour in self.y_errors:
                    syndrome ^= 1 # flip the syndrome bit if a x error is detected
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
        
        #adding in the errors to the graph        
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
                colour = 'yellowgreen'   #normal ancilla colour for the plain layout
            else:
                colour = 'white'
            
            element = Rectangle((j-0.3/2,i-0.3/2), 0.3, 0.3, color = colour, ec = 'black', linewidth = 2)
            ax.add_patch(element)
            ax.text(j, i, 'X', ha = 'center', va = 'center', fontweight = 'bold', fontsize = 10)
            
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
                colour = 'lightgreen'   #normal ancilla colour for the plain layout
            else:
                colour = 'white'
            
            element = Rectangle((j-0.3/2,i-0.3/2), 0.3, 0.3, color = colour, ec = 'black', linewidth = 2)
            ax.add_patch(element)
            ax.text(j, i, 'Z', ha = 'center', va = 'center', fontweight = 'bold', fontsize = 10)
            
            if syndrome_type == 'Z' and coord in self.z_syndrome:
                syndrome_val = self.z_syndrome[coord]
                ax.text(j+0.35, i+0.35, str(syndrome_val), ha='center', va='center', 
                       fontweight='bold', fontsize=12, color='purple' if syndrome_val else 'green',
                       bbox=dict(facecolor='black', alpha=0.2, edgecolor='white'))   
            
        if syndrome_type != 'X':
            self._add_elements('connection', 'grey', '--', self.x_stabilisers.items(), ax)
        if syndrome_type != 'Z':
            self._add_elements('connection', 'grey', '-', self.z_stabilisers.items(), ax)

class Error_Correction(Error_Injection):
    """This class is designed to implement the error correction procedure for the surface code.
    It inherits from Error_Injection, allowing it to access the error injection and syndrome calculation methods.
    The class provides methods to perform error correction based on the calculated syndromes."""
    
    def __init__(self, d):
        super().__init__(d)
        self.data_q, self.z_anc, self.x_anc, self.qr = self.layout(d)
        self.qc = QuantumCircuit(self.qr)
        self.syndrome_extractor = SyndromeExtraction(d)
        self.x_syndrome, self.z_syndrome = self._calculate_syndromes()
        
    def calculate_syndromes_from_circuit(self):
        combined_x = self.x_errors | self.y_errors
        combined_z = self.z_errors | self.y_errors
        
        x_bits = self.syndrome_extractor.create_x_syndrome_circuit(
            x_pos=combined_x, z_pos=combined_z)
        z_bits = self.syndrome_extractor.create_z_syndrome_circuit(
            x_pos=combined_x, z_pos=combined_z)
        
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
                
        G = nx.Graph()
        size = 2*self.d-1
        
        boundary_ids = []
        for error_coord in errors:
            G.add_node(error_coord)
            
            boundary_id = ('boundary', error_coord)
            boundary_ids.append(boundary_id)
            
            i, j = error_coord
            if syndrome_type == 'X': #x goes left to right
                boundary_distance = min(j, size - 1 - j) + 0.5
            elif syndrome_type == 'Z': #z goes up and down
                boundary_distance = min(i, size - 1 - i) + 0.5
            G.add_edge(error_coord, boundary_id, weight=-boundary_distance)
                 
        #pair them up
        for i in range(len(errors)):
            for j in range(i + 1, len(errors)):
                G.add_edge(errors[i], errors[j], weight=-self._manhattan_distance(errors[i], errors[j]))

        for i in range(len(boundary_ids)):
            for j in range(i + 1, len(boundary_ids)):
                G.add_edge(boundary_ids[i], boundary_ids[j], weight=0) # so all nodes are matched

        self.matches = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)
        self._touched = set()
        self._generate_coordinate_pairs(self.matches, syndrome_type)
        
        for coord in self._touched:
            self._apply_correction(coord, syndrome_type)
        
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
            
            self._flip_bits(value_1, value_2, syndrome_type)
    
    def _get_boundary_coord(self, other_value, syndrome_type):
        i_2, j_2 = other_value
        size = 2 * self.d - 1
        
        if syndrome_type == 'X':
            j_1 = j_2 #in same column
            if (size - 1 - i_2) < i_2:
                i_1 = size - 1
            else:
                i_1 = 0
        elif syndrome_type == 'Z':
            i_1 = i_2
            if (size - 1 - j_2) < j_2:
                j_1 = size - 1
            else:
                j_1 = 0
        
        return (i_1, j_1)
            
    def _flip_bits(self, value_1, value_2, syndrome_type):
            path_cells = set() #done so paths do not cross and repeat
        
            #traverse x axis first
            left, right = self._sort_pairs(value_1, value_2, 'j') #sort them inb ascending order
            
            row = left[0]
            for j in range(left[1], right[1] + 1):
                coord = (row,j)
                if coord in self.data_q:
                    path_cells.add(coord) #if present, remove, if not present, add
                        
            #traverse y axis
            top, bottom = self._sort_pairs(left, right, 'i') #sort them inb ascending order
            col = right[1]
            
            for i in range(top[0], bottom[0] + 1):
                coord = (i, col)
                if coord in self.data_q:
                    path_cells.add(coord)
                    
            self._touched.symmetric_difference_update(path_cells)
    
    def _apply_correction(self, coord, syndrome_type):
        if coord not in self.data_q:
            return
        qubit = self.data_q[coord]
        if syndrome_type == 'X':
            if coord in self.z_errors:
                self.z_errors.discard(coord)
            elif coord in self.y_errors:
                self.y_errors.discard(coord)
                self.x_errors.add(coord)
            else:
                self.z_errors.add(coord)
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
                
                  
    def _sort_pairs(self, val_1, val_2, basis):
        if basis == 'i':
            op_1, op_2 = val_1[0], val_2[0]
        else:
            op_1, op_2 = val_1[1], val_2[1]
            
        if op_1 > op_2:
            return val_2, val_1
        else:
            return val_1, val_2
            
    def visualise_correction(self, resync = 'quantum'):
        #calculate syndrome each time
        #show path taken -> needs pairs
        #show flipped qubits
        self.visualise_errors()
    
    def Execute(self, x_errors, z_errors, resync = 'quantum', visual = False):      
        if visual: self.visualise_layout()
        
        #validate the inputs
        for error_list in (x_errors or [], z_errors or []): #iterate through the lists together
            for error in error_list:
                if error not in self.data_q:
                    raise ValueError(f"{error} is not a data qubit, error cannot be inserted")
        
        self.inject_errors(x_errors, z_errors)
        
        if visual: self.visualise_errors()
        
        if resync == 'quantum':
            self.x_syndrome, self.z_syndrome = self.calculate_syndromes_from_circuit()
        
        self.minimum_weight_matching('X', self.x_syndrome, resync = None)
        self.minimum_weight_matching('Z', self.z_syndrome, resync = None)
        
        if resync == 'quantum':
            self.x_syndrome, self.z_syndrome = self.calculate_syndromes_from_circuit()
            
        if visual: self.visualise_correction()        
                

if __name__ == "__main__":
    d = 3
    error_correction = Error_Correction(d)
    error_correction.Execute(x_errors=[(1, 3), (1, 1), (2, 0), (3, 1), (4, 2), (3, 3), (2, 4)],
    z_errors=[(2, 2), (2, 0), (3, 1), (4, 2), (3, 3), (2, 4)], visual=True)
    
    
