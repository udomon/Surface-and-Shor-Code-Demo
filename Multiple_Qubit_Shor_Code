import qiskit
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import SamplerV2 as AerSampler
from qiskit import transpile








class ShorCode:


  '''
  
  This class is almost identical to Shor_Code.py ShorCode class, however there has been some restructuring of the quantum circuit 
  so that the MultipleLogicalQubitsShorCode(ShorCode) class creates one entire quantum circuit for all the logical qubits.
  
  '''
    def __init__(self, index=0, external_qc=None):
        tag = f"L{index}_"  # avoid name collisions once multiple ShorCodes share a circuit
        self.b1_qr = QuantumRegister(3, tag+'XQ1')
        self.b1_ar = QuantumRegister(2, tag+'XA1')
        self.b1_cr = ClassicalRegister(2, tag+'XC1')

        self.b2_qr = QuantumRegister(3, tag+'XQ2')
        self.b2_ar = QuantumRegister(2, tag+'XA2')
        self.b2_cr = ClassicalRegister(2, tag+'XC2')

        self.b3_qr = QuantumRegister(3, tag+'XQ3')
        self.b3_ar = QuantumRegister(2, tag+'XA3')
        self.b3_cr = ClassicalRegister(2, tag+'XC3')

        self.phase_ar = QuantumRegister(2, tag+'ZA')
        self.phase_cr = ClassicalRegister(2, tag+'ZC')
        self.final_cr = ClassicalRegister(1, tag+'final_output')

        self.qubit_regs = [self.b1_qr, self.b2_qr, self.b3_qr]
        self.classical_regs = [self.b1_cr, self.b2_cr, self.b3_cr]
        self.ancilla_regs = [self.b1_ar, self.b2_ar, self.b3_ar]

        if external_qc is not None:
            # multi-qubit mode: don't build a private circuit, use the shared one
            self.qc = external_qc
        else:
            # single-qubit standalone mode: unchanged from your original
            self.qc = QuantumCircuit(self.b1_qr, self.b1_cr, self.b1_ar,
                                      self.b2_qr, self.b2_cr, self.b2_ar,
                                      self.b3_qr, self.b3_cr, self.b3_ar,
                                      self.phase_cr, self.phase_ar,
                                      self.final_cr)

    # Encoding, Error, BitFlip, PhaseFlip, Decoding, Sampling, Execute
    
    
    def Encoding(self):
      
        
        
    #encoding for phase flip
        self.qc.cx(self.b1_qr[0], self.b2_qr[0])
        self.qc.cx(self.b1_qr[0], self.b3_qr[0])

        #applying hadamard for phase flip
        for reg in self.qubit_regs:
            self.qc.h(reg[0])
            
        for qr in self.qubit_regs:
            self.qc.cx(qr[0], qr[1])
            self.qc.cx(qr[0], qr[2])
            
    def Error(self, x_pos = None, z_pos = None):

        # x_pos format expected from JS: {'block': 1, 'idx': 1}
        if x_pos:
            block = int(x_pos['block']) - 1  # Convert block 1,2,3 to index 0,1,2
            idx = int(x_pos['idx'])
            qr = self.qubit_regs[block]
            self.qc.x(qr[idx])

        # z_pos format expected from JS: block number (e.g., 2)
        if z_pos:
            block = int(z_pos) - 1
            # A phase flip behaves like a Z gate. 
            self.qc.z(self.qubit_regs[block][0])
            
    def BitFlip(self, qr, cr, ar):
        #syndrome measurement and parity check
        
        #parity check (if r0 != r1, cr1 = 0, if r1 != r2, cr1 = 1)
        self.qc.cx(qr[0], ar[0])
        self.qc.cx(qr[1], ar[0])
        self.qc.cx(qr[1], ar[1])
        self.qc.cx(qr[2], ar[1])
            
        self.qc.measure(ar, cr)  
        
        with self.qc.if_test((cr, 1)):
            self.qc.x(qr[2])
        with self.qc.if_test((cr, 2)):
            self.qc.x(qr[0])
        with self.qc.if_test((cr, 3)):
            self.qc.x(qr[1])
            
    def PhaseFlip(self):
        for qr in self.qubit_regs:
            for qubit in qr:
                self.qc.h(qubit)
        
        #block 1  - 3  
        for block_clbit in range(3):
            for phase_anclbit in range(2):
                self.qc.cx(self.b1_qr[block_clbit], self.phase_ar[phase_anclbit])
                
        for block_clbit in range(3):
            self.qc.cx(self.b2_qr[block_clbit], self.phase_ar[0])
            
        for block_clbit in range(3):
            self.qc.cx(self.b3_qr[block_clbit], self.phase_ar[1])
        
        self.qc.measure(self.phase_ar, self.phase_cr)   
        with self.qc.if_test((self.phase_cr, 1)):
            self.qc.x(self.b2_qr[0])
        with self.qc.if_test((self.phase_cr, 2)):
            self.qc.x(self.b3_qr[0])
        with self.qc.if_test((self.phase_cr, 3)):
            self.qc.x(self.b1_qr[0])
        
        
        for qr in self.qubit_regs:
            for qubit in qr:
                self.qc.h(qubit)
                
    def Decoding(self):
        # undo bit-flip encoding
        for qr in self.qubit_regs:
            self.qc.cx(qr[0], qr[2])
            self.qc.cx(qr[0], qr[1])

        # undo phase encoding
        for reg in self.qubit_regs:
            self.qc.h(reg[0])

        # undo block entanglement
        self.qc.cx(self.b1_qr[0], self.b3_qr[0])
        self.qc.cx(self.b1_qr[0], self.b2_qr[0])
        
    def Sampling(self):
        sampler = AerSampler()
        sampling_job = sampler.run([self.qc], shots=64)
        result = sampling_job.result()

        counts = result[0].data.final_output.get_counts()
        return counts
        
    def Execute(self, x_pos = None, z_pos = None):
        self.Encoding()
        self.Error(x_pos, z_pos)
            
        for i in range(3):
            self.BitFlip(self.qubit_regs[i], self.classical_regs[i], self.ancilla_regs[i])
            
        self.PhaseFlip()
        self.Decoding()
        
        self.qc.measure(self.qc.qubits[0], self.final_cr)
        
        #fig = (self.qc.draw(output = "mpl", fold = -1))
        #plt.show()
        
        output = self.Sampling()
                
        print(f'Probability of 1: {output.get("1", 0)/64 * 100:.2f}%')
        print(f'Probability of 0: {output.get("0", 0)/64 * 100:.2f}%')
    
    



class MultipleLogicalQubitsShorCode(ShorCode):



  '''
  __init__ -> initalises self.logical_qubits which is a list of ShorCode objects for each logical qubit in the circuit.
  all_regs combines all the registers for each indiviual logical qubit together.
  self.total_circuit is the entire quantum circuit made from all registers.

  Within __init__, each logical qubit is encoded.

  self.logical_x applies the logical X gate to the logical qubit at a specified index.

  self.logical_z applies the logical Z gate to the logical qubit at a specified index.

  self.logical_CNOT applies the logical CNOT gate to a control logical qubit and target logical qubit.
  
  
  self.error_injection_block inject X and Z errors to a specificed logical qubit using the ShorCode method for error injection.

  self.error_correction_block corrects the erros using the ShorCode error correction methods of self.BitFlip() and self.PhaseFlip()

  self.full_decoding applies the decoding procedure identically to the ShorCode class and fully decodes all the encoded logical qubits so each 
  qubit represented by the logical qubits can be sampled

  self.sampling() samples all the qubits which have been epresented as logical qubits.

  
  
  
  
  
  '''
    def __init__(self, num_logical_qubits):
        self.num_logical_qubits = num_logical_qubits
        self.logical_qubits = [ShorCode(index=i) for i in range(num_logical_qubits)]



        all_regs = []
        
        for lq in self.logical_qubits:
            all_regs += [lq.b1_qr, lq.b1_ar, lq.b1_cr,
                         lq.b2_qr, lq.b2_ar, lq.b2_cr,
                         lq.b3_qr, lq.b3_ar, lq.b3_cr,
                         lq.phase_ar, lq.phase_cr,
                         lq.final_cr]
            
        

        self.total_circuit = QuantumCircuit(*all_regs)

        self.qc = self.total_circuit  # Set the total circuit as the main circuit for this class



        for lq in self.logical_qubits:
            lq.qc = self.total_circuit



        for log_qubit in self.logical_qubits:
            log_qubit.Encoding()
        



    def logical_x(self, logical_qubit_index):

        logical_qubit = self.logical_qubits[logical_qubit_index]

        for qr in logical_qubit.qubit_regs:
            self.total_circuit.z(qr[0])



    def logical_z(self, logical_qubit_index):
        logical_qubit = self.logical_qubits[logical_qubit_index]

        for qubit in logical_qubit.b1_qr:
            logical_qubit.qc.x(qubit)

        for qubit in logical_qubit.b2_qr:
            logical_qubit.qc.x(qubit)
        
        for qubit in logical_qubit.b3_qr:
            logical_qubit.qc.x(qubit)
    
    def logical_CNOT(self, control_index, target_index):
        target_qubit = self.logical_qubits[control_index]
        control_qubit = self.logical_qubits[target_index]

        for i in range(3):
            control_qubit.qc.cx(control_qubit.b1_qr[i], target_qubit.b1_qr[i])
        for i in range(3):
            control_qubit.qc.cx(control_qubit.b2_qr[i], target_qubit.b2_qr[i])
        for i in range(3):
            control_qubit.qc.cx(control_qubit.b3_qr[i], target_qubit.b3_qr[i])
    def error_injection_block(self, logical_qubit_index, x_pos=None, z_pos=None):
        logical_qubit = self.logical_qubits[logical_qubit_index]
        logical_qubit.Error(x_pos, z_pos)

    def error_correction_block(self, logical_qubit_index):
        logical_qubit = self.logical_qubits[logical_qubit_index]
        for i in range(3):
            logical_qubit.BitFlip(logical_qubit.qubit_regs[i], logical_qubit.classical_regs[i], logical_qubit.ancilla_regs[i])
        logical_qubit.PhaseFlip()

    def full_decoding(self):
        for logical_qubit in self.logical_qubits:
            logical_qubit.Decoding()

    def sampling(self):
        for qb in self.logical_qubits:
            qb.qc.measure(qb.b1_qr[0], qb.final_cr[0])

        sampler = AerSampler()
        sampling_job = sampler.run([self.total_circuit], shots=100)
        result = sampling_job.result()  
        pub_result = result[0]
        res = []
        for i in range(self.num_logical_qubits):
            counts = getattr(pub_result.data, f"L{i}_final_output").get_counts()
            print(f"Logical qubit {i}: {counts}")
            res.append(counts)
        return res

        
   

#Example testing


multiple_logical_qubits_shor_code = MultipleLogicalQubitsShorCode(num_logical_qubits=2)

multiple_logical_qubits_shor_code.logical_x(logical_qubit_index=0)



multiple_logical_qubits_shor_code.error_injection_block(0, {'block': 1, 'idx': 1}, 0)

multiple_logical_qubits_shor_code.error_correction_block(0)

multiple_logical_qubits_shor_code.full_decoding()



res = multiple_logical_qubits_shor_code.sampling()


print("output of the total circuit:", res)




    
        
    



