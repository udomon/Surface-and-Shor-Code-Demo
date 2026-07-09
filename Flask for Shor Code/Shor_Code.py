import qiskit
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import SamplerV2 as AerSampler
from qiskit import transpile

class ShorCode:
    def __init__(self):
        self.b1_qr = QuantumRegister(3,'XQ1')
        self.b1_ar = QuantumRegister(2,'XA1')
        self.b1_cr = ClassicalRegister(2, 'XC1')

        self.b2_qr = QuantumRegister(3,'XQ2')
        self.b2_ar = QuantumRegister(2, 'XA2')
        self.b2_cr = ClassicalRegister(2,'XC2')


        self.b3_qr = QuantumRegister(3,'XQ3')
        self.b3_ar = QuantumRegister(2, 'XA3')
        self.b3_cr = ClassicalRegister(2,'XC3')

        self.phase_ar = QuantumRegister(2, 'ZA')
        self.phase_cr = ClassicalRegister(2,'ZC')

        self.final_cr = ClassicalRegister(1, 'final_output')

        self.qc = QuantumCircuit(self.b1_qr, self.b1_cr, self.b1_ar, 
                            self.b2_qr, self.b2_cr, self.b2_ar, 
                            self.b3_qr, self.b3_cr, self.b3_ar, 
                            self.phase_cr, self.phase_ar,
                            self.final_cr)


        self.qubit_regs = [self.b1_qr, self.b2_qr, self.b3_qr]
        self.classical_regs = [self.b1_cr, self.b2_cr, self.b3_cr]
        self.ancilla_regs = [self.b1_ar, self.b2_ar, self.b3_ar]
    
    def Encoding(self):
        qc = self.qc
        
        
    #encoding for phase flip
        qc.cx(self.b1_qr[0], self.b2_qr[0])
        qc.cx(self.b1_qr[0], self.b3_qr[0])

        #applying hadamard for phase flip
        for reg in self.qubit_regs:
            qc.h(reg[0])
            
        for qr in self.qubit_regs:
            qc.cx(qr[0], qr[1])
            qc.cx(qr[0], qr[2])
            
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
    
    
if __name__ == "__main__": 
    code = ShorCode()
    code.Execute({'block': 1, 'idx': 1}, 2)
    