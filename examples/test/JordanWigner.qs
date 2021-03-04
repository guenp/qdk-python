// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Samples.Chemistry.SimpleVQE.JordanWigner {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.ChemUtils;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Format of data passed from C# to Q# to represent terms of the Hamiltonian.
    /// Hamiltonian terms, tuple of orbital indices and coefficients
    newtype JWOptimizedHTerms = (HTerm[], HTerm[], HTerm[], HTerm[]);
    
    /// # Summary
    /// Format of data passed from C# to Q# to represent preparation of the initial state
    /// Complex number, excitation indices
    newtype JordanWignerInputState = ((Double, Double), Int[]);
    
    /// # Summary
    /// Format of data passed from C# to Q# to represent all information for Hamiltonian simulation.
    /// The meaning of the data represented is determined by the algorithm that receives it.
    newtype JordanWignerEncodingData = (Int, JWOptimizedHTerms, (Int, JordanWignerInputState[]), Double);

    function _JordanWignerClusterOperatorGeneratorIndex(data: JordanWignerInputState): GeneratorIndex {
        let ((real, imaginary), idxFermions) = data!;
        if (Length(idxFermions) == 2) {
            // PQ term
            return GeneratorIndex(([0],[real]),idxFermions);
        }
        elif(Length(idxFermions) == 4){
            // PQRS term
            return GeneratorIndex(([2],[real]),idxFermions);
        }
        else{
            // Any other term in invalid
            return GeneratorIndex(([-1],[0.0]),[0]);
        }
    }

    function _JordanWignerClusterOperatorGeneratorSystemImpl(data : JordanWignerInputState[], idx: Int) : GeneratorIndex {
        return _JordanWignerClusterOperatorGeneratorIndex(data[idx]);
    }

    function JordanWignerClusterOperatorGeneratorSystem (data : JordanWignerInputState[]) : GeneratorSystem {
        return GeneratorSystem(Length(data), _JordanWignerClusterOperatorGeneratorSystemImpl(data, _));
    }

    function _ComputeJordanWignerBitString(nFermions: Int,  idxFermions: Int[]) : Bool[] {
        if (Length(idxFermions) % 2 != 0) {
            fail $"ComputeJordanWignerString failed. `idxFermions` must contain an even number of terms.";
        }

        mutable zString = new Bool[nFermions];
        for (fermionIdx in idxFermions) {
            if (fermionIdx >= nFermions) {
                fail $"ComputeJordanWignerString failed. fermionIdx {fermionIdx} out of range.";
            }
            for (idx in 0..fermionIdx) {
                set zString w/= idx <- not zString[idx];
            }
        }

        for (fermionIdx in idxFermions){
            set zString w/= fermionIdx <- false;
        }
        return zString;
    }

    function _ComputeJordanWignerPauliZString(nFermions: Int,  idxFermions: Int[]) : Pauli[] {
        let bitString = _ComputeJordanWignerBitString(nFermions, idxFermions);
        return BoolArrayAsPauli (PauliZ, true, bitString);
    }

    function _ComputeJordanWignerPauliString(nFermions: Int,  idxFermions: Int[], pauliReplacements : Pauli[]) : Pauli[] {
        mutable pauliString = _ComputeJordanWignerPauliZString(nFermions, idxFermions);

        for (idx in IndexRange(idxFermions)) {
            let idxFermion = idxFermions[idx];
            let op = pauliReplacements[idx];
            set pauliString w/= idxFermion <- op;
        }

        return pauliString;
    }

    operation _ApplyJordanWignerClusterOperatorPQTerm (term : GeneratorIndex, stepSize: Double, qubits : Qubit[]) : Unit is Adj + Ctl {
        let ((idxTermType, coeff), idxFermions) = term!;
        let p = idxFermions[0];
        let q = idxFermions[1];
        if(p == q){
            fail($"Unitary coupled-cluster PQ failed: indices {p}, {q} must be distinct");
        }
        let angle = 0.5 * coeff[0] * stepSize;
        let ops = [[PauliX, PauliY], [PauliY, PauliX]];
        let signs = [+1.0, -1.0];

        for ((op, sign) in Zipped(ops, signs)) {
            let pauliString = _ComputeJordanWignerPauliString(Length(qubits), idxFermions, op);
            Exp(pauliString, sign * angle, qubits);
        }
    }

    function _JordanWignerClusterOperatorPQRSTermSigns(indices: Int[]) : (Int[], Double[], Double) {
        let p = indices[0];
        let q = indices[1];
        let r = indices[2];
        let s = indices[3];
        mutable sorted = new Int[4];
        mutable signs = new Double[8];
        mutable sign = 1.0;

        if(p>q){
            set sign = sign * -1.0;
        }
        if(r>s){
            set sign = sign * -1.0;
        }
        if( Min([p,q]) > Min([r,s]) ){
            set sign = sign * -1.0;
            set sorted = [Min([r,s]), Max([r,s]), Min([p,q]), Max([p,q])];
        }
        else{
            set sorted = [Min([p,q]), Max([p,q]), Min([r,s]), Max([r,s])];
        }
        // sorted is now in the order
        // [p`,q`,r`,s`], where p`<q`; r`<s`, and Min(p`,q`) is smaller than Min(r`,s`).

        let p1 = sorted[0];
        let q1 = sorted[1];
        let r1 = sorted[2];
        let s1 = sorted[3];
        // Case (p,q) < (r,s) and (p,q) > (r,s)
        if(q1 < r1){
            // p1 < q1 < r1 < s1
            return ([p1,q1,r1,s1],[1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0], sign);
        }
        // Case interleaved
        elif(q1 > r1 and q1 < s1){
            // p1 < r1 < q1 < s1
            return ([p1,r1,q1,s1],[-1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0], sign);
        }
        // Case contained
        elif(q1 > r1 and q1 > s1){
            // p1 < r1 < s1 < q1
            return ([p1,r1,s1,q1],[1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0], sign);
        }
        else{
            fail("Completely invalid cluster operator specified.");
        }
    }

    operation _ApplyJordanWignerClusterOperatorPQRSTerm (term : GeneratorIndex, stepSize: Double, qubits : Qubit[]) : Unit is Adj + Ctl {
        let ((idxTermType, coeff), idxFermions) = term!;
        let p = idxFermions[0];
        let q = idxFermions[1];
        let r = idxFermions[2];
        let s = idxFermions[3];
        let angle = 0.125 * coeff[0] * stepSize;
        
        if(p == q or p==r or p==s or q==r or q==s or r==s){
            fail($"Unitary coupled-cluster PQRS failed: indices {p}, {q}, {r}, {s} must be distinct");
        }

        let x = PauliX;
        let y = PauliY;

        let ops = [[y,y,x,y],[x,x,x,y],[x,y,y,y],[y,x,y,y],[x,y,x,x],[y,x,x,x],[y,y,y,x],[x,x,y,x]];
        let (sortedIndices, signs, globalSign) = _JordanWignerClusterOperatorPQRSTermSigns([p,q,r,s]);

        for ((op, sign) in Zipped(ops, signs)) {
            let pauliString = _ComputeJordanWignerPauliString(Length(qubits), sortedIndices, op);
            Exp(pauliString, globalSign * sign * angle, qubits);
        }
    }

    operation _JordanWignerClusterOperatorImpl(generatorIndex : GeneratorIndex, stepSize : Double, qubits : Qubit[]) : Unit is Adj + Ctl {
        let ((idxTermType, idxDoubles), idxFermions) = generatorIndex!;
        let termType = idxTermType[0];
            
        if (termType == 0) {
            _ApplyJordanWignerClusterOperatorPQTerm (generatorIndex, stepSize, qubits);
        }
        elif (termType == 2) {
            _ApplyJordanWignerClusterOperatorPQRSTerm (generatorIndex, stepSize, qubits);
        }
    }

    function _JordanWignerClusterOperatorFunction (generatorIndex : GeneratorIndex) : EvolutionUnitary {   
        return EvolutionUnitary(_JordanWignerClusterOperatorImpl(generatorIndex, _, _));
    }

    function JordanWignerClusterOperatorEvolutionSet () : EvolutionSet {
        return EvolutionSet(_JordanWignerClusterOperatorFunction(_));
    }
}
