from pymol import cmd

def measure_distances(nuclease_object='chain A', dna_object='chain C'):
    
    protein_selection = f'sele and {nuclease_object}'
    
    # Create a dictionary to store the minimum distances for each DNA atom
    distances = {}
    
    # Get the list of atoms in the DNA object
    dna_atoms = cmd.get_model(dna_object).atom
    for dna_atom in dna_atoms:
        min_dist = float('inf')  # Initialize minimum distance to a large value
        for protein_atom in cmd.get_model(protein_selection).atom:
            # Calculate the distance between the DNA atom and the protein atom
            dist = cmd.distance(f'dist_{dna_atom.index}_{protein_atom.index}', dna_atom.index, f'{nuclease_object} and index {protein_atom.index}')
            min_dist = min(min_dist, dist)  # Update minimum distance
        
        # Store the minimum distance for the current DNA atom
        distances[dna_atom] = min_dist
    
    print("Distance between DNA nucleotides and closest protein atom:")
    for dna_atom, dist in distances.items():
        # Get the DNA residue number and atom name
        residue_number = dna_atom.resi_number
        atom_name = dna_atom.name
        
        print(f"DNA residue {residue_number} {atom_name}: {dist:.3f} Å")

cmd.extend("measure_distances", measure_distances)
